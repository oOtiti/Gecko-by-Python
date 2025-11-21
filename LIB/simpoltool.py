# ======================================================================
# Purpose : Compute the vapor pressures and enthalpies of vaporization
# of the species provided as input based on the simpol SAR.
# Ref. for the SAR: Pankow, J. F. and Asher, W. E.: SIMPOL.1: a simple 
# group contribution method for predicting vapor pressures and 
# enthalpies of vaporization of multifunctional organic compounds, 
# Atmos. Chem. Phys., 2773-2796, 2008. 
# https://doi.org/10.5194/acp-8-2773-2008
# https://www.atmos-chem-phys.net/8/2773/2008/
# ======================================================================
from database import psimpgrp298, hsimpgrp298, p0simpgrp298, h0simpgrp298

def simpol(chem, bond, group, nring):
    temp = 298.0  # ref temperature for Psat and latent heat
    simpgroup = [0] * 30  # simpgroup(30) in Fortran
    
    # get simpol groups
    simpgroup = simpolgrp(chem, group, bond, nring)
    
    # add constant term
    logPvap = p0simpgrp298
    latentheat = h0simpgrp298
    
    # add group contribution
    for j in range(30):  # Fortran: j=1~30
        if simpgroup[j] != 0:
            logPvap += simpgroup[j] * psimpgrp298[j]
            latentheat += simpgroup[j] * hsimpgrp298[j]
    
    latentheat = -2.303 * 8.314 * latentheat
    return logPvap, latentheat

# ======================================================================
# Purpose: read and store the values of the groups for the Simpol SAR 
# ======================================================================
from keyparameter import tfu1, dirgecko
import math

def load_simpoldat():
    temp = 298.0  # ref temperature for Psat and latent heat
    bk0 = [0.0] * 4  # parameter of the simpol SAR ("constant" group)
    bkgr = [[0.0 for _ in range(4)] for _ in range(30)]  # parameter of the simpol SAR
    simpdat = [[0.0 for _ in range(4)] for _ in range(31)]  # SIZE(bkgr,1)+1=31, parameter of the simpol SAR
    
    filename = dirgecko.strip() + 'DATA/simpol.dat'
    try:
        with open(filename, 'r') as f:
            ilin = 0
            while True:
                line = f.readline()
                if not line:
                    break
                if line[0] == '!':
                    continue
                if line[:3] == 'END':
                    break
                
                ilin += 1
                if ilin > 31:
                    break
                # read line into simpdat
                parts = list(map(float, line.strip().split()))
                for j in range(4):
                    simpdat[ilin-1][j] = parts[j]  # Fortran: simpdat(ilin,j), Python: 0-based
    except IOError as ierr:
        print(f'--error--, in load_simpoldat while trying to open file:{filename.strip()}')
        raise SystemExit("in load_simpoldat") from ierr
    
    if ilin != 31:
        print(f'--error--, unexpected # of data in :{filename.strip()}')
        raise SystemExit("in load_simpoldat")
    
    # assign constant term and group terms
    bk0[:] = simpdat[0]  # Fortran: bk0(:)=simpdat(1,:)
    for j in range(30):
        bkgr[j][:] = simpdat[j+1]  # Fortran: bkgr(1:,:)=simpdat(2:,:)
    
    # store the overall group contribution @ the reference T (here 298K)
    global p0simpgrp298, h0simpgrp298, psimpgrp298, hsimpgrp298
    p0simpgrp298 = (bk0[0]/temp) + bk0[1] + (bk0[2]*temp) + bk0[3]*math.log(temp)
    h0simpgrp298 = bk0[0] - bk0[2]*temp*temp - bk0[3]*temp
    
    for j in range(30):
        psimpgrp298[j] = (bkgr[j][0]/temp) + bkgr[j][1] + (bkgr[j][2]*temp) + bkgr[j][3]*math.log(temp)
        hsimpgrp298[j] = bkgr[j][0] - bkgr[j][2]*temp*temp - bkgr[j][3]*temp

# ======================================================================
# Purpose: provide the simpol groups for the simpol SAR to estimate
# the vapor pressure of organic species.
# ======================================================================
from keyparameter import mxtrk, mxlcd, mxlest
from ringtool import findring
from mapping import estertrack
from toolbox import stoperr, countstring

def simpolgrp(chem, group, bond, nring):
    simpgroup = [0] * 30  # output: simpgroup(30) in Fortran
    n_bond = len(bond)
    # internal variables
    neigh = [['' for _ in range(4)] for _ in range(n_bond)]  # neigh(SIZE(bond,1),4)
    tgroup = group.copy()  # tgroup(SIZE(group))
    ring = [0] * n_bond  # ring(SIZE(bond,1)): 0=no ring, 1=yes ring
    rngflg = 0  # 0 = 'no ring', 1 = 'yes ring'
    ngr = 0  # number of nodes
    ogr = 0  # number of -O- nodes
    ncd = 0  # number of Cd node
    nbnei = [0] * n_bond  # nbnei(SIZE(bond,1))
    nei_ind = [[0 for _ in range(4)] for _ in range(n_bond)]  # nei_ind(SIZE(bond,1),4)
    netrack = 0  # # of ester tracks
    etrack = [[0 for _ in range(mxlest)] for _ in range(mxtrk)]  # etrack(mxtrk,mxlest)
    etracklen = [0] * mxtrk  # etracklen(mxtrk)
    progname = 'simpolgrp'
    
    # initialize
    for i in range(n_bond):
        ring[i] = 0
        nbnei[i] = 0
        for j in range(4):
            neigh[i][j] = ''
            nei_ind[i][j] = 0
    for i in range(30):
        simpgroup[i] = 0
    
    # count the number of nodes
    ngr = sum(1 for g in group if g.strip() != '')
    # count number of -O- nodes
    ogr = sum(1 for g in group[:ngr] if g.strip() == '-O-')
    # count number of Cd nodes
    ncd = sum(1 for g in group[:ngr] if g.strip()[:2] == 'Cd')
    
    # if rings exist, find nodes belonging to rings
    if nring > 0:
        # call findring(1,2,ngr,bond,rngflg,ring) - Fortran 1-based parameters
        findring(1, 2, ngr, bond, rngflg, ring)
    
    # find neighbours
    for i in range(ngr):  # Fortran: i=1~ngr (1-based), Python: 0~ngr-1
        for j in range(ngr):  # Fortran: j=1~ngr
            if bond[i][j] != 0:  # bond(i,j) in Fortran is bond[i-1][j-1] in Python? No: bond is input as Fortran's bond(:,:), so use directly
                nbnei[i] += 1  # number of neighbours
                neigh[i][nbnei[i]-1] = group[j]  # Fortran: neigh(i,nbnei(i))=group(j), Python: 0-based
                nei_ind[i][nbnei[i]-1] = j + 1  # store Fortran 1-based index
    
    # carbon number (GROUP 1) 
    # -----------------------
    simpgroup[0] = ngr - ogr  # Fortran: simpgroup(1)
    
    # Rings (GROUP 3 & 4) - need revision but should be fine if aromatic has only 1 ring
    # --------------------
    if nring != 0:
        if 'c1' in chem:
            simpgroup[2] = nring  # GROUP 3 (Fortran: simpgroup(3))
        else:
            simpgroup[3] = nring  # GROUP 4 (Fortran: simpgroup(4))
    
    # C=C bonds (GROUP 5 >C=C< and 6 C=C-C=O in a ring)
    # --------------------
    if ncd != 0:
        simpgroup[4] = ncd // 2  # GROUP 5 (Fortran: simpgroup(5))
        for i in range(ngr):  # Fortran: i=1~ngr
            if group[i].strip()[:2] == 'Cd':
                ico = 0
                icd = 0
                for j in range(nbnei[i]):  # Fortran: j=1~nbnei(i)
                    # check neighbour groups
                    if neigh[i][j].strip()[:2] == 'CO':
                        ico = j + 1  # store Fortran 1-based j
                    if neigh[i][j].strip()[:2] == 'Cd':
                        icd = j + 1  # store Fortran 1-based j
                if ico != 0 and icd != 0:
                    # check if all nodes are in ring (ring is 0-based, nei_ind stores 1-based)
                    node_i_ring = ring[i] == 1
                    node_ico_ring = ring[nei_ind[i][ico-1]-1] == 1  # nei_ind[i][ico-1] is 1-based
                    node_icd_ring = ring[nei_ind[i][icd-1]-1] == 1
                    if node_i_ring and node_ico_ring and node_icd_ring:
                        simpgroup[5] += 1  # GROUP 6 (Fortran: simpgroup(6))
    
    # Ester (GROUP 11 & 30) - formate considered as ester
    # --------------------
    if ogr != 0:
        # make ester tracks, i.e. (CO-O)x tracks
        estertrack(chem, bond, group, ngr, netrack, etracklen, etrack)
        
        for i in range(netrack):  # Fortran: i=1~netrack
            track_idx = i  # Fortran i -> Python track_idx (0-based)
            et_len = etracklen[track_idx]
            nester = et_len // 2  # integer division
            ogr -= nester  # keep record of remaining -O-
            
            # "simple" CO-O, CO-O-CO-O, etc structures (even # of nodes in track)
            if et_len % 2 == 0:
                # check for nitroester >C(NO2)-CO-O- special case
                nitroester = 0
                # find inod (Fortran 1-based)
                if tgroup[etrack[track_idx][0]-1].strip()[:2] == 'CO':  # etrack(i,1) in Fortran
                    inod = etrack[track_idx][0]
                else:
                    inod = etrack[track_idx][et_len-1]  # etrack(i,etracklen(i)) in Fortran
                
                # check neighbours of inod
                for j in range(nbnei[inod-1]):  # inod is 1-based, nbnei[inod-1] is 0-based count
                    if '(NO2)' in neigh[inod-1][j]:
                        nitroester = 1
                
                if nitroester == 0:
                    simpgroup[10] += nester  # GROUP 11 (Fortran: simpgroup(11))
                else:
                    simpgroup[29] += 1  # GROUP 30 (Fortran: simpgroup(30))
                    simpgroup[10] += nester - 1  # GROUP 11
                
                # "erase" all groups in the track
                for j in range(et_len):  # Fortran: j=1~etracklen(i)
                    node = etrack[track_idx][j]  # 1-based node index
                    tgroup[node-1] = ' '  # tgroup(node) = ' ' in Fortran
            
            # "Incomplete" ester structure O-CO-O, CO-O-CO etc (odd # of nodes in track)
            else:
                simpgroup[10] += nester  # GROUP 11 (ignore possible nitroester)
                ifirst = etrack[track_idx][0]  # etrack(i,1) in Fortran (1-based)
                
                # -O-CO-O- structure  - carbonate
                if tgroup[ifirst-1].strip() == '-O-':  # tgroup(ifirst) in Fortran
                    ipero = 0
                    # check neighbours of ifirst (1-based)
                    for j in range(nbnei[ifirst-1]):  # j=1~nbnei(ifirst) in Fortran
                        if neigh[ifirst-1][j].strip()[:3] == '-O-':
                            ipero = 1
                    if ipero == 0:
                        # kill from 1st, leave last
                        for j in range(et_len - 1):  # j=1~etracklen(i)-1 in Fortran
                            node = etrack[track_idx][j]
                            tgroup[node-1] = ' '
                    else:
                        # leave first (peroxide), kill others
                        for j in range(1, et_len):  # j=2~etracklen(i) in Fortran
                            node = etrack[track_idx][j]
                            tgroup[node-1] = ' '
                
                # -CO-O-CO- structure - Anhydre
                else:
                    if tgroup[ifirst-1].strip()[:3] == 'CHO':  # keep aldehyde, kill others
                        for j in range(1, et_len):  # j=2~etracklen(i) in Fortran
                            node = etrack[track_idx][j]
                            tgroup[node-1] = ' '
                    else:
                        # keep last, kill others
                        for j in range(et_len - 1):  # j=1~etracklen(i)-1 in Fortran
                            node = etrack[track_idx][j]
                            tgroup[node-1] = ' '
    
    # Peroxide R-O-O-R (GROUP 26)
    # ---------------------------
    if ogr != 0:
        pero_i = 0
        while pero_i < ngr:  # Fortran: i=1~ngr
            i = pero_i
            if tgroup[i].strip()[:3] == '-O-':  # tgroup(i) in Fortran
                for j in range(nbnei[i]):  # j=1~nbnei(i) in Fortran
                    inod = nei_ind[i][j]  # 1-based node index
                    if tgroup[inod-1].strip()[:3] == '-O-':  # tgroup(inod) in Fortran
                        tgroup[i] = ' '
                        tgroup[inod-1] = ' '
                        ogr -= 2
                        simpgroup[25] += 1  # GROUP 26 (Fortran: simpgroup(26))
                        pero_i += 1
                        continue  # CYCLE pero
            pero_i += 1
    
    # Ether (GROUP 12, 13, 14)
    # ---------------------------
    if ogr != 0:
        ether_i = 0
        while ether_i < ngr:  # Fortran: i=1~ngr
            i = ether_i
            if tgroup[i].strip()[:3] == '-O-':  # tgroup(i) in Fortran
                if ring[i] == 0:  # ring(i) in Fortran
                    simpgroup[11] += 1  # GROUP 12 (Fortran: simpgroup(12))
                    tgroup[i] = ' '
                    ogr -= 1
                else:
                    # check if both neighbours are 'c' (aromatic)
                    neigh1 = neigh[i][0].strip()[:1] if nbnei[i] >=1 else ''
                    neigh2 = neigh[i][1].strip()[:1] if nbnei[i] >=2 else ''
                    if neigh1 == 'c' and neigh2 == 'c':
                        simpgroup[13] += 1  # GROUP 14 (Fortran: simpgroup(14))
                    else:
                        simpgroup[12] += 1  # GROUP 13 (Fortran: simpgroup(13))
                    tgroup[i] = ' '
                    ogr -= 1
            ether_i += 1
    
    # check that all -O- group were considered
    if ogr != 0:
        mesg = "expect no more -O- groups remains but did not happen"
        stoperr(progname, mesg, chem)
    
    # loop over ending - unique nodes
    # -------------------------------
    for i in range(ngr):  # Fortran: i=1~ngr
        if tgroup[i].strip() == '':
            continue
        if tgroup[i].strip() == 'CHO':
            simpgroup[7] += 1  # GROUP 8 (Fortran: simpgroup(8))
            tgroup[i] = ' '
        elif tgroup[i].strip() == 'CO ':
            simpgroup[8] += 1  # GROUP 9 (Fortran: simpgroup(9))
            tgroup[i] = ' '
        elif tgroup[i].strip() == 'CO(OH)':
            simpgroup[9] += 1  # GROUP 10 (Fortran: simpgroup(10))
            tgroup[i] = ' '
        elif tgroup[i].strip() == 'CO(OOH)':
            simpgroup[27] += 1  # GROUP 28 (Fortran: simpgroup(28))
            tgroup[i] = ' '
        elif tgroup[i].strip() == 'CO(OONO2)':
            simpgroup[24] += 1  # GROUP 25 (Fortran: simpgroup(25))
            tgroup[i] = ' '
        elif tgroup[i].strip() == 'CO(OOOH)':
            simpgroup[27] += 1  # GROUP 28 (Fortran: simpgroup(28)) - peracid like
            tgroup[i] = ' '
        elif tgroup[i].strip() == 'CO(ONO2)':
            simpgroup[24] += 1  # GROUP 25 (Fortran: simpgroup(25)) - PAN like
            tgroup[i] = ' '
        elif tgroup[i].strip() == 'CO(NO2)':
            simpgroup[8] += 1  # GROUP 9 (Fortran: simpgroup(9))
            simpgroup[15] += 1  # GROUP 16 (Fortran: simpgroup(16))
            tgroup[i] = ' '
    
    # loop over other moieties
    # ------------------------
    aronitro = 0
    for i in range(ngr):  # Fortran: i=1~ngr
        if tgroup[i].strip() == '':
            continue
        
        # count (OH) groups
        nn = countstring(tgroup[i], '(OH)')
        if nn != 0:
            if tgroup[i][0] == 'C':
                simpgroup[6] += nn  # GROUP 7 (Fortran: simpgroup(7))
            elif tgroup[i][0] == 'c':
                simpgroup[16] += nn  # GROUP 17 (Fortran: simpgroup(17))
        
        # count (OOH) groups
        nn = countstring(tgroup[i], '(OOH)')
        if nn != 0:
            simpgroup[26] += nn  # GROUP 27 (Fortran: simpgroup(27))
        
        # count (OOOH) groups
        nn = countstring(tgroup[i], '(OOOH)')
        if nn != 0:
            simpgroup[26] += nn  # GROUP 27 (Fortran: simpgroup(27)), OOH like
        
        # count (ONO2) groups
        nn = countstring(tgroup[i], '(ONO2)')
        if nn != 0:
            simpgroup[14] += nn  # GROUP 15 (Fortran: simpgroup(15))
        
        # count (NO2) groups
        nn = countstring(tgroup[i], '(NO2)')
        if nn != 0:
            simpgroup[15] += nn  # GROUP 16 (Fortran: simpgroup(16))
            if tgroup[i][0] == 'c':
                aronitro = 1  # raise flag for nitrophenols, GROUP 30
    
    # check nitrophenol - switch phenol (17) to nitrophenol (29)
    if simpgroup[15] != 0 and simpgroup[16] != 0:
        simpgroup[28] = simpgroup[16]  # GROUP 29 (Fortran: simpgroup(29))
        simpgroup[16] = 0  # reset GROUP 17 (Fortran: simpgroup(17))
    
    return simpgroup