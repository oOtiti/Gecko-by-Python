import math
from keyparameter import mxnode, mxlgr, mxlfo, mxring, mxcp, mxhyd, mxhiso
from rjtool import rjgrm
from searching import srh5
from stdgrbond import grbond
from atomtool import getatoms
from mapping import abcde_map, chemmap
from khydtool import get_hydrate
from database import nhlcdb, hlcdb_chem, hlcdb_dat
from toolbox import stoperr

# ======================================================================
# PURPOSE: This subroutine returns the effective Henry's law coefficient 
#          (Keff) for the species given as input (chem).
# 
# For details about the Gromhe SAR, see T. Raventos-Duran et al., 
# Atmospheric Chemistry and Physics, 7643-7654, 2010.
# ======================================================================
def gromhe(chem):
    import math
    from keyparameter import mxnode, mxlgr, mxlfo, mxring, mxcp, mxhyd, mxhiso
    from rjtool import rjgrm
    from searching import srh5
    from stdgrbond import grbond
    from atomtool import getatoms
    from mapping import abcde_map, chemmap
    from khydtool import get_hydrate
    from database import nhlcdb, hlcdb_chem, hlcdb_dat
    from toolbox import stoperr
    # chem tables
    tchem = chem
    if chem[:3] == '#mm':
        tchem = chem[3:]
    
    group = [""] * mxnode
    bond = [[0 for _ in range(mxnode)] for _ in range(mxnode)]
    dbflg = 0
    nring = 0
    node = 0
    rjg = [[0 for _ in range(2)] for _ in range(mxring)]  # ring-join group pairs

    # group contributors :
    kest = 0.0           # estimated intrinsic Henry's law coef.
    sigma = 0.0
    caoxa = 0.0
    caoxb = 0.0
    cato = 0
    hato = 0
    onitrofol = 0
    haloica = 0
    nogrp = 0
    noh16 = 0
    noh15 = 0

    # info returned by chemmap
    ngrp = 0
    nodetype = [""] * mxnode
    alifun = [0.0 for _ in range(21)]
    cdfun = [0.0 for _ in range(21)]
    arofun = [0.0 for _ in range(21)]
    mapfun = [[[0.0 for _ in range(21)] for _ in range(3)] for _ in range(mxnode)]
    funflg = [0 for _ in range(mxnode)]
    tabester = [[0 for _ in range(2)] for _ in range(4)]  # 1= -O- side, 2= CO side
    nfcd = 0
    nfaro = 0
    ierr = 0

    # tracks (deep=9)
    nabcde = [0 for _ in range(9)]
    tabcde = [[[0 for _ in range(mxnode)] for _ in range(mxcp)] for _ in range(9)]

    # returned by hydration (to compute Keff)
    nwa = 0
    chemhyd = [["" for _ in range(mxhiso)] for _ in range(mxhyd)]  # new formula with hydrates
    yhyd = [[0.0 for _ in range(mxhiso)] for _ in range(mxhyd)]
    nhyd = [0 for _ in range(mxhyd)]
    khydstar = 0.0

    mult = 0.0
    tsig = 0.0
    dist = 0.0
    i = 0
    j = 0
    k = 0
    l = 0
    ipos = 0
    it = 0
    i1 = 0
    i2 = 0
    ncd = 0
    inte = 0
    nato = 0
    oato = 0
    ir = 0
    is_ = 0  # 避免与Python关键字冲突
    flato = 0
    brato = 0
    clato = 0

    progname = 'gromhe '
    mesg = ""

    # Tafta sigma values, alifatic compounds:
    #  1= -OH    ;  2= -NO2     ;  3= -ONO2   ; 4= -OOH    ;  5= -F     
    #  6= -Cl    ;  7=  -Br     ;  8= -I      ; 9= -CHO    ; 10= -CO-   
    # 11= -COOH  ; 12= -CO(OOH) ; 13= -PAN    ; 14= -O-    ; 15= R-COO-R
    # 16 = HCO-O-R; 17= -CO(F) ; 18= -CO(Cl) ;  19= -CO(Br) ; 20= -CO(I)
    # 21 = -CO(O-)
    # Tafta sigma values for aliphatic compounds. sigma for ester depends
    # whether ester is connected from the -CO- side or the -O- side. 
    # Sigma=2.56 for -O- side and 2.00 for CO side.
    # CMV 16/06/14 : add sigma taft for CO(O-), n°21
    alisig = [
        0.62, 1.47, 1.38, 0.62, 1.10, 0.94, 1.00, 1.00, 2.15, 1.81,
        2.08, 2.08, 2.00, 1.81, 2.56, 2.90, 2.44, 2.37, 2.37, 2.37,
        -1.06
    ]
    sigester = 2.00  # sigma for ester, CO side

    # Check if henry's law constant is already known in the database
    #-------------------------------------------------------------
    ipos = srh5(tchem, hlcdb_chem, nhlcdb)
    if ipos > 0:
        Keff = math.log10(hlcdb_dat[ipos-1][0])
        return Keff

    # build the group and bond matrix for chem
    grbond(tchem, group, bond, dbflg, nring)
    for i in range(1, mxnode + 1):
        bond[i-1][i-1] = 0
    rjgrm(nring, group, rjg)  # rm ring index and get ring closing nodes

    node = 0
    for i in range(1, mxnode + 1):
        if group[i-1][:1] != ' ':
            node += 1

    # Get hydration constant
    # ----------------------
    get_hydrate(tchem, nwa, nhyd, chemhyd, yhyd, khydstar)

    # ---------------------------------
    # Count functionalities
    # ---------------------------------

    # get the number of atom - C and H counters
    getatoms(tchem, cato, hato, nato, oato, ir, is_, flato, brato, clato)
        
    ierr = 0
    chemmap(tchem, node, group, bond, ngrp, nodetype,
            alifun, cdfun, arofun, mapfun, funflg,
            tabester, nfcd, nfaro, ierr)
    if ierr != 0:
        mesg = "chemmap returned problems"
        stoperr(progname, mesg, chem)

    # ---------------------------------
    # compute taft sigmas
    # ---------------------------------
    sigma = 0.0

    # 1-1 interaction
    # ----------------
    for i in range(1, node + 1):
        if funflg[i-1] > 1:
            for k in range(1, 9):
                for it in range(1, 3):  # check for both saturated and unsaturated aliphatic
                    # self interaction
                    if mapfun[i-1][it-1][k-1] != 0.0:
                        if mapfun[i-1][it-1][k-1] > 1.0:
                            sigma += (alisig[k-1] / 0.4) * mapfun[i-1][it-1][k-1] * (mapfun[i-1][it-1][k-1] - 1.0)
                    # cross interaction
                    for l in range(1, 9):
                        if (mapfun[i-1][it-1][k-1] > 0.0) and (l != k):
                            sigma += (alisig[k-1] / 0.4) * mapfun[i-1][it-1][k-1] * mapfun[i-1][it-1][l-1]

    # multiple node interaction
    # -------------------------
    for i in range(1, node + 1):    # start loop among node
        if funflg[i-1] != 0:   # manage node i
            abcde_map(bond, i, node, nabcde, tabcde)

            # ester functionality is on 2 nodes - remove the path that would lead
            # to count some interaction twice. This is performed by "killing" the
            # "second" node of the functionality.
            for k in range(1, 5):
                if tabester[k-1][0] != 0:
                    for ipos in range(2, 10):
                        for it in range(1, nabcde[ipos-1] + 1):
                            for j in range(1, ipos):
                                if tabcde[ipos-1][it-1][j-1] == tabester[k-1][0]:
                                    if tabcde[ipos-1][it-1][j] == tabester[k-1][1]:
                                        tabcde[ipos-1][it-1][j] = 0  # set the node to 0
                                if tabcde[ipos-1][it-1][j-1] == tabester[k-1][1]:
                                    if tabcde[ipos-1][it-1][j] == tabester[k-1][0]:
                                        tabcde[ipos-1][it-1][j] = 0  # set the node to 0

            # start loop for sigma effect
            for j in range(1, 21):
                if (mapfun[i-1][0][j-1] != 0.0) or (mapfun[i-1][1][j-1] != 0.0):  # start sigma - taft
                    for inte in range(2, 10):
                        for k in range(1, nabcde[inte-1] + 1):
                            ncd = 0
                            if tabcde[inte-1][k-1][inte-1] == 0:
                                continue
                            if funflg[tabcde[inte-1][k-1][inte-1]-1] == 0:
                                continue

                            # CMV! rule out case when tabcde(int,k,int) == 0 because
                            # CMV! associated funflg(0) (invisible seg fault!)
                            # CMV! can be /= 0, depending on memory issues

                            for l in range(1, inte):  # start the loop to see if Cd are in between
                                i1 = tabcde[inte-1][k-1][l-1]
                                if i1 == 0:
                                    continue
                                i2 = tabcde[inte-1][k-1][l]
                                if i2 == 0:
                                    continue
                                if bond[i1-1][i2-1] == 2:
                                    ncd += 1
                            mult = float(funflg[tabcde[inte-1][k-1][inte-1]-1])
                            if tabcde[inte-1][k-1][1] == 0:
                                mult = 0.0  # kill second node on ester
                            tsig = alisig[j-1]
                            if j == 15:  # check ester
                                if group[i-1][:2] == 'CO':
                                    tsig = sigester  # CO instead of -O-
                            dist = float(inte) - float(ncd) - 2.0
                            if group[tabcde[inte-1][k-1][inte-1]-1][:3] == '-O-':
                                dist -= 1.0  # decrease the distance
                            if nodetype[i-1] == 'o':
                                dist -= 1.0
                            sigma += (tsig * (0.4 ** dist) * mult)

    # lump group on Cd and saturated aliphatic
    for i in range(1, 21):
        if cdfun[i-1] > 0.0:
            alifun[i-1] += cdfun[i-1]

    # lump group on aromatic and saturated aliphatic
    for i in range(1, 21):
        if arofun[i-1] > 0.0:
            alifun[i-1] += arofun[i-1]

    # -------------------------------------
    # compute other group-group descriptor
    # -------------------------------------

    # compute correcting factors
    nogrp = 1
    if nogrp != 0:
        nogrp = 0

    # check for ortho nitro phenols
    onitrofol = nitrofol(node, group, bond)

    # check for halogen next to a carboxylic group
    haloica = haloic(node, group, bond)

    # check for H bounding from an alcohol
    noh15 = 0
    noh16 = 0
    if alifun[0] > 0.0:
        noh15, noh16 = hydol(node, group, bond)

    # check for -CO-C(X) and -CO-C-C(X) structure
    caoxa, caoxb = caox(node, group, bond)

    # ALL DATA
    kest = -1.51807             \
           - 0.13764 * sigma       \
           - 2.66308 * onitrofol   \
           + 0.973418 * haloica    \
           - 1.76731 * caoxa       \
           - 1.09538 * caoxb       \
           - 1.02636 * noh16       \
           - 0.596601 * noh15      \
           + 0.499757 * cato       \
           - 0.31011 * hato        \
           - 0.27637 * nogrp       \
           - 0.524482 * nfcd       \
           - 1.12707 * nfaro       \
           + 4.56593 * alifun[0]   \
           + 3.01781 * alifun[1]   \
           + 2.03719 * alifun[2]   \
           + 4.8687 * alifun[3]    \
           + 0.59747 * alifun[4]   \
           + 0.870225 * alifun[5]  \
           + 1.05894 * alifun[6]   \
           + 1.22831 * alifun[7]   \
           + 2.59268 * alifun[8]   \
           + 3.16535 * alifun[9]   \
           + 5.09321 * alifun[10]  \
           + 4.683 * alifun[11]    \
           + 1.93179 * alifun[12]  \
           + 2.4008 * alifun[13]   \
           + 2.78457 * alifun[14]  \
           + 2.35694 * alifun[15]  

    # keff=10**kest
    # keff=keff*(1+khydstar)
    # keff=log10(keff)
    # Keff=kest! + log10(1+khydstar)
    Keff = kest + math.log10(1 + khydstar)
    return Keff


#=======================================================================
# Purpose: Count the number of ortho nitro-phenols group in a molecule
# ======================================================================
def nitrofol(ng, group, bond):
    onitrofol = 0
    tempoh16 = 0
    i = 0
    k = 0

    # only hydroxy group (whether alkohol or carboxylic) is seek of H bonding
    for i in range(1, ng + 1):
        if 'OH' in group[i-1]:  # OPEN 'OH'
            if 'CO(OH)' in group[i-1]:
                continue  # EXCLUDE ACID H
            if group[i-1][:1] != 'c':
                continue  # EXCLUDE non aromatic OH

            # search for nitro group on the alpha carbon
            tempoh16 = 0     # a phenol is found
            for k in range(1, ng + 1):
                if group[k-1][:1] != 'c':
                    continue  # EXCLUDE non aromatic groups
                if bond[i-1][k-1] != 0:
                    if 'NO2' in group[k-1]:
                        tempoh16 += 1

            # increment is maximum 1 per OH (a given OH can only be involded in a single bond)
            if tempoh16 > 0:
                onitrofol += 1

    return onitrofol


#=======================================================================
# Purpose: Count the # of -C(X)-CO(OH) structure, where X is an halogen 
#=======================================================================
def haloic(ng, group, bond):
    haloica = 0
    tempoh15 = 0
    i = 0
    k = 0

    # start loop - search for carboxylic group
    for i in range(1, ng + 1):
        if 'CO(OH)' in group[i-1]:
            tempoh15 = 0
            
            # H-bond with a functional group on the alpha carbon
            for k in range(1, ng + 1):
                if bond[i-1][k-1] != 0:
                    if '(F)' in group[k-1]:
                        tempoh15 += 1
                    if '(Cl)' in group[k-1]:
                        tempoh15 += 1
                    if '(Br)' in group[k-1]:
                        tempoh15 += 1
                    if '(I)' in group[k-1]:
                        tempoh15 += 1

            # increment is maximum 1 per CO(OH) (a given OH can only be involded in a
            # single bond)
            if tempoh15 > 0:
                haloica += 1

    return haloica


#=======================================================================
# Purpose: count the # of -CO-C(X) -and CO-C-C(X) structures, where X is
# an oxygenated group (carbonyl,alkohol, nitro, hydroperoxide, ether) 
#=======================================================================
def caox(ng, group, bond):
    caoxa = 0.0  # caoxa : # of -CO-C(X)
    caoxb = 0.0  # caoxb : # of -CO-C-C(X)
    tempoh15 = 0
    tempoh16 = 0
    i = 0
    k = 0
    l = 0
    m = 0
    ncoco15 = 0.0
    ncoco16 = 0.0

    # start loop - search for carbonyl
    for i in range(1, ng + 1):
        if (group[i-1][:2] == 'CO') or (group[i-1][:3] == 'CHO'):
            tempoh15 = 0
            tempoh16 = 0

            # look for group on alpha position 
            for k in range(1, ng + 1):
                if bond[i-1][k-1] != 0:
                    if '(NO2)' in group[k-1]:
                        tempoh15 += 1
                    if group[k-1][:2] == 'CO':
                        tempoh15 += 1
                        ncoco15 += 1
                    if '(OH)' in group[k-1]:
                        tempoh15 += 1
                    if 'CHO' in group[k-1]:
                        tempoh15 += 1
                        ncoco15 += 1
                    if '(OOH)' in group[k-1]:
                        tempoh15 += 1

                    # look for group on beta position 
                    for l in range(1, ng + 1):
                        if (group[i-1][:1] == 'c') and (group[l-1][:1] == 'c'):
                            continue
                        if (bond[k-1][l-1] != 0) and (l != i):
                            if '(NO2)' in group[l-1]:
                                tempoh16 += 1
                            if group[l-1][:2] == 'CO':
                                tempoh16 += 1
                                ncoco16 += 1
                            if group[l-1][:3] == '-O-':
                                tempoh15 += 1
                            if '(OH)' in group[l-1]:
                                tempoh16 += 1
                            if '(OOH)' in group[l-1]:
                                tempoh16 += 1
                            if 'CHO' in group[l-1]:
                                tempoh16 += 1
                                ncoco16 += 1

                            # check for ether on position gamma
                            for m in range(1, ng + 1):
                                if (bond[l-1][m-1] != 0) and (m != k):
                                    if group[m-1][:3] == '-O-':
                                        tempoh16 += 1

            if tempoh15 > 0:
                caoxa += 1.0
            if tempoh16 > 0:
                caoxb += 1.0

    # correct the value for dicarbonyl (otherwise counted twice)
    if ncoco15 > 1:
        ncoco15 /= 2.0
        caoxa -= ncoco15
    if ncoco16 > 1:
        ncoco16 /= 2.0
        caoxb -= ncoco16

    return caoxa, caoxb


#=======================================================================
# This subroutine count the number of alcohol moiety leading to 
# intramolecular H bounding thrue a 6 or 5 member ring. Priority
# is given to 6 member ring. A alcohol moiety can only be counted once.
#=======================================================================
def hydol(ng, group, bond):
    from keyparameter import mxcp
    from mapping import abcde_map
    from ringtool import findring

    noh15 = 0
    noh16 = 0
    nohoh15 = 0
    nohoh16 = 0
    tempoh15 = 0
    tempoh16 = 0
    i = 0
    j = 0
    k = 0
    l = 0
    i1 = 0
    i2 = 0
    tohoh15 = 0
    tohoh16 = 0

    ndeep = 4
    nabcde = [0 for _ in range(ndeep)]
    tabcde = [[[0 for _ in range(len(group))] for _ in range(mxcp)] for _ in range(ndeep)]

    rngflg = 0       # 0 = 'no', 1 = 'yes'
    ring = [0 for _ in range(len(group))]    # =1 if node participates in current ring
    tring = 0

    # start loop - only hydroxy group (whether alkohol or carboxylic) is seek of H bonding
    for i in range(1, ng + 1):
        if '(OH)' in group[i-1]:              # OPEN 'OH'
            if 'CO(OH)' in group[i-1]:
                continue  # EXCLUDE ACID H 

            tohoh16 = 0
            tohoh15 = 0
            tempoh16 = 0
            tempoh15 = 0
            abcde_map(bond, i, ng, nabcde, tabcde)

            # gamma position - H bonding occurs only if the nodes does not belong to a cycle
            for k in range(1, nabcde[3] + 1):
                l = tabcde[3][k-1][3]
                if group[l-1][:3] == '-O-':
                    i1 = i
                    i2 = tabcde[3][k-1][1]
                    findring(i1, i2, ng, bond, rngflg, ring)
                    tring = 0
                    if rngflg != 0:
                        for j in range(1, 5):
                            if ring[tabcde[3][k-1][j-1]-1] != 0:
                                tring += 1
                    if tring <= 2:
                        tempoh16 += 1

            # beta position
            for k in range(1, nabcde[2] + 1):
                l = tabcde[2][k-1][2]
                i1 = i
                i2 = tabcde[2][k-1][1]
                findring(i1, i2, ng, bond, rngflg, ring)
                tring = 0
                if rngflg != 0:
                    for j in range(1, 4):
                        if ring[tabcde[2][k-1][j-1]-1] != 0:
                            tring += 1

                if tring <= 2:  # exit if nodes belong to the same ring
                    if '(ONO2)' in group[l-1]:
                        tempoh16 += 1
                    if '(F)' in group[l-1]:
                        tempoh16 += 1
                    if '(Cl)' in group[l-1]:
                        tempoh16 += 1
                    if '(Br)' in group[l-1]:
                        tempoh16 += 1
                    if '(I)' in group[l-1]:
                        tempoh16 += 1
                    if '(OOH)' in group[l-1]:
                        tempoh16 += 1
                    if 'CHO' in group[l-1]:
                        tempoh16 += 1
                    if group[l-1][:3] == 'CO ':
                        tempoh16 += 1
                    if '(OH)' in group[l-1]:
                        tempoh16 += 1
                        tohoh16 += 1       # 2 hydroxy can only make 1 bond 
                    if group[l-1][:3] == '-O-':
                        tempoh15 += 1

            # alpha position
            for k in range(1, nabcde[1] + 1):
                l = tabcde[1][k-1][1]
                if '(NO2)' in group[l-1]:
                    tempoh16 += 1
                if '(ONO2)' in group[l-1]:
                    tempoh15 += 1
                if '(F)' in group[l-1]:
                    tempoh15 += 1
                if '(Cl)' in group[l-1]:
                    tempoh15 += 1
                if '(Br)' in group[l-1]:
                    tempoh15 += 1
                if '(I)' in group[l-1]:
                    tempoh15 += 1
                if '(OOH)' in group[l-1]:
                    tempoh15 += 1
                if 'CHO' in group[l-1]:
                    tempoh15 += 1
                if group[l-1][:3] == 'CO ':
                    tempoh15 += 1
                if '(OH)' in group[l-1]:
                    tempoh15 += 1
                    tohoh15 += 1

            # curent position
            if '(OOH)' in group[i-1]:
                tempoh15 += 1

            # Analyse H bonds. A given -OH make only one bond. Priority is given to 
            # 6 member ring, then 5 member ring. Care must be taken for dihydroxy 
            # species, which can only make one H-bond
            if tempoh16 > 0:
                noh16 += 1
                if (tohoh16 > 0) and (tempoh16 == 1):
                    nohoh16 += 1   # Hbond through dihydroxy group only
            elif tempoh15 > 0:
                noh15 += 1
                if (tohoh15 > 0) and (tempoh15 == 1):
                    nohoh15 += 1   # Hbond through dihydroxy group only

    nohoh16 = nohoh16 // 2       # integer division on purpose
    if nohoh16 != 0:
        noh16 -= nohoh16
    nohoh15 = nohoh15 // 2       # integer division on purpose
    if nohoh15 != 0:
        noh15 -= nohoh15

    return noh15, noh16