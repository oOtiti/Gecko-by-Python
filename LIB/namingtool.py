from keyparameter import mxnode
from database import  nfn,namfn,chemfn
from atomtool import cnum
from searching import srh5
from dictstackdb import nrec,namlst
from toolbox import  nfn,namfn,chemfn
from atomtool import getatoms
from toolbox import stoperr
from rjtool import rjsrm
def naming(chem, iloc, name, fgrp):
    nord = 20
    ord_str = '1234XTEPAGHNUDVKORCZ'
    nalfa = 61
    alfa = '0123456789ABCDEFGHIJKLMNOPQRSTUWXYZabcdefghijklmnopqrstuvwxyz'
    
    ii = 0
    jj = 0
    kk = 0

    ic = cnum(chem)
    if ic < 2 or ic > mxnode:
        mesg = "Illegal species:"
        stoperr('naming', mesg, chem)

    lfg = 0
    codefg(chem, fgrp, lfg)
    
    if nfn > 0:
        iloc = srh5(chem, chemfn, nfn)
        if iloc > 0:
            name = namfn[iloc]
            iloc = srh5(name, namlst, nrec)
            if iloc > 0:
                print('--error-- in naming. Species: ', chem)
                print('is provided with a fixed name, which was')
                print('already used. Name = ', name)
                raise Exception("in naming")
            iloc = -iloc
            return

    il = 0
    for i in range(1, nord + 1):
        for j in range(1, lfg + 1):
            if fgrp[j - 1] == ord_str[i - 1]:
                name = name[:il] + fgrp[j - 1] + name[il + 1:]
                il += 1
                if il == 2:
                    break
        if il == 2:
            break
    
    if il == 0:
        mesg = " Cannot find 1st character"
        stoperr('naming', mesg, chem)
    
    if il == 1:
        name = name[:1] + '0' + name[2:]
    
    j = ic % 10
    name = name[:2] + str(j) + name[3:]
    
    name = name[:3] + '000'
    iloc = srh5(name, namlst, nrec)
    
    lockfn = True
    if iloc <= 0:
        iloc = -iloc
        if name[:3] != namlst[iloc][:3]:
            for l in range(1, nfn + 1):
                if name == namfn[l - 1]:
                    ii = 1
                    jj = 1
                    kk = 1
                    lockfn = False
                    break
            if lockfn:
                return
    
    if lockfn:
        for k in range(iloc + 1, nrec + 2):
            if namlst[k - 1][:3] != name[:3]:
                iloc = k
                break
        iloc = iloc - 1
        
        if namlst[iloc - 1][0] == ' ' and name[3:6] == '000':
            ii = 1
            jj = 1
            kk = 1
        else:
            ii = alfa.find(namlst[iloc - 1][3]) + 1
            jj = alfa.find(namlst[iloc - 1][4]) + 1
            kk = alfa.find(namlst[iloc - 1][5]) + 1
        
        if ii == 0 or jj == 0 or kk == 0:
            print('--error-- in naming. One of the character is not set. ')
            print('name in the list: ', namlst[iloc - 1])
            print('name search: ', name)
            print('chem: ', chem)
            mesg = 'One of the character is not set.'
            stoperr('naming', mesg, chem)
            raise Exception("in naming")
    
    while True:
        kk += 1
        if kk > nalfa:
            kk = 2
            jj += 1
        if jj > nalfa:
            jj = 1
            ii += 1
        if ii > nalfa:
            print('--error-- in naming. All the slot are occupied.')
            print('name:', name, ' chem: ', chem)
            raise Exception("in naming")
        
        name = name[:3] + alfa[ii - 1] + alfa[jj - 1] + alfa[kk - 1]
        
        name_exists = False
        for l in range(1, nfn + 1):
            if name == namfn[l - 1]:
                name_exists = True
                break
        
        if not name_exists:
            break

def codefg(chem, fgrp, lfg):
    global mxring
    
    tchem = chem
    lfg = 1
    fgrp = ' '
    locriegee = False
    
    ic, ih, in_, io, ir, is_, ifl, ibr, icl = 0, 0, 0, 0, 0, 0, 0, 0, 0
    getatoms(tchem, ic, ih, in_, io, ir, is_, ifl, ibr, icl)
    
    if ir != 0:
        if ir > 1:
            if '.(OO.)' in tchem:
                locriegee = True
            if '.(ZOO.)' in tchem:
                locriegee = True
            if '.(EOO.)' in tchem:
                locriegee = True
            if not locriegee:
                mesg = "Illegal di-radical "
                stoperr('codefg', mesg, chem)
        
        fgrp = '0'
        if locriegee:
            fgrp = '4'
        elif 'CO(OO.)' in tchem:
            fgrp = '3'
        elif '(OO.)' in tchem:
            fgrp = '2'
        elif '(O.)' in tchem:
            fgrp = '1'
        
        fgrp = fgrp + '.'
        lfg = 3
    
    nring = 0
    rjs = [[0, 0] for _ in range(mxring)]
    if '12' in tchem or 'C2' in tchem or '-O2' in tchem:
        fgrp = fgrp[:lfg] + 'TT'
        lfg += 2
        nring = 2
    elif 'C1' in tchem or '-O1' in tchem:
        fgrp = fgrp[:lfg] + 'T'
        lfg += 1
        nring = 1
    
    if 'c2' in tchem:
        fgrp = fgrp[:lfg] + 'RR'
        lfg += 2
        nring = 2
    elif 'c1' in tchem:
        fgrp = fgrp[:lfg] + 'R'
        lfg += 1
        nring = 1
    
    if nring > 0:
        rjsrm(nring, tchem, rjs)
    
    if ifl > 0:
        fgrp = fgrp[:lfg] + 'F'
        lfg += 1
    if ibr > 0:
        fgrp = fgrp[:lfg] + 'B'
        lfg += 1
    if icl > 0:
        fgrp = fgrp[:lfg] + 'L'
        lfg += 1
    
    i = tchem.find('CdO')
    if i != -1:
        fgrp = fgrp[:lfg] + 'X'
        lfg += 1
    
    nc = len(tchem.rstrip())
    
    np = 0
    for i in range(1, nc - 8 + 1):
        if tchem[i - 1:i + 8] == 'CO(OONO2)':
            fgrp = fgrp[:lfg] + 'P'
            lfg += 1
            np += 1
    
    na = 0
    for i in range(1, nc - 5 + 1):
        if tchem[i - 1:i + 5] == 'CO(OH)':
            fgrp = fgrp[:lfg] + 'A'
            lfg += 1
            na += 1
    
    ng = 0
    for i in range(1, nc - 6 + 1):
        if tchem[i - 1:i + 6] == 'CO(OOH)':
            fgrp = fgrp[:lfg] + 'G'
            lfg += 1
            ng += 1
    
    for i in range(3, nc - 4 + 1):
        if tchem[i - 1:i + 4] == '(OOH)' and tchem[i - 3:i - 1] != 'CO':
            fgrp = fgrp[:lfg] + 'H'
            lfg += 1
    
    for i in range(3, nc - 5 + 1):
        if tchem[i - 1:i + 5] == '(OOOH)' and tchem[i - 3:i - 1] != 'CO':
            fgrp = fgrp[:lfg] + 'Z'
            lfg += 1
    
    for i in range(3, nc - 5 + 1):
        if tchem[i - 1:i + 5] == '(ONO2)' or (tchem[i - 1:i + 6] == '(OONO2)' and tchem[i - 3:i - 1] != 'CO'):
            fgrp = fgrp[:lfg] + 'N'
            lfg += 1
    
    for i in range(1, nc + 1):
        if tchem[i - 1] == '=':
            fgrp = fgrp[:lfg] + 'U'
            lfg += 1
    
    for i in range(1, nc - 2 + 1):
        if tchem[i - 1:i + 2] == 'CHO':
            fgrp = fgrp[:lfg] + 'D'
            lfg += 1
    
    for i in range(1, nc - 4 + 1):
        if tchem[i - 1:i + 4] == '(NO2)':
            fgrp = fgrp[:lfg] + 'V'
            lfg += 1
    
    nf = 0
    for i in range(1, nc - 1 + 1):
        if tchem[i - 1:i + 1] == 'CO':
            nf += 1
    nf = nf - np - na - ng
    if nf > 0:
        for i in range(1, nf + 1):
            fgrp = fgrp[:lfg] + 'K'
            lfg += 1
    
    for i in range(1, nc - 3 + 1):
        if tchem[i:i + 4] == '(OH)' and tchem[i - 1:i + 4] != 'O(OH)':
            fgrp = fgrp[:lfg] + 'O'
            lfg += 1
    
    for i in range(1, nc - 2 + 1):
        if tchem[i - 1:i + 2] == '-O-':
            fgrp = fgrp[:lfg] + 'E'
            lfg += 1
    
    hcstring = 'CH1234().'
    lohc = True
    for i in range(1, nc + 1):
        if tchem[i - 1] not in hcstring:
            lohc = False
            break
    
    if lohc and lfg == 1:
        fgrp = 'C'
    
    if fgrp.strip() == '':
        mesg = "No letter code found."
        stoperr('codefg', mesg, chem)
    
    if any(c in fgrp for c in '056789FBLWSM'):
        mesg = "Letter code currently not allowed."
        stoperr('codefg', mesg, chem)