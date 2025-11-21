def isomer(chem, brtio, chg, tabinfo):
    from dictstackdb import nrec, dict, stabl, dbrch, mxcri, mxiso, diccri
    from keyparameter import mxlfo, mxlfl
    from namingtool import codefg
    from toolbox import stoperr
    
    chg = 0
    lfg = 1
    fgrp = ' '
    dicflg = [0] * mxiso
    score = [0] * mxiso
    dictinfo = [0] * mxiso
    tabinfo = [0] * mxcri
    progname='isomer'
    interacgrp(chem, tabinfo)
    
    if stabl < 2:
        return
    
    codefg(chem, fgrp, lfg)
    nlen = fgrp.index(' ')
    
    niso = 0
    for i in range(1, nrec + 1):
        if tabinfo[0] != diccri[i-1][0]:
            continue
            
        match = True
        for j in range(nlen):
            if fgrp[j] != dict[i-1][mxlfo + 11 + j:mxlfo + 11 + j + 1]:
                match = False
                break
        if not match:
            continue
            
        niso += 1
        if niso > mxiso:
            mesg = "isomer # > mxiso"
            stoperr(progname, mesg, chem)
        dictinfo[niso-1] = i
    
    if niso == 0:
        return
    
    tiso = niso
    
    if brtio >= 3E-2:
        return
    
    if brtio < 5E-2:
        if niso > 1:
            for kc in range(6, 10):
                seliso(kc, tiso, dictinfo, tabinfo, dicflg, niso)
                
        if niso > 1:
            seliso(10, tiso, dictinfo, tabinfo, dicflg, niso)
            
        if niso > 1:
            seliso(11, tiso, dictinfo, tabinfo, dicflg, niso)
            
        if niso > 1:
            seliso(2, tiso, dictinfo, tabinfo, dicflg, niso)
            
        if niso > 1:
            for kc in range(17, 14, -1):
                seliso(kc, tiso, dictinfo, tabinfo, dicflg, niso)
                
        if niso > 1:
            seliso(14, tiso, dictinfo, tabinfo, dicflg, niso)
            
        if niso > 1:
            for kc in range(13, 11, -1):
                seliso(kc, tiso, dictinfo, tabinfo, dicflg, niso)
                
        if niso > 1:
            for kc in range(18, 22):
                seliso(kc, tiso, dictinfo, tabinfo, dicflg, niso)
                
        if niso > 1:
            for kc in range(3, 6):
                seliso(kc, tiso, dictinfo, tabinfo, dicflg, niso)
    
    if niso == 1:
        for i in range(tiso):
            if dicflg[i] == 0:
                j = dictinfo[i]
                chem = dict[j-1][9:mxlfo + 10]
                dbrch[j-1] = dbrch[j-1] + brtio
                chg = 1
                return
                
    elif niso > 1:
        imax = 0
        ymax = 0.0
        for i in range(tiso):
            if dicflg[i] == 0:
                j = dictinfo[i]
                if dbrch[j-1] > ymax:
                    ymax = dbrch[j-1]
                    imax = j
        chem = dict[imax-1][9:mxlfo + 10]
        dbrch[imax-1] = dbrch[imax-1] + brtio
        chg = 1
        return


def interacgrp(chem, tabinfo):
    from keyparameter import mxnode, mxlgr, mxring
    from rjtool import rjgrm
    from stdgrbond import grbond
    
    group = [''] * mxnode
    nodetype = [''] * mxnode
    tabgrp = [1] * mxnode
    bond = [[0] * mxnode for _ in range(mxnode)]
    dbflg = 0
    nring = 0
    rjg = [[0] * 2 for _ in range(mxring)]
    
    grbond(chem, group, bond, dbflg, nring)
    
    for i in range(len(bond)):
        bond[i][i] = 0
        
    if nring > 0:
        rjgrm(nring, group, rjg)
        
    ng = sum(1 for g in group if g != '')
    tabinfo[0] = ng
    
    for i in range(ng):
        if group[i][:2] == 'CO':
            nodetype[i] = 'y'
        elif group[i][:3] == 'CHO':
            nodetype[i] = 'y'
        elif group[i][:1] == 'c':
            nodetype[i] = 'r'
        elif group[i][:3] == '-O-':
            nodetype[i] = 'o'
        else:
            nodetype[i] = 'n'
            
    for i in range(ng):
        if group[i][:4] == 'CH3 ':
            tabgrp[i] = 0
            tabinfo[1] += 1
        elif group[i][:4] == 'CH2 ':
            tabgrp[i] = 0
            tabinfo[2] += 1
        elif group[i][:3] == 'CH ':
            tabgrp[i] = 0
            tabinfo[3] += 1
        elif group[i][:2] == 'C ':
            tabgrp[i] = 0
            tabinfo[4] += 1
            
    for i in range(ng):
        isum = 0
        for j in range(ng):
            if bond[i][j] != 0:
                isum += 1
        if isum == 1:
            tabinfo[5] += 1
        if isum == 2:
            tabinfo[6] += 1
        if isum == 3:
            tabinfo[7] += 1
        if isum == 4:
            tabinfo[8] += 1
            
    for i in range(ng - 1):
        if nodetype[i] == 'y':
            for j in range(i + 1, ng):
                if bond[i][j] != 0:
                    if nodetype[j] == 'y':
                        tabinfo[9] += 1
                        if group[j] == 'CHO ':
                            tabinfo[17] += 1
                        if group[j] == 'CO(OONO2) ':
                            tabinfo[18] += 1
                        if group[j] == 'CO(OH) ':
                            tabinfo[19] += 1
                        if group[j] == 'CO(OOH) ':
                            tabinfo[20] += 1
    tabinfo[13] = tabinfo[17] + tabinfo[18] + tabinfo[19] + tabinfo[20]
    
    if 'CH2CH3' in chem:
        tabinfo[14] += 1
    if 'CH3CH2' in chem:
        tabinfo[14] += 1
    if 'CH2CH2CH3' in chem:
        tabinfo[15] += 1
    if 'CH3CH2CH2' in chem:
        tabinfo[15] += 1
    if 'CH2CH2CH2CH3' in chem:
        tabinfo[16] += 1
    if 'CH3CH2CH2CH2' in chem:
        tabinfo[16] += 1
        
    for i in range(ng):
        memo1 = 0
        memo2 = 0
        if tabgrp[i] != 0:
            for j in range(ng):
                if bond[i][j] != 0:
                    if tabgrp[j] != 0:
                        tabinfo[10] += 1
                        
                    for k in range(ng):
                        if bond[j][k] != 0 and k != i:
                            if tabgrp[k] != 0:
                                if group[i][:3] == '-O-' or group[k][:3] == '-O-':
                                    tabinfo[10] += 1
                                else:
                                    tabinfo[11] += 1
                                    
                            for l in range(ng):
                                if bond[k][l] != 0 and l != j:
                                    if tabgrp[l] != 0:
                                        if group[i][:3] == '-O-' or group[l][:3] == '-O-':
                                            tabinfo[11] += 1
                                        else:
                                            if i != memo1 and l != memo2:
                                                tabinfo[12] += 1
                                                memo1 = i
                                                memo2 = l
            tabgrp[i] = 0


def seliso(kc, tiso, dictinfo, tabinfo, dicflg, niso):
    from dictstackdb import diccri
    
    if niso == 1:
        return
        
    tempniso = niso
    tempdicflg = dicflg[:]
    
    for i in range(tiso):
        if dicflg[i] != 1:
            ij = dictinfo[i]
            if tabinfo[kc-1] != diccri[ij-1][kc-1]:
                dicflg[i] = 1
                niso -= 1
                
    if niso == 0:
        niso = tempniso
        for i in range(tiso):
            dicflg[i] = tempdicflg[i]
        return