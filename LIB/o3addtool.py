from keyflag import wrtsarinfo
from keyparameter import saru, mxcp
from ringtool import findring, ring_data
from mapping import gettrack

def o3rate_mono(bond, zebond, group, ngr, cdtable, cdsub, arrhc):
    arrhc[:] = [0.0] * len(arrhc)
    mult = 1.0
    ifct = 0
    Fring = 1.0
    k298 = 0.0
    cd1 = cdtable[0]
    cd2 = cdtable[1]
    
    if wrtsarinfo:
        print('  ')
        print(' ======= O3RATE MONO ====== ')
    
    ring = [0] * len(group)
    rngflg = 0
    findring(cd1, cd2, ngr, bond, rngflg, ring)
    
    if rngflg == 0:
        ifct = 0
        arrhc[0] = 2E-10
        for i in range(1, ngr+1):
            for j in range(1, 3):
                if ((bond[cdtable[j-1]][i-1] == 1 and (group[i-1][:2] == 'CO' or group[i-1][:3] == 'CHO')) or 
                    bond[cdtable[j-1]][i-1] == 3):
                    if j == 1:
                        Ci = cd1
                        Cf = cd2
                    else:
                        Ci = cd2
                        Cf = cd1
                    ifct = i
                    tmparrhc = [0.0] * 3
                    ref_rate_mono(bond, zebond, group, ngr, cdsub, tmparrhc, Cf, Ci, ifct)
                    arrhc[0] = min(arrhc[0], tmparrhc[0])
            if ifct != 0:
                return
    
    Ci = cd1
    Cf = cd2
    ref_rate_mono(bond, zebond, group, ngr, cdsub, arrhc, Cf, Ci, ifct)
    
    subs_fact(group, bond, ngr, ring, cdtable, mult)
    if rngflg == 1:
        ring_fact_mono(group, bond, cd1, ngr, Fring)
        for i in range(1, ngr+1):
            for j in range(1, 3):
                if (bond[i-1][cdtable[j-1]-1] == 1 and 
                    (group[i-1] == 'CO' or group[i-1] == 'CHO')):
                    Fring = Fring / 50.0
    
    arrhc[0] = arrhc[0] * mult * Fring
    
    if wrtsarinfo:
        print('=> Non-conjugated')
        print(f'mult={mult}')
        if rngflg == 1:
            print(f'Fring={Fring}')
        print(f'arrhc={arrhc[0:3]}')
        print(f'k298={k298}')
        print(' ')
        print(' ')

def o3rate_conj(bond, group, ngr, cdtable, cdsub, arrhc):
    arrhc[:] = [0.0] * len(arrhc)
    mult = 1.0
    Fring = 1.0
    
    if wrtsarinfo:
        print('  ')
        print(' ====== O3RATE CONJ ====== ')
    
    ref_rate_conj(cdtable, cdsub, arrhc)
    
    ring = [0] * len(group)
    rngflg = 0
    findring(cdtable[0], cdtable[1], ngr, bond, rngflg, ring)
    
    nbF = 0
    for i in range(1, ngr+1):
        for j in range(1, 5):
            if ((bond[i-1][cdtable[j-1]-1] == 1 and 
                 (group[i-1][:2] == 'CO' or group[i-1] == 'CHO')) or 
                bond[i-1][cdtable[j-1]-1] == 3):
                nbF += 1
    
    subs_fact(group, bond, ngr, ring, cdtable, mult)
    arrhc[0] = arrhc[0] * mult**0.5
    
    if nbF != 0:
        for i in range(1, ngr+1):
            for j in range(1, 5):
                if (bond[i-1][cdtable[j-1]-1] == 1 and 
                    (group[i-1] == 'CO' or group[i-1] == 'CHO')):
                    arrhc[0] = arrhc[0] / 50.0
    
    if (ring[cdtable[0]-1] == 1 and ring[cdtable[1]-1] == 1 and 
        ring[cdtable[2]-1] == 1 and ring[cdtable[3]-1] == 1):
        ring_fact_conj(group, bond, cdtable, ngr, Fring)
    elif ring[cdtable[0]-1] == 1 and ring[cdtable[1]-1] == 1:
        ring_fact_mono(group, bond, cdtable[0], ngr, Fring)
    elif ring[cdtable[2]-1] == 1 and ring[cdtable[3]-1] == 1:
        ring_fact_mono(group, bond, cdtable[2], ngr, Fring)
    
    arrhc[0] = arrhc[0] * Fring
    
    if wrtsarinfo:
        print('=> conjugated')
        for i in range(1, ngr+1):
            print(f'cdsub({i})={cdsub[i-1]}')
        print(f'mult={mult}')
        if rngflg == 1:
            print(f'Fring={Fring}')
        print(f'arrhc={arrhc[0:3]}')
        print(' ')
        print(' ')

def ref_rate_mono(bond, zebond, group, ngr, cdsub, arrhc, Cf, Ci, ifct):
    arrhc[:] = [0.0] * len(arrhc)
    
    if ifct != 0:
        if group[ifct-1][:3] == 'CO ':
            ester = False
            for i in range(1, ngr+1):
                if bond[ifct-1][i-1] == 3 and ifct != Ci:
                    ester = True
                    break
            
            if ester:
                if group[Ci-1][:3] == 'CdH' and group[Cf-1][:4] == 'CdH2':
                    arrhc[0] = 0.15E-17
                else:
                    arrhc[0] = 0.65E-17
            else:
                if group[Cf-1][:4] == 'CdH2':
                    if group[Ci-1][:3] == 'CdH':
                        arrhc[0] = 0.52E-17
                    else:
                        arrhc[0] = 1.20E-17
                elif group[Cf-1][:3] == 'CdH':
                    arrhc[0] = 3.90E-17
                else:
                    arrhc[0] = 0.83E-17
        
        elif group[ifct-1][:3] == 'CHO':
            if group[Ci-1][:3] == 'CdH':
                if group[Cf-1][:3] == 'CdH':
                    arrhc[0] = 0.14E-17
                else:
                    arrhc[0] = 0.18E-17
            else:
                if group[Cf-1][:4] == 'CdH2':
                    arrhc[0] = 0.12E-17
                else:
                    arrhc[0] = 0.57E-17
        
        elif group[ifct-1][:3] == '-O-':
            ester = False
            for i in range(1, ngr+1):
                if bond[ifct-1][i-1] == 3 and group[i-1][:3] == 'CO ':
                    ester = True
                    break
            
            if ester:
                if group[Ci-1][:3] == 'CdH':
                    arrhc[0] = 0.32E-17
                else:
                    arrhc[0] = 0.054E-17
            else:
                if group[Ci-1][:4] == 'CdH ':
                    if group[Cf-1][:4] == 'CdH2':
                        arrhc[0] = 17E-17
                    else:
                        arrhc[0] = 42E-17
                else:
                    arrhc[0] = 1.3E-17
        
        elif group[ifct-1][:6] == 'CO(OH)':
            arrhc[0] = 0.23E-17
        
        elif (group[ifct-1][:7] == 'CO(OOH)' or 
              group[ifct-1][:9] == 'CO(OONO2)'):
            if group[Ci-1][:3] == 'CdH' and group[Cf-1][:4] == 'CdH2':
                arrhc[0] = 0.14E-17
            else:
                arrhc[0] = 0.65E-17
        
        if group[ifct-1] == 'CO ' or group[ifct-1] == 'CHO':
            for i in range(1, ngr+1):
                if (bond[Cf-1][i-1] == 1 and 
                    (group[i-1] == 'CO ' or group[i-1] == 'CHO')):
                    arrhc[0] = 0.50E-17
        
        if group[ifct-1] == '-O-':
            for i in range(1, ngr+1):
                if (bond[Cf-1][i-1] == 3 or 
                    (bond[Ci-1][i-1] == 3 and i != ifct)):
                    arrhc[0] = 48E-17
        
        alpha = 0.19
        
        fact = 1.0
        for j in range(1, ngr+1):
            ncarb = 0
            if bond[Cf-1][j-1] == 1 or bond[Cf-1][j-1] == 3:
                ntr = 0
                track = [[0] * len(bond) for _ in range(mxcp)]
                trlen = [0] * mxcp
                gettrack(bond, j, ngr, ntr, track, trlen)
                ch_eff = [False] * len(group)
                for i_idx in range(1, ntr+1):
                    if track[i_idx-1][1] == Cf:
                        continue
                    for k in range(1, trlen[i_idx-1]+1):
                        if group[track[i_idx-1][k-1]-1][:3] != '-O-':
                            ch_eff[track[i_idx-1][k-1]-1] = True
                    if group[track[i_idx-1][0]-1][:3] == 'CO ':
                        ch_eff[track[i_idx-1][0]-1] = False
                ncarb = sum(ch_eff)
                
                if ncarb >= 1:
                    fact = fact + alpha * (ncarb - 1)
                
                if wrtsarinfo and ncarb != 0:
                    print(f'group(ifct)= {group[ifct-1]}')
                    print(f'chain length effect of substituent of Cf:{group[Cf-1]}')
                    print(f'arrhc(1)={arrhc[0]}')
                    print(f'fact={fact}')
                    print(f'alpha={alpha}')
                    print(f'ncarb={ncarb}')
        arrhc[0] = arrhc[0] * fact
        
        fact = 1.0
        for j in range(1, ngr+1):
            if j == ifct:
                continue
            ncarb = 0
            if bond[Ci-1][j-1] == 1:
                ntr = 0
                track = [[0] * len(bond) for _ in range(mxcp)]
                trlen = [0] * mxcp
                gettrack(bond, j, ngr, ntr, track, trlen)
                ch_eff = [False] * len(group)
                for i_idx in range(1, ntr+1):
                    if track[i_idx-1][1] == Ci:
                        continue
                    for k in range(1, trlen[i_idx-1]+1):
                        ch_eff[track[i_idx-1][k-1]-1] = True
                ncarb = sum(ch_eff)
                
                if ncarb >= 1:
                    fact = fact + alpha * (ncarb - 1)
                
                if wrtsarinfo and ncarb != 0:
                    print(f'chain length effect of substituent of Ci:{group[Ci-1]}')
                    print(f'fact={fact}')
                    print(f'alpha={alpha}')
                    print(f'ncarb={ncarb}')
        arrhc[0] = arrhc[0] * fact
        
        fact = 1.0
        for j in range(1, ngr+1):
            if j == Ci:
                continue
            ncarb = 0
            if bond[ifct-1][j-1] == 1 or bond[ifct-1][j-1] == 3:
                ntr = 0
                track = [[0] * len(bond) for _ in range(mxcp)]
                trlen = [0] * mxcp
                gettrack(bond, j, ngr, ntr, track, trlen)
                ch_eff = [False] * len(group)
                for i_idx in range(1, ntr+1):
                    if track[i_idx-1][1] == ifct:
                        continue
                    for k in range(1, trlen[i_idx-1]+1):
                        if group[track[i_idx-1][k-1]-1][:3] != '-O-':
                            ch_eff[track[i_idx-1][k-1]-1] = True
                ncarb = sum(ch_eff)
                
                if bond[ifct-1][j-1] == 3 and group[j-1][:3] == 'CO ':
                    ncarb = ncarb - 1
                if ncarb >= 1:
                    fact = fact + alpha * (ncarb - 1)
                
                if wrtsarinfo and ncarb != 0:
                    print(f'chain length effect of substituent of {group[ifct-1]}')
                    print(f'arrhc(1)={arrhc[0]}')
                    print(f'fact={fact}')
                    print(f'alpha={alpha}')
                    print(f'ncarb={ncarb}')
        arrhc[0] = arrhc[0] * fact
    
    else:
        nb = cdsub[Ci-1] + cdsub[Cf-1]
        if nb == 0:
            arrhc[0] = 9.14E-15
            arrhc[2] = 2580.0
        elif nb == 1:
            arrhc[0] = 2.91E-15
            arrhc[2] = 1690.0
        elif nb == 2:
            if cdsub[Ci-1] == 1:
                if zebond[Ci-1][Cf-1] == 1:
                    arrhc[0] = 3.39E-15
                    arrhc[2] = 995.0
                elif zebond[Ci-1][Cf-1] == 2:
                    arrhc[0] = 7.29E-15
                    arrhc[2] = 1120.0
                else:
                    arrhc[0] = 1.50E-16
                    arrhc[2] = 0.0
            else:
                arrhc[0] = 4.00E-15
                arrhc[2] = 1685.0
        elif nb == 3:
            arrhc[0] = 7.61E-15
            arrhc[2] = 830.0
        elif nb == 4:
            arrhc[0] = 3.00E-15
            arrhc[2] = 300.0
        else:
            print('-error1-- in ref_rate, no rate constant found')
            print(f'ref_rate, nb={nb}')
            raise Exception("in ref_rate")

def subs_fact(group, bond, ngr, ring, cdtable, mult):
    mult = 1.0
    for i in range(1, len(cdtable)+1):
        if cdtable[i-1] == 0:
            continue
        for j in range(1, ngr+1):
            if bond[cdtable[i-1]-1][j-1] == 1:
                if '(OH)' in group[j-1] and group[j-1] != 'CO(OH)':
                    mult = mult * 1.4
                if '(ONO2)' in group[j-1]:
                    mult = mult * 0.044
                nalkyl = 0
                nalkyl_cyc = 0
                for k in range(1, ngr+1):
                    if bond[k-1][j-1] == 1 and k != cdtable[i-1]:
                        if (group[k-1][:4] == 'CH3 ' or group[k-1][:3] == 'CH2' or 
                            group[k-1][:3] == 'CH ' or group[k-1][:3] == 'CH(' or 
                            group[k-1][:2] == 'C '):
                            if ring[k-1] == 0:
                                nalkyl += 1
                            else:
                                nalkyl_cyc += 1
                        elif group[k-1][:3] == 'CO ':
                            mult = mult * 0.32
                        elif group[k-1][:3] == 'CHO':
                            mult = mult * 0.32
                    elif bond[k-1][j-1] == 3:
                        for l in range(1, ngr+1):
                            if bond[l-1][k-1] == 3 and group[l-1][:3] == 'CO ':
                                mult = mult * 0.25
                                break
                        mult = mult * 0.6
                if nalkyl >= 1:
                    if nalkyl_cyc == 0:
                        mult = mult * (0.54 ** (nalkyl - 1))
                    else:
                        mult = mult * (0.54 ** nalkyl)

def ring_fact_mono(group, bond, cd1, ngr, Fring):
    mxirg = 6
    nring_ind = 0
    trackrg = [[0] * len(group) for _ in range(mxirg)]
    ring_ind = [[False] * len(group) for _ in range(mxirg)]
    
    ring_data(cd1, ngr, bond, group, nring_ind, ring_ind, trackrg)
    
    Fring = 1.0
    for i in range(1, nring_ind+1):
        mult = 1.0
        rgord = sum(ring_ind[i-1])
        if rgord == 5:
            mult = 3.9
        elif rgord == 6:
            mult = 0.52
        elif rgord == 7:
            mult = 2.0
        elif rgord == 8:
            mult = 2.8
        elif rgord == 9:
            mult = 2.1
        elif rgord == 10:
            mult = 0.24
        elif rgord == 11:
            mult = 12.0
        Fring = Fring * mult

def ref_rate_conj(cdtable, cdsub, arrhc):
    arrhc[0] = 1E-14
    arrhc[1] = 0.0
    arrhc[2] = 0.0
    
    ic1 = cdtable[0]
    ic2 = cdtable[1]
    ic3 = cdtable[2]
    ic4 = cdtable[3]
    nb = cdsub[ic1-1] + cdsub[ic2-1] + cdsub[ic3-1] + cdsub[ic4-1]
    
    if nb == 3:
        if cdsub[ic1-1] == 1 or cdsub[ic4-1] == 1:
            arrhc[2] = 1677.0
        elif cdsub[ic1-1] == 0 and cdsub[ic4-1] == 0:
            arrhc[2] = 1980.0
    elif nb == 4:
        if ((cdsub[ic1-1] == 0 and cdsub[ic4-1] == 0) or 
            (cdsub[ic1-1] == 0 and cdsub[ic3-1] == 0) or 
            (cdsub[ic2-1] == 0 and cdsub[ic3-1] == 0) or 
            (cdsub[ic2-1] == 0 and cdsub[ic4-1] == 0)):
            arrhc[2] = 1774.0
        if ((cdsub[ic2-1] == 2 and cdsub[ic3-1] == 1 and cdsub[ic4-1] == 1) or 
            (cdsub[ic1-1] == 1 and cdsub[ic2-1] == 1 and cdsub[ic3-1] == 2) or 
            (cdsub[ic2-1] == 1 and cdsub[ic3-1] == 2 and cdsub[ic4-1] == 1) or 
            (cdsub[ic1-1] == 1 and cdsub[ic2-1] == 2 and cdsub[ic1-1] == 1) or 
            (cdsub[ic1-1] == 0 and cdsub[ic2-1] == 1 and cdsub[ic3-1] == 1 and cdsub[ic4-1] == 2)):
            arrhc[2] = 1439.0
        if (cdsub[ic1-1] == 1 and cdsub[ic2-1] == 1 and 
            cdsub[ic3-1] == 1 and cdsub[ic4-1] == 1):
            arrhc[2] = 1008.0
    elif nb == 5:
        if cdsub[ic1-1] == 0 or cdsub[ic4-1] == 0:
            arrhc[2] = 1214.0
        if cdsub[ic1-1] != 0 and cdsub[ic4-1] != 0:
            arrhc[2] = 686.0
    elif nb == 6:
        if cdsub[ic1-1] == 0 or cdsub[ic4-1] == 0:
            arrhc[2] = 887.0
        if cdsub[ic1-1] != 0 and cdsub[ic4-1] != 0:
            arrhc[2] = 359.0
    elif nb == 7:
        arrhc[2] = 142.0
    elif nb == 8:
        arrhc[2] = 0.0

def ring_fact_conj(group, bond, cdtable, ngr, Fring):
    mxirg = 6
    nring_ind = 0
    trackrg = [[0] * len(group) for _ in range(mxirg)]
    ring_ind = [[False] * len(group) for _ in range(mxirg)]
    
    ring_data(cdtable[0], ngr, bond, group, nring_ind, ring_ind, trackrg)
    
    Fring = 1.0
    for i in range(1, nring_ind+1):
        mult = 1.0
        rgord = sum(ring_ind[i-1])
        if rgord == 6:
            mult = 4.50
        elif rgord == 7:
            mult = 0.44
        elif rgord == 8:
            mult = 0.06
        Fring = Fring * mult