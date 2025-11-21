def changephase():
    from keyparameter import mxlfo, mxlco, mecu, prtu, walu, tfu1, dirout
    from dictstackdb import ncha, chatab, nwpspe
    from keyflag import g2pfg, g2wfg, chafg
    from dictstackdb import nrec, dict
    from atomtool import cnum, molweight

    i = 0
    j = 0
    chem = ''
    idnam = ''
    molmass = 0.0

    if chafg:
        f = open(dirout + 'chaname.dat', 'w')

    for i in range(1, nrec):
        chem = dict[i][9:129].strip()
        if '.' in chem:
            continue
        if cnum(chem) < 2:
            continue

        if chem[0] != 'C' and chem[0] != 'c' and chem[:2] != '-O':
            if chem[:3] == '#mm':
                chem = chem[3:]
            elif chem[0] == '#':
                if 'CH2OO' in chem:
                    continue
                if '.(OO.)' in chem:
                    continue
                chem = chem[1:]
            else:
                continue
        
        molweight(chem, molmass)
        idnam = dict[i][:6].strip()

        if chafg:
            found_cha = False
            for j in range(ncha):
                if idnam == chatab[j]:
                    if chafg:
                        print(f"G{chatab[j]:6s}", file=f)
                    if g2pfg:
                        nwpspe = nwpspe + 1
                        print(f"A{idnam:6s}          /{molmass:6.1f}/", file=prtu)
                    if g2wfg:
                        nwpspe = nwpspe + 1
                        print(f"W{idnam:6s}          /{molmass:6.1f}/", file=walu)
                    found_cha = True
                    break
            if found_cha:
                continue

        if g2pfg:
            gas2aero(mecu, prtu, idnam, molmass)
        if g2wfg:
            gas2wall(mecu, walu, chem, idnam, molmass)

    if chafg:
        f.close()


def gas2aero(lrea, ldic, idnam, mweight):
    from rxwrttool import rxinit, rxwrit_dyn
    from dictstackdb import nwpspe

    mxprod = 4
    mxreac = 3
    r = [''] * mxreac
    p = [''] * mxprod
    s = [0.0] * mxprod
    arrh = [0.0] * 3
    idreac = 0
    auxinfo = [0.0] * 9
    charfrom = ''
    charto = ''

    for i in range(len(r)):
        r[i] = ' '
    for i in range(len(s)):
        s[i] = 0.0
    for i in range(len(p)):
        p[i] = ' '
    for i in range(len(arrh)):
        arrh[i] = 0.0
    idreac = 0
    for i in range(len(auxinfo)):
        auxinfo[i] = 0.0

    nwpspe = nwpspe + 1
    print(f"A{idnam:6s}          /{mweight:6.1f}/", file=ldic)

    r[0] = idnam
    r[1] = 'AIN '
    p[0] = idnam
    s[0] = 1.0
    arrh[0] = 1.0
    idreac = 1
    charfrom = 'G'
    charto = 'A'
    rxwrit_dyn(lrea, r, s, p, arrh, idreac, auxinfo, charfrom, charto)

    r[0] = idnam
    r[1] = 'AOU '
    p[0] = idnam
    s[0] = 1.0
    arrh[0] = 1.0
    idreac = 2
    charfrom = 'A'
    charto = 'G'
    rxwrit_dyn(lrea, r, s, p, arrh, idreac, auxinfo, charfrom, charto)


def gas2wall(lrea, ldic, chem, idnam, mweight):
    from dictstackdb import nwpspe
    from rxwrttool import rxinit, rxwrit_dyn

    mxprod = 4
    mxreac = 3
    r = [''] * mxreac
    p = [''] * mxprod
    s = [0.0] * mxprod
    arrh = [0.0] * 3
    idreac = 0
    auxinfo = [0.0] * 9
    charfrom = ''
    charto = ''

    for i in range(len(r)):
        r[i] = ' '
    for i in range(len(s)):
        s[i] = 0.0
    for i in range(len(p)):
        p[i] = ' '
    for i in range(len(arrh)):
        arrh[i] = 0.0
    idreac = 0
    for i in range(len(auxinfo)):
        auxinfo[i] = 0.0

    nwpspe = nwpspe + 1
    print(f"W{idnam:6s}          /{mweight:6.1f}/", file=ldic)

    r[0] = idnam
    r[1] = 'WIN '
    p[0] = idnam
    s[0] = 1.0
    arrh[0] = 1.111E-3
    idreac = 3
    charfrom = 'G'
    charto = 'W'
    if 'O' not in chem:
        auxinfo[0] = 2.E-5
    else:
        auxinfo[0] = 1.2E-4
    rxwrit_dyn(lrea, r, s, p, arrh, idreac, auxinfo, charfrom, charto)

    r[0] = idnam
    r[1] = 'WOU '
    p[0] = idnam
    s[0] = 1.0
    arrh[0] = 1.111E-3
    idreac = 4
    charfrom = 'W'  
    charto = 'G'
    if 'O' not in chem:
        auxinfo[0] = 2.E-5
    else:
        auxinfo[0] = 1.2E-4
    rxwrit_dyn(lrea, r, s, p, arrh, idreac, auxinfo, charfrom, charto)