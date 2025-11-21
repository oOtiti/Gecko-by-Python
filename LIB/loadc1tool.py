def loadc1mch():
    from keyflag import pamfg
    from keyparameter import mxlco, mxlfo, mxlfl, dirgecko
    from dictstackdb import nrec, ninorg, inorglst, dict, namlst
    from sortstring import sort_string

    mxtprec = 200
    oname = [''] * mxtprec
    ochem = [''] * mxtprec
    ofgrp = [''] * mxtprec
    iname = [''] * mxtprec
    ichem = [''] * mxtprec
    ifgrp = [''] * mxtprec
    
    nrec = 1
    oname[0] = '######'
    ochem[0] = '###################'

    filename = dirgecko + "DATA/dic_inorg.dat"
    rddict(filename, ninorg, iname, ichem, ifgrp, nrec, oname, ochem, ofgrp)

    filename = dirgecko + "DATA/dic_c1.dat"
    rddict(filename, ninorg, iname, ichem, ifgrp, nrec, oname, ochem, ofgrp)

    if pamfg:
        filename = dirgecko + "DATA/dic_OFR.dat"
        rddict(filename, ninorg, iname, ichem, ifgrp, nrec, oname, ochem, ofgrp)

    if nrec > len(dict):
        print('error in loadc1mch, number of records exceeds mni')
        raise Exception("in rddict")
    if ninorg > len(inorglst):
        print('error in loadc1mch, number of records exceeds mni')
        raise Exception("in loadc1mch")

    for i in range(ninorg):
        inorglst[i] = f"{iname[i]:6s}   {ichem[i]:120s}  {ifgrp[i]:15s}"

    for i in range(nrec):
        dict[i] = f"{oname[i]:6s}   {ochem[i]:120s}  {ofgrp[i]:15s}"

    for i in range(nrec):
        namlst[i] = oname[i]
    sort_string(namlst[0:nrec])

    print('  ...reading inorganic reactions')
    filename = dirgecko + "DATA/mch_inorg.dat"
    rdfixmch(filename)

    if pamfg:
        filename = dirgecko + 'DATA/mch_OFR.dat'
        rdfixmch(filename)
    
    print('  ...reading CH4 chemistry')
    filename = dirgecko + "DATA/mch_singlec.dat"
    rdfixmch(filename)

    print('  ...writing CH3O2+counters reactions')
    rxmeo2ro2()


def rddict(filename, ninorg, iname, ichem, ifgrp, nrec, oname, ochem, ofgrp):
    from keyparameter import tfu1
    from atomtool import cnum
    from sortstring import sort_string
    from searching import srh5

    lenlin = 200
    line = ''
    tpline = ''
    i = 0
    ilin = 0
    ierr = 0
    ipos = 0
    loerr = 0
    nca = 0
    tpname = ''
    tpchem = ''
    tpfgrp = ''
    tponame = [''] * len(oname)
    tpochem = [''] * len(ochem)
    tpofgrp = [''] * len(ofgrp)
    tpiname = [''] * len(iname)
    tpichem = [''] * len(ichem)
    tpifgrp = [''] * len(ifgrp)

    for i in range(len(tponame)):
        tponame[i] = oname[i]
    for i in range(len(tpochem)):
        tpochem[i] = ochem[i]
    for i in range(len(tpofgrp)):
        tpofgrp[i] = ofgrp[i]
    for i in range(len(tpiname)):
        tpiname[i] = iname[i]
    for i in range(len(tpichem)):
        tpichem[i] = ichem[i]
    for i in range(len(tpifgrp)):
        tpifgrp[i] = ifgrp[i]

    try:
        f = open(filename, 'r')
    except:
        print('--error--, while reading file: ', filename)
        print('keyword "END" missing ?')
        raise Exception("in rddict")

    ilin = 0
    rdloop = True
    while rdloop:
        ilin = ilin + 1
        line = f.readline()
        if not line:
            print('--error--, while reading file: ', filename)
            print('at line number: ', ilin)
            print('keyword "END" missing ?')
            f.close()
            raise Exception("in rddict")
        if line[:3] == 'END':
            rdloop = False
            break
        if line[0] == '!':
            continue

        loerr = 0
        tpline = line.strip()
        parts = tpline.split()
        if len(parts) < 2:
            loerr = 1
        else:
            tpname = parts[0]
            tpchem = parts[1]
            if len(parts) > 2:
                tpfgrp = parts[2]
            else:
                tpfgrp = ' '
        
        if loerr != 0:
            print('--error--, while reading file: ', filename)
            print('at line number: ', ilin)
            print('can not read parameters in :', line)
            f.close()
            raise Exception("in rddict")

        nca = cnum(tpchem)
        if nca > 1:
            print('--error--, while reading file: ', filename)
            print('following species have more than 1 C:')
            print(line)
            f.close()
            raise Exception("in rddict")

        if nca == 1:
            if tpchem[0] != 'C':
                print('-error-, while reading file: ', filename)
                print('organic species does not start with C:', tpchem)
                print(line)
                f.close()
                raise Exception("in rddict")
            nrec = nrec + 1
            tponame[nrec-1] = tpname
            tpochem[nrec-1] = tpchem
            tpofgrp[nrec-1] = tpfgrp
        else:
            ninorg = ninorg + 1
            tpiname[ninorg-1] = tpname
            tpichem[ninorg-1] = tpchem
            tpifgrp[ninorg-1] = tpfgrp

    f.close()

    if nrec > len(oname) or ninorg > len(iname):
        print('--error--, while reading file: ', filename)
        print('too many species in the input file (nrec,ninorg > mxtprec)')
        raise Exception("in rddict")

    for i in range(len(ochem)):
        ochem[i] = tpochem[i]
    sort_string(ochem[0:nrec])

    ierr = 0
    for i in range(nrec - 1):
        if ochem[i] == ochem[i+1]:
            print('--error--, while reading file: ', filename)
            print('Following species identified 2 times: ', ochem[i])
            ierr = 1
    if ierr != 0:
        raise Exception("in rddict")

    for i in range(nrec):
        ilin = srh5(tpochem[i], ochem, nrec)
        if ilin <= 0:
            print('--error--, while sorting species in: ', filename)
            print('species "lost" after sorting the list: ', tpochem[i])
            raise Exception("in rddict")
        oname[ilin-1] = tponame[i]
        ofgrp[ilin-1] = tpofgrp[i]

    for i in range(len(ichem)):
        ichem[i] = tpichem[i]
    sort_string(ichem[0:ninorg])

    ierr = 0
    for i in range(ninorg - 1):
        if ichem[i] == ichem[i+1]:
            print('--error--, while reading file: ', filename)
            print('Following species identified 2 times: ', ichem[i])
            ierr = 1
    if ierr != 0:
        raise Exception("in rddict")

    for i in range(ninorg):
        ilin = srh5(tpichem[i], ichem, ninorg)
        if ilin <= 0:
            print('--error--, while sorting species in: ', filename)
            print('species "lost" after sorting the list: ', tpichem[i])
            raise Exception("in rddict")
        iname[ilin-1] = tpiname[i]
        ifgrp[ilin-1] = tpifgrp[i]


def rdfixmch(filename):
    from keyparameter import tfu1, mecu, refu
    from keyflag import wrtref
    from references import mxlcod
    from rxwrttool import count4rxn

    lenlin = 100
    line = ''
    line1 = ''
    line2 = ''
    reaction = ''
    info = ''
    comline = ''
    comtab = [''] * 3
    label = 0
    ilin = 0
    A = 0.0
    n = 0.0
    E = 0.0
    m = 0.0
    fact = 0.0
    one = 0.0
    zero = 0.0
    F0_300 = 0.0
    Finf_300 = 0.0
    E0 = 0.0
    Einf = 0.0
    Fc1 = 0.0
    Fc2 = 0.0
    Fc3 = 0.0
    Fc4 = 0.0
    n1 = 0
    n2 = 0
    loerr = 0
    ierr = 0
    idrx = 0

    one = 1.0
    zero = 0.0
    loerr = 0

    try:
        f = open(filename, 'r')
    except:
        print('--error--, while reading file: ', filename)
        print('keyword "END" missing ?')
        raise Exception("in rdfixmch")

    ilin = 0
    rdloop = True
    while rdloop:
        ilin = ilin + 1
        line = f.readline()
        if not line:
            print('--error--, while reading file: ', filename)
            print('at line number: ', ilin)
            print('keyword "END" missing ?')
            f.close()
            raise Exception("in rdfixmch")
        if line[:3] == 'END':
            rdloop = False
            break
        if line[0] == '!':
            continue

        n1 = line.find("!")
        if n1 > 0:
            line = line[:n1]

        n1 = line.find(':')
        n2 = line.find(';')
        if n1 == -1 or n2 == -1:
            loerr = 1
        if loerr != 0:
            print('--error--, missing ":" or ";" while reading: ', filename)
            print('at line number: ', ilin)
            print(line)
            f.close()
            raise Exception("in rdfixmch")

        if 'FALLOFF' in line:
            try:
                values = line[n1+1:n2].split()
                Fc1 = float(values[0])
                Fc2 = float(values[1])
                Fc3 = float(values[2])
                Fc4 = float(values[3])
            except:
                print('--error--, while reading in ', filename)
                print('at line number: ', ilin)
                print('while reading # in: ', line)
                f.close()
                raise Exception("in rdfixmch")
            
            reaction = line[:n1]
            comline = line[n2:]
            line_list = list(line)
            line_list[n2] = '/'
            line = ''.join(line_list)
            info = line[n1+1:n2+1]
            getcomcod(filename, comline, comtab)

            ilin = ilin + 1
            line1 = f.readline()
            try:
                values = line1.split()
                F0_300 = float(values[0])
                n = float(values[1])
                E0 = float(values[2])
            except:
                print('--error--, while reading in ', filename)
                print('at line number: ', ilin)
                print('while reading # in: ', line1)
                f.close()
                raise Exception("in rdfixmch")

            ilin = ilin + 1
            line2 = f.readline()
            try:
                values = line2.split()
                Finf_300 = float(values[0])
                m = float(values[1])
                Einf = float(values[2])
            except:
                print('--error--, while reading:', filename)
                print('at line number: ', ilin)
                print('while reading # in: ', line2)
                f.close()
                raise Exception("in rdfixmch")

            line1 = f'  FALLOFF /{F0_300:10.3e}{-n:4.1f}{E0:7.0f}{Fc1:6.1f}{Fc2:6.1f}{Fc3:6.1f}{Fc4:6.1f}/'

            count4rxn(3)
            print(f"{reaction:70s}                   {Finf_300:10.3e} {-m:4.1f} {Einf:7.0f}", file=mecu)
            print(line1, file=mecu)
            if wrtref:
                print(f"{reaction:70s}                   {Finf_300:10.3e} {-m:4.1f} {Einf:7.0f}", file=refu)
                print(line1, file=refu)

        elif 'EXTRA' in line:
            reaction = line[:n1]
            try:
                values = line[n1+1:n2].split()
                A = float(values[0])
                n = float(values[1])
                E = float(values[2])
            except:
                print('--error--, while reading:', filename)
                print('at line number: ', ilin)
                print('while reading # in: ', line)
                f.close()
                raise Exception("in rdfixmch")
            comline = line[n2:]
            getcomcod(filename, comline, comtab)

            ilin = ilin + 1
            line1 = f.readline()
            n1 = line1.find(':')
            n2 = line1.find(';')
            if n1 == -1 or n2 == -1:
                loerr = 1
            if loerr != 0:
                print('--error--, missing ":" or ";" while reading: ', filename)
                print(line1)
                f.close()
                raise Exception("in rdfixmch")
            try:
                label = int(line1[:n1])
            except:
                print('--error--, while reading:', filename)
                print('at line number: ', ilin)
                print('while reading label in: ', line1)
                f.close()
                raise Exception("in rdfixmch")
            line1_list = list(line1)
            line1_list[n2] = '/'
            line1 = ''.join(line1_list)
            info = line1[n1+1:n2+1]

            count4rxn(2)
            print(f"{reaction:70s}                   {A:10.3e} {n:4.1f} {E:7.0f}", file=mecu)
            print(f"  EXTRA /{label:4d}{info:50s}", file=mecu)
            if wrtref:
                print(f"{reaction:70s}                   {A:10.3e} {n:4.1f} {E:7.0f}", file=refu)
                print(f"  EXTRA /{label:4d}{info:50s}", file=refu)

        elif 'HV' in line:
            reaction = line[:n1]
            try:
                values = line[n1+1:n2].split()
                label = int(values[0])
                fact = float(values[1])
            except:
                print('--error--, while reading:', filename)
                print('while reading label & fact in: ', line)
                f.close()
                raise Exception("in rdfixmch")
            comline = line[n2:]
            getcomcod(filename, comline, comtab)

            count4rxn(1)
            print(f"{reaction:70s}                   {one:10.3e} {zero:4.1f} {zero:7.0f}", file=mecu)
            print(f"  HV / {label:5d} {fact:5.2f} /", file=mecu)
            if wrtref:
                print(f"{reaction:70s}                   {one:10.3e} {zero:4.1f} {zero:7.0f}", file=refu)
                print(f"  HV / {label:5d} {fact:5.2f} /", file=refu)

        else:
            reaction = line[:n1]
            try:
                values = line[n1+1:n2].split()
                A = float(values[0])
                n = float(values[1])
                E = float(values[2])
            except:
                print('--error--, while reading:', filename)
                print('while reading # in: ', line)
                f.close()
                raise Exception("in rdfixmch")
            comline = line[n2:]
            getcomcod(filename, comline, comtab)

            idrx = 0
            if 'TBODY' in reaction:
                idrx = 4
            count4rxn(idrx)
            print(f"{reaction:70s}                   {A:10.3e} {n:4.1f} {E:7.0f}", file=mecu)
            if wrtref:
                print(f"{reaction:70s}                   {A:10.3e} {n:4.1f} {E:7.0f}", file=refu)

    f.close()


def rxmeo2ro2():
    from keyparameter import mxlco, mxpd, mecu
    from keyflag import rx_ro2_multiclass
    from rxwrttool import rxwrit, rxinit

    p = [''] * mxpd
    r = [''] * 3
    xlabel = 0.0
    folow = [0.0] * 3
    fotroe = [0.0] * 4
    s = [0.0] * mxpd
    arrh = [0.0] * 3
    rrad = 0.0
    rmol1 = 0.0
    rmol2 = 0.0
    itype = 0
    idreac = 0
    nlabel = 0

    kwdclass = ['PERO1 ', 'PERO2 ', 'PERO3 ', 'PERO4 ', 'PERO5 ',
                'MEPERO', 'PERO7 ', 'PERO8 ', 'PERO9 ']

    ro2dat = [
        [1.00E-13,  0.,  627.],
        [1.00E-13,  0.,  341.],
        [1.00E-13,  0., -257.],
        [1.00E-13,  0., -289.],
        [1.00E-13,  0., -635.],
        [2.06E-13,  0., -365.],
        [1.00E-13,  0., -702.],
        [1.00E-13,  0., -936.],
        [2.00E-12,  0., -508.]
    ]

    stoi_m = [
        [0.8 , 0.2 , 0.0],
        [0.6 , 0.2 , 0.2],
        [0.8 , 0.2 , 0.0],
        [0.6 , 0.2 , 0.2],
        [0.8 , 0.2 , 0.0],
        [0.370, 0.315, 0.315],
        [0.6 , 0.2 , 0.2],
        [0.6 , 0.2 , 0.2],
        [0.8 , 0.2 , 0.0]
    ]

    for i in range(len(r)):
        r[i] = ' '
    for i in range(len(folow)):
        folow[i] = 0.0
    for i in range(len(fotroe)):
        fotroe[i] = 0.0
    for i in range(len(s)):
        s[i] = 0.0
    for i in range(len(p)):
        p[i] = ' '
    idreac = 0
    nlabel = 0

    if rx_ro2_multiclass:
        for itype in range(9):
            rrad = stoi_m[itype][0]
            rmol1 = stoi_m[itype][1]
            rmol2 = stoi_m[itype][2]

            rxinit(r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)

            r[0] = 'CH3O2'
            r[1] = kwdclass[itype]

            arrh[0] = ro2dat[itype][0]
            arrh[1] = ro2dat[itype][1]
            arrh[2] = ro2dat[itype][2]

            s[0] = rrad
            p[0] = 'CH3O '
            s[1] = rmol1
            p[1] = 'CH2O '
            s[2] = rmol2
            p[2] = 'CH3OH'

            rxwrit(mecu, r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)
    else:
        itype = 5
        rrad = stoi_m[itype][0]
        rmol1 = stoi_m[itype][1]
        rmol2 = stoi_m[itype][2]

        rxinit(r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)

        r[0] = 'CH3O2'
        r[1] = 'PERO1'

        arrh[0] = ro2dat[5][0]
        arrh[1] = ro2dat[5][1]
        arrh[2] = ro2dat[5][2]

        s[0] = rrad
        p[0] = 'CH3O '
        s[1] = rmol1
        p[1] = 'CH2O '
        s[2] = rmol2
        p[2] = 'CH3OH'

        rxwrit(mecu, r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)


def getcomcod(filename, comline, comtab):
    from searching import srh5
    from references import nreac_info, code

    maxcom = len(comtab)
    ncpbeg = 1
    ierr = 0
    i = 0
    n1 = 0

    for i in range(len(comtab)):
        comtab[i] = ' '

    if comline[0] != ";":
        print('--error--, while reading file: ', filename)
        print('expected 1st char (;) not found in:', comline)
        raise Exception("in getcomcod")

    parts = comline[1:].split(';')
    for i in range(min(len(parts), maxcom)):
        if parts[i].strip():
            comtab[i] = parts[i].strip()

    ierr = 0
    for i in range(maxcom):
        if comtab[i][0] != ' ':
            n1 = srh5(comtab[i], code, nreac_info)
            if n1 <= 0:
                print('--error--, while reading file ', filename)
                print('unknow comment: ', comtab[i])
                ierr = 1
    if ierr != 0:
        raise Exception("in getcomcod")