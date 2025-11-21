def nodmap(bond, top, nca, npos, nnod, tnod):
    from keyparameter import mxcp

    track = [[0] * len(bond) for _ in range(mxcp)]
    trlen = [0] * mxcp
    ntr = 0
    i = 0
    j = 0

    nnod = 0
    for i in range(len(tnod)):
        tnod[i] = 0

    gettrack(bond, top, nca, ntr, track, trlen)

    if ntr > 1:
        for i in range(1, ntr):
            if track[i-1][npos-1] == 0:
                continue
            for j in range(i+1, ntr+1):
                if track[j-1][npos-1] != 0:
                    if track[i-1][npos-1] == track[j-1][npos-1]:
                        track[j-1][npos-1] = 0

    for i in range(1, ntr + 1):
        if track[i-1][npos-1] != 0:
            nnod = nnod + 1
            tnod[nnod-1] = track[i-1][npos-1]


def gettrack(bond, top, nca, ntr, track, trlen):
    ptr = 0
    niv = 0
    nod = 0
    i = 0
    j = 0
    k = 0
    memo = [0] * len(bond)
    slope = 0
    mxtr = len(track)

    ptr = 0
    for i in range(len(trlen)):
        trlen[i] = 0
    for i in range(len(track)):
        for j in range(len(track[0])):
            track[i][j] = 0
    for i in range(len(memo)):
        memo[i] = 0

    track[0][0] = top
    memo[0] = top
    niv = 1
    ptr = 0
    ntr = 1
    nod = top
    slope = 1

    nextnode = True
    while nextnode:
        ptr = ptr + 1

        if ptr > nca:
            if niv == 1:
                nextnode = False
                break
            else:
                ptr = memo[niv-1]
                memo[niv-1] = 0
                niv = niv - 1
                nod = memo[niv-1]

                if slope != 0:
                    ntr = ntr + 1
                    if ntr > mxtr:
                        print('--error-- in gettrack. ntr is greater than mxtr')
                        raise Exception("in gettrack")
                    for j in range(niv):
                        track[ntr-1][j] = track[ntr-2][j]
                track[ntr-1][niv] = 0
                slope = 0
                continue

        if bond[ptr-1][nod-1] == 0:
            continue

        if bond[ptr-1][nod-1] != 0:
            for k in range(niv):
                if ptr == memo[k]:
                    continue
            niv = niv + 1
            track[ntr-1][niv-1] = ptr
            nod = ptr
            memo[niv-1] = nod
            ptr = 0
            slope = 1
            continue

    for j in range(len(track[0])):
        track[ntr-1][j] = 0
    ntr = ntr - 1

    for i in range(1, ntr + 1):
        for j in range(nca, 0, -1):
            if track[i-1][j-1] != 0:
                trlen[i-1] = j
                break


def abcde_map(bond, top, nca, nabcde, tabcde):
    mxdeep = len(tabcde)
    mxtr = len(tabcde[0])
    track = [[0] * len(tabcde[0][0]) for _ in range(mxtr)]
    trlen = [0] * mxtr
    ntr = 0
    ltabc = [[[0] * len(tabcde[0][0]) for _ in range(mxtr)] for _ in range(mxdeep)]
    i = 0
    j = 0
    k = 0
    ii = 0

    for i in range(len(nabcde)):
        nabcde[i] = 0
    for i in range(len(tabcde)):
        for j in range(len(tabcde[0])):
            for k in range(len(tabcde[0][0])):
                tabcde[i][j][k] = 0
    for i in range(len(ltabc)):
        for j in range(len(ltabc[0])):
            for k in range(len(ltabc[0][0])):
                ltabc[i][j][k] = 0

    gettrack(bond, top, nca, ntr, track, trlen)

    if ntr > mxtr:
        print('--error-- in abcde_map, ntr > mxtr')
        raise Exception("in abcde_map")

    for k in range(2, mxdeep + 1):
        for i in range(1, ntr + 1):
            if trlen[i-1] >= k:
                for j in range(1, k + 1):
                    ltabc[k-1][i-1][j-1] = track[i-1][j-1]

    for k in range(2, mxdeep + 1):
        for i in range(1, ntr):
            if ltabc[k-1][i-1][k-1] == 0:
                continue
            for ii in range(i+1, ntr + 1):
                for j in range(1, k + 1):
                    if ltabc[k-1][i-1][j-1] != ltabc[k-1][ii-1][j-1]:
                        break
                else:
                    for j in range(1, k + 1):
                        ltabc[k-1][ii-1][j-1] = 0

    for k in range(2, mxdeep + 1):
        for i in range(1, ntr + 1):
            if ltabc[k-1][i-1][k-1] != 0:
                nabcde[k-1] = nabcde[k-1] + 1
                for j in range(1, k + 1):
                    ii = nabcde[k-1]
                    tabcde[k-1][ii-1][j-1] = ltabc[k-1][i-1][j-1]

    nabcde[0] = 1
    tabcde[0][0][0] = top


def ciptree(chem, bond, group, ngr, nabcde, tabcde, ciptr):
    from searching import search_ipos
    from toolbox import countstring, stoperr

    nciptr = [0] * len(ciptr)
    mxbd = 5
    nC = [0] * mxbd
    nO = [0] * mxbd
    nH = [0] * mxbd
    nN = [0] * mxbd
    i = 0
    idepth = 0
    ig = 0
    nnod = 0
    ipos = 0
    pvnod = 0
    nfun = 0
    ii = 0
    nel = 0
    irk = 0
    ilast = 0

    progname = 'ciptree'
    mesg = ''

    for i in range(len(ciptr)):
        for j in range(len(ciptr[0])):
            ciptr[i][j] = 0
    for i in range(len(nciptr)):
        nciptr[i] = 0

    for idepth in range(1, ngr + 1):
        nnod = nabcde[idepth-1]
        if nnod == 0:
            break

        for i in range(len(nC)):
            nC[i] = 0
            nH[i] = 0
            nO[i] = 0
            nN[i] = 0

        for i in range(1, nnod + 1):
            ig = tabcde[idepth-1][i-1][idepth-1]

            if group[ig-1][:3] == 'CH3':
                nC[0] = nC[0] + 1
                nH[1] = nH[1] + 3
            elif group[ig-1][:3] == 'CH2':
                nC[0] = nC[0] + 1
                nH[1] = nH[1] + 2
            elif group[ig-1][:3] == 'CHO':
                nC[0] = nC[0] + 1
                nH[1] = nH[1] + 1
                nO[1] = nO[1] + 2
            elif group[ig-1][:2] == 'CO':
                nC[0] = nC[0] + 1
                nO[1] = nO[1] + 2
            elif group[ig-1][:2] == 'CH':
                nC[0] = nC[0] + 1
                nH[1] = nH[1] + 1
            elif group[ig-1][:4] == 'CdH2':
                nC[0] = nC[0] + 1
                nH[1] = nH[1] + 2
            elif group[ig-1][:3] == 'CdH':
                nC[0] = nC[0] + 1
                nH[1] = nH[1] + 1
            elif group[ig-1][:3] == 'CdO':
                nC[0] = nC[0] + 1
                nO[1] = nO[1] + 2
            elif group[ig-1][:1] == 'C':
                nC[0] = nC[0] + 1
            elif group[ig-1][:3] == '-O-':
                nO[0] = nO[0] + 1
            elif group[ig-1][:2] == 'cH':
                nC[0] = nC[0] + 1
                nH[1] = nH[1] + 1
            elif group[ig-1][:1] == 'c':
                nC[0] = nC[0] + 1

            if group[ig-1][:2] == 'Cd':
                if idepth > 1:
                    pvnod = tabcde[idepth-1][i-1][idepth-2]
                    if bond[pvnod-1][ig-1] == 2:
                        nC[0] = nC[0] + 1

            nfun = countstring(group[ig-1], '(OH)')
            if nfun != 0:
                nO[1] = nO[1] + nfun
                nH[2] = nH[2] + nfun

            nfun = countstring(group[ig-1], '(OOH)')
            if nfun != 0:
                nO[1] = nO[1] + nfun
                nO[2] = nO[2] + nfun
                nH[3] = nH[3] + nfun

            nfun = countstring(group[ig-1], '(OOOH)')
            if nfun != 0:
                nO[1] = nO[1] + nfun
                nO[2] = nO[2] + nfun
                nO[3] = nO[3] + nfun
                nH[4] = nH[4] + nfun

            nfun = countstring(group[ig-1], '(ONO2)')
            if nfun != 0:
                nO[1] = nO[1] + nfun
                nN[2] = nN[2] + nfun
                nO[3] = nO[3] + 2 * nfun

            nfun = countstring(group[ig-1], '(OONO2)')
            if nfun != 0:
                nO[1] = nO[1] + nfun
                nO[2] = nO[2] + nfun
                nN[3] = nN[3] + nfun
                nO[4] = nO[4] + 2 * nfun

            nfun = countstring(group[ig-1], '(NO2)')
            if nfun != 0:
                nN[1] = nN[1] + nfun
                nO[2] = nO[2] + 2 * nfun

            nfun = countstring(group[ig-1], '(OO.)')
            if nfun != 0:
                nO[1] = nO[1] + nfun
                nO[2] = nO[2] + nfun

            nfun = countstring(group[ig-1], '(O.)')
            if nfun != 0:
                nO[1] = nO[1] + nfun

        for ii in range(1, mxbd + 1):
            if nC[ii-1] == 0:
                continue
            nel = nC[ii-1]
            irk = idepth + ii - 1
            ipos = search_ipos(12, ciptr[irk-1])
            ilast = nciptr[irk-1]
            if ipos <= ilast:
                for j in range(ilast, ipos-1, -1):
                    ciptr[irk-1][j+nel-1] = ciptr[irk-1][j-1]
            for j in range(ipos, ipos + nel):
                ciptr[irk-1][j-1] = 12
            nciptr[irk-1] = nciptr[irk-1] + nel

        for ii in range(1, mxbd + 1):
            if nO[ii-1] == 0:
                continue
            nel = nO[ii-1]
            irk = idepth + ii - 1
            ipos = search_ipos(16, ciptr[irk-1])
            ilast = nciptr[irk-1]
            if ipos <= ilast:
                for j in range(ilast, ipos-1, -1):
                    ciptr[irk-1][j+nel-1] = ciptr[irk-1][j-1]
            for j in range(ipos, ipos + nel):
                ciptr[irk-1][j-1] = 16
            nciptr[irk-1] = nciptr[irk-1] + nel

        for ii in range(1, mxbd + 1):
            if nN[ii-1] == 0:
                continue
            nel = nN[ii-1]
            irk = idepth + ii - 1
            ipos = search_ipos(14, ciptr[irk-1])
            ilast = nciptr[irk-1]
            if ipos <= ilast:
                for j in range(ilast, ipos-1, -1):
                    ciptr[irk-1][j+nel-1] = ciptr[irk-1][j-1]
            for j in range(ipos, ipos + nel):
                ciptr[irk-1][j-1] = 14
            nciptr[irk-1] = nciptr[irk-1] + nel

        for ii in range(2, mxbd + 1):
            if nH[ii-1] == 0:
                continue
            nel = nH[ii-1]
            irk = idepth + ii - 1
            ipos = search_ipos(1, ciptr[irk-1])
            ilast = nciptr[irk-1]
            if ipos <= ilast:
                for j in range(ilast, ipos-1, -1):
                    ciptr[irk-1][j+nel-1] = ciptr[irk-1][j-1]
            for j in range(ipos, ipos + nel):
                ciptr[irk-1][j-1] = 1
            nciptr[irk-1] = nciptr[irk-1] + nel

    for i in range(1, len(nciptr) + 1):
        if nciptr[i-1] == 0:
            break
        if nciptr[i-1] > len(ciptr[0]):
            mesg = "too many elements added in a line of CIP tree"
            stoperr(progname, mesg, chem)


def alkenetrack(chem, bond, group, ngr, ncdtrack, cdtracklen, cdtrack):
    from keyparameter import mxcp
    from toolbox import stoperr

    cdbond = [[0] * len(group) for _ in range(len(group))]
    cdpst = [0] * len(group)
    i = 0
    j = 0
    ncd = 0
    nprim = 0
    mxcd = len(cdtrack[0])

    track = [[0] * len(group) for _ in range(mxcp)]
    trlen = [0] * mxcp
    ntr = 0

    progname = 'alkenetrack '
    mesg = ''

    ncdtrack = 0
    for i in range(len(cdtrack)):
        for j in range(len(cdtrack[0])):
            cdtrack[i][j] = 0
    for i in range(len(cdtracklen)):
        cdtracklen[i] = 0

    ncd = 0
    for i in range(1, ngr + 1):
        if group[i-1][:2] == 'Cd':
            ncd = ncd + 1
            for j in range(1, ngr + 1):
                if j == i:
                    continue
                if bond[i-1][j-1] != 0:
                    if group[j-1][:2] == 'Cd':
                        cdbond[i-1][j-1] = 1
                        cdbond[j-1][i-1] = 1
    if ncd == 0:
        return

    for i in range(1, ngr + 1):
        cdpst[i-1] = sum(cdbond[i-1])

    nprim = 0
    for i in range(1, ngr + 1):
        if cdpst[i-1] == 1:
            nprim = nprim + 1
        if cdpst[i-1] > 2:
            mesg = "tertiary Cd structure identified (and not allowed)"
            stoperr(progname, mesg, chem)
    if nprim == 0:
        mesg = "unexpected cyclic Cd structure identified"
        stoperr(progname, mesg, chem)

    for i in range(1, ngr + 1):
        if cdpst[i-1] == 1:
            gettrack(cdbond, i, ngr, ntr, track, trlen)
            if ntr > 1:
                mesg = "unexpected branching identified on Cd structure"
                stoperr(progname, mesg, chem)
            if trlen[0] > 4:
                mesg = "More than 4 Cd identified (and not allowed)"
                stoperr(progname, mesg, chem)
            if trlen[0] == 3:
                mesg = ">C=C=C< structure identified and not allowed"
                stoperr(progname, mesg, chem)
            ncdtrack = ncdtrack + 1
            if ncdtrack > len(cdtrack):
                mesg = "maximum number of Cd tracks reached."
                stoperr(progname, mesg, chem)
            for j in range(1, mxcd + 1):
                cdtrack[ncdtrack-1][j-1] = track[0][j-1]
            cdtracklen[ncdtrack-1] = trlen[0]
            cdpst[track[0][trlen[0]-1]] = 0


def estertrack(chem, bond, group, ngr, netrack, etracklen, etrack):
    from keyparameter import mxcp
    from toolbox import stoperr

    ebond = [[0] * len(bond) for _ in range(len(bond))]
    i = 0
    j = 0
    nest = 0
    ntr = 0
    nprim = 0
    mxenod = len(etrack[0])

    track = [[0] * len(group) for _ in range(mxcp)]
    trlen = [0] * mxcp
    progname = 'estertrack'
    mesg = ''

    netrack = 0
    for i in range(len(etrack)):
        for j in range(len(etrack[0])):
            etrack[i][j] = 0
    for i in range(len(etracklen)):
        etracklen[i] = 0

    nest = 0
    for i in range(1, ngr + 1):
        if group[i-1][:3] == '-O-':
            for j in range(1, ngr + 1):
                if j == i:
                    continue
                if bond[i-1][j-1] != 0:
                    if group[j-1] == 'CO':
                        ebond[i-1][j-1] = 1
                        ebond[j-1][i-1] = 1
                        nest = nest + 1
                    elif group[j-1] == 'CHO':
                        ebond[i-1][j-1] = 1
                        ebond[j-1][i-1] = 1
                        nest = nest + 1
    if nest == 0:
        return

    epst = [0] * len(bond)
    for i in range(1, ngr + 1):
        epst[i-1] = sum(ebond[i-1])

    nprim = 0
    for i in range(1, ngr + 1):
        if epst[i-1] == 1:
            nprim = nprim + 1
        if epst[i-1] > 2:
            mesg = "tertiary ester structure identified (and not allowed)"
            stoperr(progname, mesg, chem)
    if nprim == 0:
        mesg = "unexpected cyclic ester structure identified"
        stoperr(progname, mesg, chem)

    for i in range(1, ngr + 1):
        if epst[i-1] == 1:
            gettrack(ebond, i, ngr, ntr, track, trlen)
            if ntr > 1:
                mesg = "unexpected branching identified on ester structure"
                stoperr(progname, mesg, chem)
            if trlen[0] > mxenod:
                mesg = "More than mx possible ester identified (and not allowed)"
                stoperr(progname, mesg, chem)
            netrack = netrack + 1
            if netrack > len(etrack):
                mesg = "maximum number of Cd tracks reached."
                stoperr(progname, mesg, chem)
            for j in range(1, mxenod + 1):
                etrack[netrack-1][j-1] = track[0][j-1]
            etracklen[netrack-1] = trlen[0]
            epst[track[0][trlen[0]-1]] = 0


def chemmap(chem, node, group, bond, ngrp, nodetype, alifun, cdfun, arofun, mapfun, funflg, tabester, nfcd, nfcr, ierr):
    i = 0
    j = 0
    k = 0
    lgr = 0
    nnod = 0
    tnod = [0] * len(bond)
    nf = 0
    ialpha = 0
    ialpha2 = 0
    iy = 0
    rflg = 0
    dflg = 0
    yflg = 0
    ichecko = 0
    ichecky = 0
    ytab = [0] * 2
    nester = 0

    ierr = 0
    nfcd = 0
    nfcr = 0
    ngrp = 0
    for i in range(len(alifun)):
        alifun[i] = 0.0
    for i in range(len(cdfun)):
        cdfun[i] = 0.0
    for i in range(len(arofun)):
        arofun[i] = 0.0
    for i in range(len(mapfun)):
        for j in range(len(mapfun[0])):
            for k in range(len(mapfun[0][0])):
                mapfun[i][j][k] = 0.0
    for i in range(len(nodetype)):
        nodetype[i] = ' '
    for i in range(len(funflg)):
        funflg[i] = 0
    nester = 0
    for i in range(len(tabester)):
        for j in range(len(tabester[0])):
            tabester[i][j] = 0
    lgr = len(group[0])

    for i in range(1, node + 1):
        if group[i-1][:2] == 'CO':
            nodetype[i-1] = 'y'
        elif group[i-1][:3] == 'CHO':
            nodetype[i-1] = 'y'
        elif group[i-1][:1] == 'c':
            nodetype[i-1] = 'r'
        elif group[i-1][:3] == '-O-':
            nodetype[i-1] = 'o'
        elif group[i-1][:2] == 'Cd':
            nodetype[i-1] = 'd'
        else:
            nodetype[i-1] = 'n'

    if '(OH)' in chem:
        for i in range(1, node + 1):
            if '(OH)' in group[i-1]:
                nf = 0
                for j in range(0, lgr - 3):
                    if group[i-1][j:j+4] == '(OH)':
                        nf = nf + 1
                if nodetype[i-1] == 'n':
                    alifun[0] = alifun[0] + nf
                    mapfun[i-1][0][0] = nf
                    ngrp = ngrp + nf
                elif nodetype[i-1] == 'r':
                    arofun[0] = arofun[0] + nf
                    mapfun[i-1][2][0] = nf
                    ngrp = ngrp + nf
                    nfcr = nfcr + nf
                elif nodetype[i-1] == 'd':
                    cdfun[0] = cdfun[0] + nf
                    mapfun[i-1][1][0] = nf
                    ngrp = ngrp + nf
                    nfcd = nfcd + nf
                elif nodetype[i-1] == 'y':
                    nodmap(bond, i, node, 2, nnod, tnod)
                    if nnod > 1:
                        for j in range(1, node + 1):
                            print('group', j, '-', group[j-1])
                        print('-- error --, a unique C is expected in')
                        print('alpha position of a carboxylic group')
                        raise Exception("in chemmap")
                    ialpha = tnod[0]
                    if nodetype[ialpha-1] == 'r':
                        arofun[10] = arofun[10] + 1
                        mapfun[i-1][2][10] = 1
                        ngrp = ngrp + 1
                        nfcr = nfcr + 1
                    elif nodetype[ialpha-1] == 'd':
                        cdfun[10] = cdfun[10] + 1
                        mapfun[i-1][1][10] = 1
                        ngrp = ngrp + 1
                        nfcd = nfcd + 1
                    else:
                        alifun[10] = alifun[10] + 1
                        mapfun[i-1][0][10] = 1
                        ngrp = ngrp + 1

    if '(NO2)' in chem:
        for i in range(1, node + 1):
            if '(NO2)' in group[i-1]:
                nf = 0
                for j in range(0, lgr - 4):
                    if group[i-1][j:j+5] == '(NO2)':
                        nf = nf + 1
                if nodetype[i-1] == 'n':
                    alifun[1] = alifun[1] + nf
                    mapfun[i-1][0][1] = nf
                    ngrp = ngrp + nf
                elif nodetype[i-1] == 'y':
                    alifun[1] = alifun[1] + nf
                    mapfun[i-1][0][1] = nf
                    ngrp = ngrp + nf
                elif nodetype[i-1] == 'd':
                    cdfun[1] = cdfun[1] + nf
                    mapfun[i-1][1][1] = nf
                    ngrp = ngrp + nf
                    nfcd = nfcd + nf
                elif nodetype[i-1] == 'r':
                    arofun[1] = arofun[1] + nf
                    mapfun[i-1][2][1] = nf
                    ngrp = ngrp + nf
                    nfcr = nfcr + nf
                else:
                    print('-- error --, a (NO2) group is borne by an ')
                    print(' unexpected group in chem :')
                    print(chem)
                    raise Exception("in chemmap")

    if '(ONO2)' in chem:
        for i in range(1, node + 1):
            if '(ONO2)' in group[i-1]:
                if 'CO(ONO2)' in group[i-1]:
                    continue
                nf = 0
                for j in range(0, lgr - 5):
                    if group[i-1][j:j+6] == '(ONO2)':
                        nf = nf + 1
                if nodetype[i-1] == 'n':
                    alifun[2] = alifun[2] + nf
                    mapfun[i-1][0][2] = nf
                    ngrp = ngrp + nf
                elif nodetype[i-1] == 'd':
                    cdfun[2] = cdfun[2] + nf
                    mapfun[i-1][1][2] = nf
                    ngrp = ngrp + nf
                    nfcd = nfcd + nf
                elif nodetype[i-1] == 'r':
                    arofun[2] = arofun[2] + nf
                    mapfun[i-1][2][2] = nf
                    ngrp = ngrp + nf
                    nfcr = nfcr + nf
                else:
                    print('-- error --, a (ONO2) group is borne by an ')
                    print('unexpected group in chem :')
                    print(chem)
                    raise Exception("in chemmap")

    if '(OOH)' in chem:
        for i in range(1, node + 1):
            if '(OOH)' in group[i-1]:
                nf = 0
                for j in range(0, lgr - 4):
                    if group[i-1][j:j+5] == '(OOH)':
                        nf = nf + 1
                if nodetype[i-1] == 'n':
                    alifun[3] = alifun[3] + nf
                    mapfun[i-1][0][3] = nf
                    ngrp = ngrp + nf
                elif nodetype[i-1] == 'r':
                    arofun[3] = arofun[3] + nf
                    mapfun[i-1][2][3] = nf
                    ngrp = ngrp + nf
                    nfcr = nfcr + nf
                elif nodetype[i-1] == 'd':
                    cdfun[3] = cdfun[3] + nf
                    mapfun[i-1][1][3] = nf
                    ngrp = ngrp + nf
                    nfcd = nfcd + nf
                elif nodetype[i-1] == 'y':
                    nodmap(bond, i, node, 2, nnod, tnod)
                    if nnod > 1:
                        for j in range(1, node + 1):
                            print('group', j, '-', group[j-1])
                        print('-- error --, a unique C is expected in')
                        print('alpha position of a peracid group')
                        print(chem)
                        raise Exception("in chemmap")
                    ialpha = tnod[0]
                    if nodetype[ialpha-1] == 'r':
                        arofun[11] = arofun[11] + 1
                        mapfun[i-1][2][11] = 1
                        ngrp = ngrp + 1
                    elif nodetype[ialpha-1] == 'd':
                        cdfun[11] = cdfun[11] + 1
                        mapfun[i-1][1][11] = 1
                        ngrp = ngrp + 1
                    else:
                        alifun[11] = alifun[11] + 1
                        mapfun[i-1][0][11] = 1
                        ngrp = ngrp + 1

    if '(F)' in chem:
        for i in range(1, node + 1):
            if '(F)' in group[i-1]:
                nf = 0
                for j in range(0, lgr - 2):
                    if group[i-1][j:j+3] == '(F)':
                        nf = nf + 1
                if nodetype[i-1] == 'n':
                    alifun[4] = alifun[4] + nf
                    mapfun[i-1][0][4] = nf
                    ngrp = ngrp + nf
                elif nodetype[i-1] == 'r':
                    arofun[4] = arofun[4] + nf
                    mapfun[i-1][2][4] = nf
                    ngrp = ngrp + nf
                    nfcr = nfcr + nf
                elif nodetype[i-1] == 'd':
                    cdfun[4] = cdfun[4] + nf
                    mapfun[i-1][1][4] = nf
                    ngrp = ngrp + nf
                    nfcd = nfcd + nf
                elif nodetype[i-1] == 'y':
                    nodmap(bond, i, node, 2, nnod, tnod)
                    if nnod > 1:
                        for j in range(1, node + 1):
                            print('group', j, '-', group[j-1])
                        print('-- error --, a unique C is expected in')
                        print('alpha position of a fluoro acyl group')
                        print(chem)
                        raise Exception("in chemmap")
                    ialpha = tnod[0]
                    if nodetype[ialpha-1] == 'r':
                        arofun[16] = arofun[16] + 1
                        mapfun[i-1][2][16] = 1
                        ngrp = ngrp + 1
                    elif nodetype[ialpha-1] == 'd':
                        cdfun[16] = cdfun[16] + 1
                        mapfun[i-1][1][16] = 1
                        ngrp = ngrp + 1
                    else:
                        alifun[16] = alifun[16] + 1
                        mapfun[i-1][0][16] = 1
                        ngrp = ngrp + 1

    if '(Cl)' in chem:
        for i in range(1, node + 1):
            if '(Cl)' in group[i-1]:
                nf = 0
                for j in range(0, lgr - 3):
                    if group[i-1][j:j+4] == '(Cl)':
                        nf = nf + 1
                if nodetype[i-1] == 'n':
                    alifun[5] = alifun[5] + nf
                    mapfun[i-1][0][5] = nf
                    ngrp = ngrp + nf
                elif nodetype[i-1] == 'r':
                    arofun[5] = arofun[5] + nf
                    mapfun[i-1][2][5] = nf
                    ngrp = ngrp + nf
                    nfcr = nfcr + nf
                elif nodetype[i-1] == 'd':
                    cdfun[5] = cdfun[5] + nf
                    mapfun[i-1][1][5] = nf
                    ngrp = ngrp + nf
                    nfcd = nfcd + nf
                elif nodetype[i-1] == 'y':
                    nodmap(bond, i, node, 2, nnod, tnod)
                    if nnod > 1:
                        for j in range(1, node + 1):
                            print('group', j, '-', group[j-1])
                        print('-- error --, a unique C is expected in')
                        print('alpha position of a chloro acyl group')
                        print(chem)
                        raise Exception("in chemmap")
                    ialpha = tnod[0]
                    if nodetype[ialpha-1] == 'r':
                        arofun[17] = arofun[17] + 1
                        mapfun[i-1][2][17] = 1
                        ngrp = ngrp + 1
                    elif nodetype[ialpha-1] == 'd':
                        cdfun[17] = cdfun[17] + 1
                        mapfun[i-1][1][17] = 1
                        ngrp = ngrp + 1
                    else:
                        alifun[17] = alifun[17] + 1
                        mapfun[i-1][0][17] = 1
                        ngrp = ngrp + 1

    if '(Br)' in chem:
        for i in range(1, node + 1):
            if '(Br)' in group[i-1]:
                nf = 0
                for j in range(0, lgr - 3):
                    if group[i-1][j:j+4] == '(Br)':
                        nf = nf + 1
                if nodetype[i-1] == 'n':
                    alifun[6] = alifun[6] + nf
                    mapfun[i-1][0][6] = nf
                    ngrp = ngrp + nf
                elif nodetype[i-1] == 'r':
                    arofun[6] = arofun[6] + nf
                    mapfun[i-1][2][6] = nf
                    ngrp = ngrp + nf
                    nfcr = nfcr + nf
                elif nodetype[i-1] == 'd':
                    cdfun[6] = cdfun[6] + nf
                    mapfun[i-1][1][6] = nf
                    ngrp = ngrp + nf
                    nfcd = nfcd + nf
                elif nodetype[i-1] == 'y':
                    nodmap(bond, i, node, 2, nnod, tnod)
                    if nnod > 1:
                        for j in range(1, node + 1):
                            print('group', j, '-', group[j-1])
                        print('-- error --, a unique C is expected in')
                        print('alpha position of a bromo acyl group')
                        print(chem[:70])
                        raise Exception("in chemmap")
                    ialpha = tnod[0]
                    if nodetype[ialpha-1] == 'r':
                        arofun[18] = arofun[18] + 1
                        mapfun[i-1][2][18] = 1
                        ngrp = ngrp + 1
                    elif nodetype[ialpha-1] == 'd':
                        cdfun[18] = cdfun[18] + 1
                        mapfun[i-1][1][18] = 1
                        ngrp = ngrp + 1
                    else:
                        alifun[18] = alifun[18] + 1
                        mapfun[i-1][0][18] = 1
                        ngrp = ngrp + 1

    if '(I)' in chem:
        for i in range(1, node + 1):
            if '(I)' in group[i-1]:
                nf = 0
                for j in range(0, lgr - 2):
                    if group[i-1][j:j+3] == '(I)':
                        nf = nf + 1
                if nodetype[i-1] == 'n':
                    alifun[7] = alifun[7] + nf
                    mapfun[i-1][0][7] = nf
                    ngrp = ngrp + nf
                elif nodetype[i-1] == 'r':
                    arofun[7] = arofun[7] + nf
                    mapfun[i-1][2][7] = nf
                    ngrp = ngrp + nf
                    nfcr = nfcr + nf
                elif nodetype[i-1] == 'd':
                    cdfun[7] = cdfun[7] + nf
                    mapfun[i-1][1][7] = nf
                    ngrp = ngrp + nf
                    nfcd = nfcd + nf
                elif nodetype[i-1] == 'y':
                    nodmap(bond, i, node, 2, nnod, tnod)
                    if nnod > 1:
                        for j in range(1, node + 1):
                            print('group', j, '-', group[j-1])
                        print('-- error --, a unique C is expected in')
                        print('alpha position of a Iodo acyl group')
                        print(chem[:70])
                        raise Exception("in chemmap")
                    ialpha = tnod[0]
                    if nodetype[ialpha-1] == 'r':
                        arofun[19] = arofun[19] + 1
                        mapfun[i-1][2][19] = 1
                        ngrp = ngrp + 1
                        nfcr = nfcr + 1
                    elif nodetype[ialpha-1] == 'd':
                        cdfun[19] = cdfun[19] + 1
                        mapfun[i-1][1][19] = 1
                        ngrp = ngrp + 1
                        nfcd = nfcd + 1
                    else:
                        alifun[19] = alifun[19] + 1
                        mapfun[i-1][0][19] = 1
                        ngrp = ngrp + 1

    if '-O' in chem:
        for i in range(1, node + 1):
            if group[i-1][:3] == '-O-':
                nodmap(bond, i, node, 2, nnod, tnod)
                rflg = 0
                dflg = 0
                yflg = 0
                for j in range(1, nnod + 1):
                    ialpha = tnod[j-1]
                    if nodetype[ialpha-1] == 'y':
                        ichecky = 0
                        for k in range(4):
                            if tabester[k][1] == ialpha:
                                ichecky = 1
                        if ichecky == 0:
                            yflg = yflg + 1
                            ytab[yflg-1] = ialpha
                if yflg == 0:
                    continue

                if yflg != 0:
                    nester = nester + 1
                    if nester > 4:
                        print("in chemmap, nester > 4")
                        raise Exception("in chemmap")
                    iy = ytab[0]
                    tabester[nester-1][0] = i
                    tabester[nester-1][1] = iy
                    for j in range(1, nnod + 1):
                        if tnod[j-1] != ytab[0]:
                            ialpha = tnod[j-1]
                    if group[iy-1][:3] == 'CHO':
                        if nodetype[ialpha-1] == 'r':
                            arofun[15] = arofun[15] + 1
                            mapfun[i-1][2][15] = 1
                            mapfun[iy-1][2][15] = 1
                            ngrp = ngrp + 1
                            nfcr = nfcr + 1
                        elif nodetype[ialpha-1] == 'd':
                            cdfun[15] = cdfun[15] + 1
                            mapfun[i-1][1][15] = 1
                            mapfun[iy-1][1][15] = 1
                            ngrp = ngrp + 1
                            nfcd = nfcd + 1
                        else:
                            alifun[15] = alifun[15] + 1
                            mapfun[i-1][0][15] = 1
                            mapfun[iy-1][0][15] = 1
                            ngrp = ngrp + 1
                    elif group[iy-1][:3] == 'CO ':
                        nodmap(bond, iy, node, 2, nnod, tnod)
                        for j in range(1, nnod + 1):
                            if tnod[j-1] != i:
                                ialpha2 = tnod[j-1]
                        rflg = 0
                        dflg = 0

                        if nodetype[ialpha-1] == 'r':
                            rflg = rflg + 1
                        if nodetype[ialpha-1] == 'd':
                            dflg = dflg + 1
                        if nodetype[ialpha2-1] == 'r':
                            rflg = rflg + 1
                        if nodetype[ialpha2-1] == 'd':
                            dflg = dflg + 1

                        if rflg != 0:
                            arofun[14] = arofun[14] + 1.0
                            mapfun[i-1][2][14] = 1.0
                            mapfun[iy-1][2][14] = 1.0
                            nfcr = nfcr + rflg
                        elif dflg != 0:
                            cdfun[14] = cdfun[14] + 1.0
                            mapfun[i-1][1][14] = 1.0
                            mapfun[iy-1][1][14] = 1.0
                            nfcd = nfcd + dflg
                        else:
                            alifun[14] = alifun[14] + 1.0
                            mapfun[i-1][0][14] = 1.0
                            mapfun[iy-1][0][14] = 1.0
                        ngrp = ngrp + 1

    if 'CHO' in chem:
        for i in range(1, node + 1):
            if group[i-1][:3] == 'CHO':
                nodmap(bond, i, node, 2, nnod, tnod)
                if nnod != 1:
                    print('-- error --, a unique C is expected in')
                    print('alpha position of a CHO  group')
                    print(chem)
                    raise Exception("in chemmap")
                ialpha = tnod[0]
                if nodetype[ialpha-1] == 'o':
                    ichecko = 0
                    ichecky = 0
                    for k in range(4):
                        if tabester[k][0] == ialpha:
                            ichecko = 1
                        if tabester[k][1] == i:
                            ichecky = 1
                    if ichecky == 1:
                        continue
                    if ichecko == 0:
                        continue

                if nodetype[ialpha-1] == 'r':
                    arofun[8] = arofun[8] + 1
                    mapfun[i-1][2][8] = 1
                    ngrp = ngrp + 1
                    nfcr = nfcr + 1
                elif nodetype[ialpha-1] == 'd':
                    cdfun[8] = cdfun[8] + 1
                    mapfun[i-1][1][8] = 1
                    ngrp = ngrp + 1
                    nfcd = nfcd + 1
                else:
                    alifun[8] = alifun[8] + 1
                    mapfun[i-1][0][8] = 1
                    ngrp = ngrp + 1

    if 'CO' in chem:
        for i in range(1, node + 1):
            if group[i-1][:3] == 'CO ':
                nodmap(bond, i, node, 2, nnod, tnod)
                if nnod != 2:
                    print('-- error --, only 2 C is expected in')
                    print('alpha position of a -CO-  group')
                    print(chem, '  nnod=', nnod)
                    print(bond[:node][:node])
                    raise Exception("in chemmap")
                rflg = 0
                dflg = 0
                for j in range(1, nnod + 1):
                    ialpha = tnod[j-1]
                    if nodetype[ialpha-1] == 'o':
                        ichecko = 0
                        ichecky = 0
                        for k in range(4):
                            if tabester[k][0] == ialpha:
                                ichecko = 1
                            if tabester[k][1] == i:
                                ichecky = 1
                        if ichecky == 1:
                            break
                        if ichecko == 0:
                            break
                    if nodetype[ialpha-1] == 'r':
                        rflg = rflg + 1
                    if nodetype[ialpha-1] == 'd':
                        dflg = dflg + 1
                else:
                    if rflg != 0:
                        arofun[9] = arofun[9] + 1.0
                        mapfun[i-1][2][9] = 1.0
                        nfcr = nfcr + rflg
                    elif dflg != 0:
                        cdfun[9] = cdfun[9] + 1.0
                        mapfun[i-1][1][9] = 1.0
                        nfcd = nfcd + dflg
                    else:
                        alifun[9] = alifun[9] + 1.0
                        mapfun[i-1][0][9] = 1.0
                    ngrp = ngrp + 1

    if 'CO(OONO2' in chem or 'CO(ONO2' in chem:
        for i in range(1, node + 1):
            if group[i-1][:9] == 'CO(OONO2)' or group[i-1][:8] == 'CO(ONO2)':
                nodmap(bond, i, node, 2, nnod, tnod)
                if nnod != 1:
                    print('-- error --, a unique C is expected in')
                    print('alpha position of a CO(OONO2)  group')
                    print(chem)
                    raise Exception("in chemmap")
                ialpha = tnod[0]
                if nodetype[ialpha-1] == 'o':
                    print('-- error --,  -O-CO(OONO2) group is unexpected')
                    print(chem)
                    nodetype[ialpha-1] = 'n'
                if nodetype[ialpha-1] == 'r':
                    arofun[12] = arofun[12] + 1
                    mapfun[i-1][2][12] = 1
                    ngrp = ngrp + 1
                    nfcr = nfcr + 1
                elif nodetype[ialpha-1] == 'd':
                    cdfun[12] = cdfun[12] + 1
                    mapfun[i-1][1][12] = 1
                    ngrp = ngrp + 1
                    nfcd = nfcd + 1
                else:
                    alifun[12] = alifun[12] + 1
                    mapfun[i-1][0][12] = 1
                    ngrp = ngrp + 1

    if '-O' in chem:
        for i in range(1, node + 1):
            if group[i-1][:3] == '-O-':
                nodmap(bond, i, node, 2, nnod, tnod)
                if nnod != 2:
                    print('-- error --, only 2 C is expected in')
                    print('alpha position of a -O-  group')
                    print(chem)
                    raise Exception("in chemmap")
                rflg = 0
                dflg = 0
                for j in range(1, nnod + 1):
                    ialpha = tnod[j-1]
                    if nodetype[ialpha-1] == 'y':
                        ichecko = 0
                        ichecky = 0
                        for k in range(4):
                            if tabester[k][0] == i:
                                ichecko = 1
                            if tabester[k][1] == ialpha:
                                ichecky = 1
                        if ichecko == 1:
                            break
                        if ichecky == 0:
                            break
                    if nodetype[ialpha-1] == 'r':
                        rflg = rflg + 1
                    if nodetype[ialpha-1] == 'd':
                        dflg = dflg + 1
                else:
                    if rflg != 0:
                        arofun[13] = arofun[13] + 1.0
                        mapfun[i-1][2][13] = 1.0
                        nfcr = nfcr + rflg

                        if rflg > 1:
                            ierr = 1
                            return
                    elif dflg != 0:
                        cdfun[13] = cdfun[13] + 1.0
                        mapfun[i-1][1][13] = 1.0
                        nfcd = nfcd + dflg
                    else:
                        alifun[13] = alifun[13] + 1.0
                        mapfun[i-1][0][13] = 1.0
                    ngrp = ngrp + 1

    for i in range(1, node + 1):
        for j in range(3):
            for k in range(20):
                if mapfun[i-1][j][k] != 0.0:
                    if mapfun[i-1][j][k] < 1.0:
                        funflg[i-1] = funflg[i-1] + 1
                    else:
                        funflg[i-1] = funflg[i-1] + int(mapfun[i-1][j][k])