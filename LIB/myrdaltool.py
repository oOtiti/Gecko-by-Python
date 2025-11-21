def myrdalprop(chem, bond, group, nring, rjg, weight, Tb, logPvap, deltaHvap):
    from ringtool import findring
    from toolbox import stoperr
    import math

    temp = 298.0
    tau = 0.0
    HBN = 0.0
    i = 0
    j = 0
    k = 0
    nca = 0
    sp2 = 0
    sp3 = 0
    indring = 0
    Jobgroup = [0] * 29
    OH = 0.0
    COOH = 0.0
    rgallpath = [[0] * len(group) for _ in range(len(rjg))]
    rgpath = [0] * len(group)
    nshare = 0
    begrg = 0
    endrg = 0
    rngflg = 0

    progname = 'myrdalprop '
    mesg = ''

    for i in range(len(group)):
        if group[i][0] != ' ':
            nca = nca + 1

    JdeltaTb = [
        23.58,  22.88,  21.74,  18.25,  18.18,  24.96,  24.14,  -0.03,
        38.13,  66.86,  92.88,  22.42,  76.75,  72.24, 169.09,  81.10,
        68.78, 112.10, 157.42,  27.15,  21.78,  21.32,  26.73,  31.01,
        31.22,  94.97,  76.34, 152.54, 115.30
    ]

    jobakgr(group, bond, nca, nring, rjg, Jobgroup)

    Tb = 0.0
    for i in range(29):
        Tb = Tb + JdeltaTb[i] * Jobgroup[i]
    Tb = Tb + 198.0

    sp3 = (Jobgroup[11] + Jobgroup[1] + Jobgroup[2] + Jobgroup[3] +
           Jobgroup[17] + Jobgroup[18] + Jobgroup[18] + Jobgroup[15] + Jobgroup[28])

    sp2 = (Jobgroup[13] + Jobgroup[14] + Jobgroup[12] + Jobgroup[5] + Jobgroup[15] +
           Jobgroup[17] + Jobgroup[18] + Jobgroup[18] + Jobgroup[6] + Jobgroup[27])

    indring = nring

    if nring > len(rjg):
        mesg = "Ring number > mri"
        stoperr(progname, mesg, chem)

    if nring == 2:
        for i in range(len(rgallpath)):
            for j in range(len(rgallpath[0])):
                rgallpath[i][j] = 0

        for i in range(nring):
            begrg = rjg[i][0]
            endrg = rjg[i][1]
            findring(begrg, endrg, nca, bond, rngflg, rgpath)
            for j in range(nca):
                if rgpath[j] == 1:
                    rgallpath[i][j] = 1

        indring = 1
        for i in range(nring - 1):
            for j in range(i + 1, nring):
                nshare = 0
                for k in range(nca):
                    if rgallpath[i][k] == 1 and rgallpath[j][k] == 1:
                        nshare = nshare + 1
                if nshare <= 1:
                    indring = indring + 1

    tau = sp3 + 0.5 * sp2 + 0.5 * indring - 1.0
    if tau < 0:
        tau = 0.0

    OH = Jobgroup[10] + Jobgroup[28]
    COOH = Jobgroup[14]
    HBN = (OH + COOH) ** 0.5 / weight

    logPvap = (-(86.0 + 0.4 * tau + 1421 * HBN) * (Tb - temp) / (19.1 * temp)
               + (-90.0 - 2.1 * tau) / 19.1 * ((Tb - temp) / temp - math.log(Tb / temp)))

    deltaHvap = (Tb * (86.0 + 0.4 * tau + 1421 * HBN)
                 + (-90.0 - 2.1 * tau) * (temp - Tb))

    return Tb, logPvap, deltaHvap


def jobakgr(group, bond, nca, nring, rjg, Jobgroup):
    from ringtool import findring

    i = 0
    j = 0
    k = 0
    begrg = 0
    endrg = 0
    ring = [0] * len(group)
    indexrg = [0] * len(group)
    rngflg = 0
    nc = 0
    ibeg = 0
    eflg = 0
    tbond = [row[:] for row in bond]

    for i in range(len(Jobgroup)):
        Jobgroup[i] = 0
    for i in range(len(indexrg)):
        indexrg[i] = 0

    if nring > len(rjg):
        print('--error-- in jobakgr. nring > mxring')
        raise Exception("in jobakgr")

    if nring > 0:
        for i in range(nring):
            begrg = rjg[i][0]
            endrg = rjg[i][1]
            findring(begrg, endrg, nca, bond, rngflg, ring)
            for j in range(nca):
                if ring[j] == 1:
                    indexrg[j] = 1

    for i in range(nca):
        if group[i][:4] == 'CHO ':
            Jobgroup[13] = Jobgroup[13] + 1
            continue
        elif group[i][:7] == 'CO(OH) ':
            Jobgroup[14] = Jobgroup[14] + 1
            continue
        elif group[i][:10] == 'CO(OONO2) ':
            Jobgroup[18] = Jobgroup[18] + 1
            continue
        elif group[i][:4] == '-O- ':
            if indexrg[i] == 0:
                eflg = 0
                for j in range(nca):
                    if tbond[i][j] == 3:
                        if group[j][:3] == 'CO ':
                            eflg = eflg + 1
                            for k in range(nca):
                                if tbond[j][k] == 3:
                                    tbond[j][k] = 0
                                    tbond[k][j] = 0
                if eflg >= 1:
                    Jobgroup[15] = Jobgroup[15] + 1
                    Jobgroup[12] = Jobgroup[12] - 1
                else:
                    Jobgroup[11] = Jobgroup[11] + 1
            else:
                Jobgroup[24] = Jobgroup[24] + 1
            continue
        elif group[i][:2] == 'CO':
            if indexrg[i] == 0:
                Jobgroup[12] = Jobgroup[12] + 1
            else:
                Jobgroup[25] = Jobgroup[25] + 1
            if group[i][:3] == 'CO ':
                continue
        elif group[i][:4] == 'CH3 ':
            Jobgroup[0] = Jobgroup[0] + 1
            continue
        elif group[i][:3] == 'CH2':
            if indexrg[i] == 0:
                Jobgroup[1] = Jobgroup[1] + 1
            else:
                Jobgroup[19] = Jobgroup[19] + 1
            if group[i][:4] == 'CH2 ':
                continue
        elif group[i][:5] == 'CdH2 ':
            Jobgroup[4] = Jobgroup[4] + 1
            continue
        elif group[i][:2] == 'CH':
            if indexrg[i] == 0:
                Jobgroup[2] = Jobgroup[2] + 1
            else:
                Jobgroup[20] = Jobgroup[20] + 1
            if group[i][:3] == 'CH ':
                continue
        elif group[i][:3] == 'CdH':
            if indexrg[i] == 0:
                Jobgroup[5] = Jobgroup[5] + 1
            else:
                Jobgroup[22] = Jobgroup[22] + 1
            if group[i][:4] == 'CdH ':
                continue
        elif group[i][:3] == 'Cd ':
            if indexrg[i] == 0:
                Jobgroup[6] = Jobgroup[6] + 1
            else:
                Jobgroup[23] = Jobgroup[23] + 1
        elif group[i][:1] == 'C':
            if indexrg[i] == 0:
                Jobgroup[3] = Jobgroup[3] + 1
            else:
                Jobgroup[21] = Jobgroup[21] + 1
            if group[i][:2] == 'C ':
                continue
        elif group[i][:2] == 'cH':
            Jobgroup[22] = Jobgroup[22] + 1
            if group[i][:3] == 'cH ':
                continue
        elif group[i][:1] == 'c':
            Jobgroup[23] = Jobgroup[23] + 1
            if group[i][:5] == 'c(OH)':
                Jobgroup[26] = Jobgroup[26] + 1
                continue

        nc = group[i].find(' ')
        if nc == -1:
            nc = len(group[i])

        ibeg = group[i].find('(OH)')
        if ibeg != -1:
            for j in range(ibeg, nc - 3):
                if group[i][j:j+4] == '(OH)' and group[i][:6] != 'CO(OH)':
                    Jobgroup[10] = Jobgroup[10] + 1

        ibeg = group[i].find('(OOH)')
        if ibeg != -1:
            for j in range(ibeg, nc - 4):
                if group[i][j:j+5] == '(OOH)':
                    Jobgroup[28] = Jobgroup[28] + 1

        ibeg = group[i].find('(OOOH)')
        if ibeg != -1:
            for j in range(ibeg, nc - 5):
                if group[i][j:j+6] == '(OOOH)':
                    Jobgroup[28] = Jobgroup[28] + 1

        ibeg = group[i].find('(ONO2)')
        if ibeg != -1:
            for j in range(ibeg, nc - 5):
                if group[i][j:j+6] == '(ONO2)':
                    Jobgroup[17] = Jobgroup[17] + 1

        ibeg = group[i].find('(NO2)')
        if ibeg != -1:
            for j in range(ibeg, nc - 4):
                if group[i][j:j+5] == '(NO2)':
                    Jobgroup[27] = Jobgroup[27] + 1