def ratings(nca, ogroup, bond, nring, rank):
    from keyparameter import digit, mxring
    from rjtool import rjgrm, rjgadd
    from primetool import primes

    group = ogroup[:]
    ncx = [0] * len(ogroup)
    rjg = [[0] * 2 for _ in range(mxring)]
    nrat = 24
    nines = '999999999999999999999999'
    maxm = ''
    rating = [''] * len(ogroup)

    for i in range(len(rank)):
        rank[i] = 0

    for i in range(1, nca + 1):
        ncx[i-1] = 0
        rating[i-1] = '0' + nines[1:]
        for j in range(1, nca + 1):
            if bond[i-1, j-1] > 0:
                ncx[i-1] = ncx[i-1] + 1
                rating[i-1] = digit(ncx[i-1]) + rating[i-1][1:]

    if nring > 0:
        rjgrm(nring, group, rjg)

    for i in range(1, nca + 1):
        rating[i-1] = rating[i-1][0] + nines[1:]
        if '.' in group[i-1]:
            rating[i-1] = rating[i-1][:1] + '1' + rating[i-1][2:]
        if 'CdO' in group[i-1]:
            rating[i-1] = rating[i-1][:3] + '1' + rating[i-1][4:]
        if 'Cd' in group[i-1]:
            rating[i-1] = rating[i-1][:4] + '3' + rating[i-1][5:]
        if 'CdH' in group[i-1]:
            rating[i-1] = rating[i-1][:4] + '2' + rating[i-1][5:]
        if 'CdH2' in group[i-1]:
            rating[i-1] = rating[i-1][:4] + '1' + rating[i-1][5:]
        if '-O-' in group[i-1]:
            rating[i-1] = rating[i-1][:5] + '1' + rating[i-1][6:]
        if 'CHO' in group[i-1]:
            rating[i-1] = rating[i-1][:6] + '1' + rating[i-1][7:]
        if 'CO(OONO2)' in group[i-1]:
            rating[i-1] = rating[i-1][:7] + '1' + rating[i-1][8:]
        if 'CO(OOH)' in group[i-1]:
            rating[i-1] = rating[i-1][:8] + '1' + rating[i-1][9:]
        if '(OOH)' in group[i-1]:
            rating[i-1] = rating[i-1][:9] + '3' + rating[i-1][10:]
        if '(OOH)(OOH)' in group[i-1]:
            rating[i-1] = rating[i-1][:9] + '2' + rating[i-1][10:]
        if '(OOH)(OOH)(OOH)' in group[i-1]:
            rating[i-1] = rating[i-1][:9] + '1' + rating[i-1][10:]
        if '(ONO2)' in group[i-1]:
            rating[i-1] = rating[i-1][:10] + '3' + rating[i-1][11:]
        if '(ONO2)(ONO2)' in group[i-1]:
            rating[i-1] = rating[i-1][:10] + '2' + rating[i-1][11:]
        if '(ONO2)(ONO2)(ONO2)' in group[i-1]:
            rating[i-1] = rating[i-1][:10] + '1' + rating[i-1][11:]
        if 'CO(OH)' in group[i-1]:
            rating[i-1] = rating[i-1][:11] + '1' + rating[i-1][12:]
        if 'CO' in group[i-1]:
            rating[i-1] = rating[i-1][:12] + '1' + rating[i-1][13:]
        if 'F' in group[i-1]:
            rating[i-1] = rating[i-1][:13] + '1' + rating[i-1][14:]
        if 'Br' in group[i-1]:
            rating[i-1] = rating[i-1][:14] + '1' + rating[i-1][15:]
        if 'Cl' in group[i-1]:
            rating[i-1] = rating[i-1][:15] + '1' + rating[i-1][16:]
        if 'S' in group[i-1]:
            rating[i-1] = rating[i-1][:16] + '1' + rating[i-1][17:]
        if 'NH' in group[i-1]:
            rating[i-1] = rating[i-1][:17] + '1' + rating[i-1][18:]
        if '(NO2)' in group[i-1]:
            rating[i-1] = rating[i-1][:18] + '3' + rating[i-1][19:]
        if '(NO2)(NO2)' in group[i-1]:
            rating[i-1] = rating[i-1][:18] + '2' + rating[i-1][19:]
        if '(NO2)(NO2)(NO2)' in group[i-1]:
            rating[i-1] = rating[i-1][:18] + '1' + rating[i-1][19:]
        if 'NO' in group[i-1]:
            rating[i-1] = rating[i-1][:19] + '1' + rating[i-1][20:]
        if '(OH)' in group[i-1]:
            rating[i-1] = rating[i-1][:20] + '3' + rating[i-1][21:]
        if '(OH)(OH)' in group[i-1]:
            rating[i-1] = rating[i-1][:20] + '2' + rating[i-1][21:]
        if '(OH)(OH)(OH)' in group[i-1]:
            rating[i-1] = rating[i-1][:20] + '1' + rating[i-1][21:]
        if 'C' in group[i-1]:
            rating[i-1] = rating[i-1][:21] + '4' + rating[i-1][22:]
        if 'CH' in group[i-1]:
            rating[i-1] = rating[i-1][:21] + '3' + rating[i-1][22:]
        if 'CH2' in group[i-1]:
            rating[i-1] = rating[i-1][:21] + '2' + rating[i-1][22:]
        if 'CH3' in group[i-1]:
            rating[i-1] = rating[i-1][:21] + '1' + rating[i-1][22:]
        if 'c' in group[i-1]:
            rating[i-1] = rating[i-1][:22] + '2' + rating[i-1][23:]
        if 'cH' in group[i-1]:
            rating[i-1] = rating[i-1][:22] + '1' + rating[i-1][23:]

    cntr = 0
    nines = '999999999999999999999999'
    grloop = True
    i = 1
    while grloop and i <= nca:
        maxm = nines
        for j in range(1, nca + 1):
            if rating[j-1] < maxm:
                maxm = rating[j-1]
        if maxm == nines:
            grloop = False
            break
        cntr = cntr + 1
        for j in range(1, nca + 1):
            if rating[j-1] <= maxm:
                rating[j-1] = nines
                rank[j-1] = cntr
        i += 1

    for i in range(1, nca + 1):
        if rank[i-1] == 0:
            print('--error--, in ratings, rank not assigned for group:', i)
            raise Exception("in ratings")

    if nring > 0:
        rjgadd(nring, group, rjg)

    primes(nca, bond, rank)