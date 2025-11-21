def uniqring(nring, nca, group, bond, rank, rjg):
    from keyparameter import prim
    from rjtool import rjgrm, rjgadd
    from stdratings import ratings

    ppb = [[0] * len(bond) for _ in range(len(bond))]
    cprim = [0] * len(group)
    tbond = [[0] * len(bond) for _ in range(len(bond))]
    i = 0
    j = 0
    ii = 0
    jj = 0
    k = 0
    maxp = 0
    n = 0
    rngflg = 0
    top = 0
    sec = 0
    parent = [0] * len(group)
    child = [[0] * 3 for _ in range(len(group))]
    nk = [0] * len(group)
    iold = [0] * len(group)
    trjg = [[0] * 2 for _ in range(nring)]
    trank = [0] * len(group)
    tgroup = [''] * len(group)
    rgpath = [0] * len(bond)

    for i in range(len(rjg)):
        for j in range(len(rjg[0])):
            rjg[i][j] = 0
    for i in range(len(tbond)):
        for j in range(len(tbond[0])):
            tbond[i][j] = 0
    for i in range(len(tgroup)):
        tgroup[i] = ''

    for i in range(nca):
        cprim[i] = prim[rank[i]]
    for i in range(len(ppb)):
        for j in range(len(ppb[0])):
            ppb[i][j] = 0

    for i in range(1, nca + 1):
        for j in range(1, nca + 1):
            if i != j and bond[i-1][j-1] > 0:
                ppb[i-1][j-1] = cprim[i-1] * cprim[j-1]
                tbond[i-1][j-1] = bond[i-1][j-1]

    rjgrm(nring, group, rjg)

    for n in range(1, nring + 1):
        grloop = True
        while grloop:
            maxp = 0
            for i in range(1, nca + 1):
                for j in range(i + 1, nca + 1):
                    if bond[i-1][j-1] != 2:
                        if ppb[i-1][j-1] > maxp:
                            ii = i
                            jj = j
                            maxp = ppb[ii-1][jj-1]
            findring(ii, jj, nca, tbond, rngflg, rgpath)
            ppb[ii-1][jj-1] = 0
            if rngflg == 1:
                tbond[ii-1][jj-1] = 0
                tbond[jj-1][ii-1] = 0
                rjg[n-1][0] = ii
                rjg[n-1][1] = jj
                for i in range(len(tgroup)):
                    tgroup[i] = group[i]
                if nring > 1:
                    ratings(nca, tgroup, tbond, nring, rank)
                for i in range(len(tgroup)):
                    tgroup[i] = ''
            else:
                continue
            grloop = False

    top = ii
    sec = 0
    nodloop = False
    for j in range(2, nca + 1):
        for i in range(1, j):
            if tbond[i-1][j-1] > 0:
                if i == ii:
                    sec = j
                if j == ii:
                    sec = i
                if sec > 0:
                    nodloop = True
                    break
        if nodloop:
            break

    findtree(tbond, top, sec, nca, parent, child)

    for i in range(len(nk)):
        nk[i] = 0
    j = 1
    iold[j-1] = top
    if parent[iold[j-1]-1] != 0:
        print("error 1 in uniqring")
        raise Exception("in uniqring")

    i = j
    xxloop = True
    while xxloop:
        k = nk[i-1] + 1
        if k <= 3:
            if child[iold[i-1]-1][k-1] > 0:
                nk[i-1] = nk[i-1] + 1
                j = j + 1
                iold[j-1] = child[iold[i-1]-1][k-1]
                i = j
                continue
        if parent[iold[i-1]-1] == 0:
            xxloop = False
            break
        for ii in range(1, nca + 1):
            if parent[iold[i-1]-1] == iold[ii-1]:
                i = ii
                break
        else:
            xxloop = False

    for ii in range(1, nca + 1):
        if iold[ii-1] != 0:
            trank[ii-1] = rank[iold[ii-1]-1]
            tgroup[ii-1] = group[iold[ii-1]-1]
        for jj in range(1, nca + 1):
            if iold[jj-1] != 0:
                tbond[ii-1][jj-1] = bond[iold[ii-1]-1][iold[jj-1]-1]
            else:
                tbond[ii-1][jj-1] = 0
        for n in range(1, nring + 1):
            for i in range(1, 3):
                if rjg[n-1][i-1] == iold[ii-1]:
                    trjg[n-1][i-1] = ii

    for ii in range(1, nca + 1):
        rank[ii-1] = trank[ii-1]
        group[ii-1] = tgroup[ii-1]
        for jj in range(1, nca + 1):
            bond[ii-1][jj-1] = tbond[ii-1][jj-1]
    for n in range(1, nring + 1):
        rjg[n-1][0] = min(trjg[n-1][0], trjg[n-1][1])
        rjg[n-1][1] = max(trjg[n-1][0], trjg[n-1][1])

    rjgadd(nring, group, rjg)


def findring(i1, i2, nca, bond, rngflg, ring):
    from keyparameter import mxcp
    from mapping import gettrack

    track = [[0] * len(bond) for _ in range(mxcp)]
    trlen = [0] * mxcp
    ntr = 0
    i = 0
    j = 0
    m = 0
    n = 0

    rngflg = 0
    for i in range(len(ring)):
        ring[i] = 0

    for m in range(1, nca):
        for n in range(m + 1, nca + 1):
            if bond[m-1][n-1] != 0:
                gettrack(bond, m, nca, ntr, track, trlen)
                for i in range(1, ntr + 1):
                    for j in range(3, trlen[i-1] + 1):
                        if track[i-1][j-1] == n:
                            if (m == i1 and n == i2) or (m == i2 and n == i1):
                                rngflg = 1
                            for k in range(1, j + 1):
                                ring[track[i-1][k-1]-1] = 1


def findtree(con, top, sec, nca, parent, child):
    left = [0] * len(con)
    right = [0] * len(con)
    center = [0] * len(con)
    ptr = 0
    knt = 0
    nknt = 0
    i = 0
    j = 0
    k = 0
    iend = 0
    tcon = [[False] * len(con) for _ in range(len(con))]

    ptr = 0
    for i in range(len(left)):
        left[i] = 0
        right[i] = 0
        center[i] = 0
    for i in range(len(parent)):
        parent[i] = 0

    for i in range(len(con)):
        for j in range(len(con)):
            if con[i][j] != 0:
                tcon[i][j] = True

    left[top-1] = sec
    parent[sec-1] = top
    tcon[top-1][sec-1] = False
    tcon[sec-1][top-1] = False

    for k in range(1, nca + 1):
        nknt = 0
        for i in range(1, nca + 1):
            knt = 0
            for j in range(1, nca + 1):
                if tcon[i-1][j-1]:
                    knt = knt + 1
                    ptr = j
            nknt = nknt + knt
            if knt == 1:
                if parent[i-1] == 0 and i != top:
                    parent[i-1] = ptr
                    tcon[i-1][ptr-1] = False
                    tcon[ptr-1][i-1] = False
                    if left[ptr-1] == 0:
                        left[ptr-1] = i
                    elif right[ptr-1] == 0:
                        right[ptr-1] = i
                    elif center[ptr-1] == 0:
                        center[ptr-1] = i
                    else:
                        for j in range(1, nca + 1):
                            line = ''
                            for iend in range(1, nca + 1):
                                line += f" {con[j-1][iend-1]}"
                            print(line)
                        print('--error-- in findtree. No possible path')
                        raise Exception("in findtree")
        if nknt == 0:
            break

    for i in range(1, nca + 1):
        child[i-1][0] = left[i-1]
        child[i-1][1] = right[i-1]
        child[i-1][2] = center[i-1]


def rejoin(nca, x, y, x1, y1, bond, group):
    tbond = [[0] * len(bond) for _ in range(len(bond))]
    tgroup = [''] * len(group)
    i = 0
    j = 0
    ii = 0
    jj = 0
    icut = 0

    for i in range(len(tbond)):
        for j in range(len(tbond[0])):
            tbond[i][j] = 0
    for i in range(len(tgroup)):
        tgroup[i] = ''

    if x > y:
        icut = x - 1
    else:
        icut = y - 1
    x1 = x - icut
    if x1 <= 0:
        x1 = x1 + nca
    y1 = y - icut
    if y1 <= 0:
        y1 = y1 + nca

    for i in range(1, icut + 1):
        ii = nca - icut + i
        tgroup[ii-1] = group[i-1]
        for j in range(1, nca + 1):
            tbond[ii-1][j-1] = bond[i-1][j-1]

    for i in range(icut + 1, nca + 1):
        ii = i - icut
        tgroup[ii-1] = group[i-1]
        for j in range(1, nca + 1):
            tbond[ii-1][j-1] = bond[i-1][j-1]

    for i in range(len(bond)):
        for j in range(len(bond[0])):
            bond[i][j] = tbond[i][j]

    for i in range(1, nca + 1):
        for j in range(1, icut + 1):
            jj = nca - icut + j
            tbond[i-1][jj-1] = bond[i-1][j-1]
        for j in range(icut + 1, nca + 1):
            jj = j - icut
            tbond[i-1][jj-1] = bond[i-1][j-1]

    for i in range(len(group)):
        group[i] = tgroup[i]
    for i in range(len(bond)):
        for j in range(len(bond[0])):
            bond[i][j] = tbond[i][j]


def ring_data(ig, nca, tbond, tgroup, ndrg, lodrg, trackrg):
    from keyparameter import mxcp
    from mapping import gettrack

    i = 0
    j = 0
    k = 0
    l = 0
    rngflg = 0
    ring = [0] * len(tgroup)
    track = [[0] * len(tgroup) for _ in range(mxcp)]
    trlen = [0] * mxcp
    ntr = 0
    ring_tmp = [False] * len(tgroup)

    ndrg = 0
    for i in range(len(lodrg)):
        for j in range(len(lodrg[0])):
            lodrg[i][j] = False
    for i in range(len(trackrg)):
        for j in range(len(trackrg[0])):
            trackrg[i][j] = 0

    for i in range(1, nca + 1):
        if tbond[ig-1][i-1] != 0:
            findring(ig, i, nca, tbond, rngflg, ring)
            gettrack(tbond, ig, nca, ntr, track, trlen)
            if rngflg == 1:
                for j in range(1, ntr + 1):
                    for k in range(len(ring_tmp)):
                        ring_tmp[k] = False
                    for k in range(3, trlen[j-1] + 1):
                        if track[j-1][k-1] == i:
                            if ndrg == 0:
                                ndrg = 1
                                for l in range(1, k + 1):
                                    lodrg[0][track[j-1][l-1]-1] = True
                                    trackrg[0][l-1] = track[j-1][l-1]
                            else:
                                for l in range(1, k + 1):
                                    ring_tmp[track[j-1][l-1]-1] = True
                                for l in range(ndrg):
                                    match = True
                                    for m in range(len(ring_tmp)):
                                        if lodrg[l][m] != ring_tmp[m]:
                                            match = False
                                            break
                                    if match:
                                        break
                                else:
                                    ndrg = ndrg + 1
                                    if ndrg > len(lodrg):
                                        print("in ring_data, ndrg > SIZE(lodrg,1)")
                                        raise Exception("in ring_data")
                                    for l in range(1, k + 1):
                                        trackrg[ndrg-1][l-1] = track[j-1][l-1]
                                        lodrg[ndrg-1][track[j-1][l-1]-1] = True