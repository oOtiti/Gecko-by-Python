from keyparameter import mxring
from stdtool import lntree, ckgrppt, mkcopy, dwrite
from stdratings import ratings
from ringtool import uniqring

def rebond(bond, group, chem, nring):
    chem = ''
    mxca = len(bond)
    nca = 0
    idbflg = 0
    last = 0
    
    for n in range(1, mxca + 1):
        if group[n-1][0] != ' ':
            nca += 1
            last = n
        if 'Cd' in group[n-1]:
            idbflg = 1
    
    tgroup = group.copy()
    tbond = [row[:] for row in bond]
    
    for i in range(1, last + 1):
        locat = tgroup[i-1].find(')(')
        if locat != -1:
            ckgrppt(locat, tgroup[i-1])
    
    if nca < 1:
        return chem, nring
    
    if nca == 1:
        for i in range(1, last + 1):
            if tgroup[i-1] != ' ':
                chem = tgroup[i-1]
                return chem, nring
    
    ncx = 0
    for i in range(1, last + 1):
        for j in range(i + 1, last + 1):
            if tbond[i-1][j-1] > 0:
                ncx += 1
    nring = ncx - nca + 1
    
    icheck = 0
    i = 1
    while i <= nca:
        if i > nca:
            break
        if tgroup[i-1] == ' ':
            for j in range(i, last):
                tgroup[j-1] = tgroup[j]
            tgroup[last-1] = ' '
            for j in range(1, last + 1):
                for k in range(i, last):
                    tbond[k-1][j-1] = tbond[k][j-1]
            for j in range(1, last + 1):
                for k in range(i, last):
                    tbond[j-1][k-1] = tbond[j-1][k]
        
        if tgroup[i-1] != ' ':
            i += 1
        
        icheck += 1
        if icheck > last:
            print('--error-- in rebond. Infinite loop when erasing blanks lines')
            raise Exception("in rebond")
    
    rank = [0] * len(tbond)
    ratings(nca, tgroup, tbond, nring, rank)
    
    if nring > 0:
        rjg = [[0, 0] for _ in range(mxring)]
        uniqring(nring, nca, tgroup, tbond, rank, rjg)
        for k in range(nring):
            i = rjg[k][0]
            j = rjg[k][1]
            tbond[i-1][j-1] = 0
            tbond[j-1][i-1] = 0
    
    lobond = [[False for _ in range(len(tbond[0]))] for _ in range(len(tbond))]
    for i in range(len(tbond)):
        for j in range(len(tbond[0])):
            if tbond[i][j] != 0:
                lobond[i][j] = True
    clngth = [[0 for _ in range(len(tbond))] for _ in range(len(tbond))]
    path = [[[0 for _ in range(len(tbond))] for _ in range(len(tbond))] for _ in range(len(tbond))]
    beg = 0
    for i in range(1, last + 1):
        if tgroup[i-1] != ' ':
            beg = i
            for j in range(i + 1, mxca + 1):
                if lobond[i-1][j-1]:
                    clngth = [[0 for _ in range(len(tbond))] for _ in range(len(tbond))]
                    path = [[[0 for _ in range(len(tbond))] for _ in range(len(tbond))] for _ in range(len(tbond))]
                    lntree(tbond, i, j, nca, clngth, path)
                    break
            break
    
    for i in range(1, nca + 1):
        if clngth[beg-1][i-1] != 0:
            leaf = path[beg-1][i-1][clngth[beg-1][i-1] - 1]
            last2 = path[beg-1][i-1][clngth[beg-1][i-1] - 2]
            clngth2 = [[0 for _ in range(len(tbond))] for _ in range(len(tbond))]
            path2 = [[[0 for _ in range(len(tbond))] for _ in range(len(tbond))] for _ in range(len(tbond))]
            lntree(tbond, leaf, last2, nca, clngth2, path2)
            for j in range(1, nca + 1):
                if clngth2[leaf-1][j-1] != 0:
                    ptr = 1
                    for k in range(1, nca + 1):
                        ig = path2[leaf-1][j-1][k-1]
                        if ig != 0:
                            pg = 0
                            ng = 0
                            if k > 1:
                                pg = path2[leaf-1][j-1][k-2]
                            if k < nca:
                                ng = path2[leaf-1][j-1][k]
                            mkcopy(lobond, tgroup, nca, rank, nring, ig, pg, ng, ptr, chem)
                    break
            break
    
    if idbflg != 0:
        dwrite(chem)
    
    return chem, nring

def swap(gold, pold, pnew, gnew):
    lengr = len(gold)
    lengr2 = len(gnew)
    if lengr != lengr2:
        print("stop in swap - wrong size")
        raise Exception("in swap")
    
    gnew_list = [''] * lengr
    lg = len(gold.strip())
    lold = len(pold.strip())
    lnew = len(pnew.strip())
    
    if lg + lnew - lold > lengr:
        print('--error-- in swap')
        print('new group more than lengr characters')
        raise Exception('in swap')
    
    ibeg = gold.find(pold.strip())
    if ibeg == -1:
        print('--error-- in swap')
        print(f'{pold} not in {gold}')
        raise Exception("in swap")
    
    iend = ibeg + lold - 1
    
    n = 0
    for i in range(lengr):
        if n >= lengr:
            break
        if i < ibeg or i > iend:
            if n < lengr:
                gnew_list[n] = gold[i]
                n += 1
        if i == ibeg:
            for j in range(lnew):
                if n < lengr:
                    gnew_list[n] = pnew[j]
                    n += 1
    
    gnew = ''.join(gnew_list)
    return gnew