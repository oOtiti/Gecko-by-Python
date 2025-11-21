# ======================================================================
# Purpose: return the (non-standardised) chemical formula "chem" 
# for the bond and group tables provided has input
# ======================================================================
from keyparameter import mxring
from stdtool import lntree, ckgrppt, mkcopy, dwrite
from stdratings import ratings
from ringtool import uniqring

def rebond(bond, group):
    # 初始化输出变量
    chem = ' '
    nring = 0
    
    # 声明并初始化内部变量
    mxca = len(bond)  # SIZE(bond,1)
    n = len(bond[0]) if mxca > 0 else 0  # SIZE(bond,2)
    last = 0
    icheck = 0
    nca = 0
    idbflg = 0
    ncx = 0
    beg = 0
    leaf = 0
    last2 = 0
    ptr = 0
    ig = 0
    pg = 0
    ng = 0
    locat = 0
    
    # 数组初始化（对应Fortran声明）
    rjg = [[0 for _ in range(2)] for _ in range(mxring)]  # rjg(mxring,2)
    rank = [0 for _ in range(mxca)]  # rank(SIZE(bond,1))
    # path(SIZE(bond,1),SIZE(bond,1),SIZE(bond,1))
    path = [[[0 for _ in range(mxca)] for _ in range(mxca)] for _ in range(mxca)]
    # clngth(SIZE(bond,1),SIZE(bond,1))
    clngth = [[0 for _ in range(mxca)] for _ in range(mxca)]
    # lobond(SIZE(bond,1),SIZE(bond,2))
    lobond = [[False for _ in range(n)] for _ in range(mxca)]
    # tgroup(SIZE(group)) - 拷贝输入group，保持1-based逻辑（Python用索引-1访问）
    tgroup = [' '] * (len(group) + 1)  # tgroup(1~len(group))
    for i in range(len(group)):
        tgroup[i+1] = group[i]
    # tbond(SIZE(bond,1),SIZE(bond,2)) - 拷贝输入bond，1-based
    tbond = [[0 for _ in range(n)] for _ in range(mxca + 1)]  # tbond(1~mxca, 1~n)
    for i in range(mxca):
        for j in range(n):
            tbond[i+1][j+1] = bond[i][j]
    
    # 初始化部分变量
    for i in range(mxca):
        for j in range(mxca):
            clngth[i][j] = 0
    for i in range(mxca):
        for j in range(n):
            lobond[i][j] = False
    for i in range(mxca):
        for j in range(mxca):
            for k in range(mxca):
                path[i][j][k] = 0
    
    # 计算nca（非空group数量）和idbflg（是否含Cd组）
    nca = 0
    idbflg = 0
    for n_idx in range(1, mxca + 1):  # Fortran: n=1,mxca
        if tgroup[n_idx][:1] != ' ':
            nca += 1
            last = n_idx
        if 'Cd' in tgroup[n_idx]:
            idbflg = 1
    
    # 检查每个group中官能团排序（是否有')('）
    for i in range(1, last + 1):  # Fortran: i=1,last
        locat = tgroup[i].find(')(')
        if locat != -1:
            # locat在Fortran中是1-based，Python中find返回0-based，需+1匹配原逻辑
            ckgrppt(locat + 1, tgroup[i])
    if nca < 1:
        return chem.strip(), nring
    
    # 若只有一个非空group，直接返回
    if nca == 1:
        for i in range(1, last + 1):  # Fortran: i=1,last
            if tgroup[i] != ' ':
                chem = tgroup[i]
                return chem.strip(), nring
    
    # 计算环数量nring
    ncx = 0
    for i in range(1, last + 1):  # Fortran: i=1,last
        for j in range(i + 1, last + 1):  # Fortran: j=i+1,last
            if tbond[i][j] > 0:
                ncx += 1
    nring = ncx - nca + 1
    
    # 删除bond和group中的空行（rmloop循环）
    icheck = 0
    i = 1
    while True:  # Fortran: rmloop DO
        if i > nca:
            break
        if tgroup[i] == ' ':
            # 移动group：j从i到last，tgroup(j)=tgroup(j+1)
            for j in range(i, last):  # Fortran: j=i,last
                tgroup[j] = tgroup[j + 1]
            tgroup[last] = ' '
            # 移动bond列：j从1到last，tbond(k,j)=tbond(k+1,j)（k从i到last-1）
            for j in range(1, last + 1):  # Fortran: j=1,last
                for k in range(i, last):  # Fortran: k=i,last-1
                    tbond[k][j] = tbond[k + 1][j]
            # 移动bond行：j从1到last，tbond(j,k)=tbond(j,k+1)（k从i到last-1）
            for j in range(1, last + 1):  # Fortran: j=1,last
                for k in range(i, last):  # Fortran: k=i,last-1
                    tbond[j][k] = tbond[j][k + 1]
        if tgroup[i] != ' ':
            i += 1
        
        icheck += 1
        if icheck > last:
            print('--error-- in rebond. Infinite loop when erasing blanks lines')
            raise SystemExit("in rebond")
    
    # 计算基团优先级，处理环结构
    ratings(nca, tgroup, tbond, nring, rank)
    if nring > 0:
        uniqring(nring, nca, tgroup, tbond, rank, rjg)
        for k in range(1, nring + 1):  # Fortran: k=1,nring
            i = rjg[k - 1][0]  # rjg是Python list，0-based
            j = rjg[k - 1][1]
            tbond[i][j] = 0
            tbond[j][i] = 0
    
    # 构建逻辑键矩阵lobond（tbond≠0时为True）
    for i in range(1, mxca + 1):  # Fortran: 遍历bond所有行
        for j in range(1, n + 1):  # Fortran: 遍历bond所有列
            if tbond[i][j] != 0:
                lobond[i - 1][j - 1] = True  # lobond是Python 0-based数组
    
    # 找到键的起始点beg，构建路径树
    beg = 0
    grloop_break = False
    for i in range(1, last + 1):  # Fortran: i=1,last（grloop DO）
        if tgroup[i] != ' ':
            beg = i
            for j in range(i + 1, mxca + 1):  # Fortran: j=i+1,mxca
                # lobond是0-based，i-1,j-1对应tbond(i,j)
                if lobond[i - 1][j - 1]:
                    lntree(tbond, i, j, nca, clngth, path)
                    grloop_break = True
                    break
            if grloop_break:
                break
    
    # 构建化学 formula（wrtloop循环）
    wrtloop_break = False
    for i in range(1, nca + 1):  # Fortran: i=1,nca（wrtloop DO）
        # clngth是0-based？原Fortran中clngth(beg,i)，beg和i是1-based，需调整
        if clngth[beg - 1][i - 1] != 0:
            # path(beg,i,clngth(beg,i)) -> 转换为Python 0-based索引
            clng_val = clngth[beg - 1][i - 1]
            leaf = path[beg - 1][i - 1][clng_val - 1]
            last2 = path[beg - 1][i - 1][clng_val - 2] if clng_val >= 2 else 0
            lntree(tbond, leaf, last2, nca, clngth, path)
            
            for j in range(1, nca + 1):  # Fortran: j=1,nca
                if clngth[leaf - 1][j - 1] != 0:
                    ptr = 1
                    for k in range(1, nca + 1):  # Fortran: k=1,nca
                        ig = path[leaf - 1][j - 1][k - 1]
                        if ig != 0:
                            if k > 1:
                                pg = path[leaf - 1][j - 1][k - 2]
                            else:
                                pg = 0
                            if k < nca:
                                ng = path[leaf - 1][j - 1][k]
                            else:
                                ng = 0
                            mkcopy(lobond, tgroup, nca, rank, nring, ig, pg, ng, ptr, chem)
                    wrtloop_break = True
                    break
            if wrtloop_break:
                break
    
    # 处理双键（添加'='）
    if idbflg != 0:
        dwrite(chem)
    
    return chem.strip(), nring


# ======================================================================
# Purpose : swaps pieces pold, pnew in group (gold, gnew)
# ======================================================================
def swap(gold, pold, pnew):
    # 初始化变量
    lengr = len(gold)
    lengr2 = len(pnew)  # gnew的长度由返回值决定，这里先获取pnew长度
    gnew = ' ' * lengr  # 初始化gnew为对应长度的空格字符串
    
    # 检查gold和gnew长度一致性（原Fortran逻辑）
    if lengr != len(gnew):
        print("stop in swap - wrong size")
        raise SystemExit("in swap")
    
    # 计算有效长度（排除末尾空格）
    lg = gold.find(' ')
    if lg == -1:
        lg = lengr
    else:
        lg -= 1  # Fortran: INDEX(gold,' ')-1
    
    lold = pold.find(' ')
    if lold == -1:
        lold = len(pold)
    else:
        lold -= 1  # Fortran: INDEX(pold,' ')-1
    
    lnew = pnew.find(' ')
    if lnew == -1:
        lnew = len(pnew)
    else:
        lnew -= 1  # Fortran: INDEX(pnew,' ')-1
    
    # 检查新基团长度是否超过限制
    if lg + lnew - lold > lengr:
        print('--error-- in swap')
        print('new group more than lengr characters')
        raise SystemExit('in swap')
    
    # 定位pold在gold中的位置（ibeg：1-based）
    ibeg = gold.find(pold[:lold]) + 1  # Fortran是1-based索引，find返回0-based+1
    if ibeg == 0:  # 未找到（find返回-1，+1后为0）
        print('--error-- in swap')
        print(f"{pold} not in {gold}")
        raise SystemExit("in swap")
    iend = ibeg + lold - 1  # Fortran: ibeg+lold-1
    
    # 构建新基团gnew
    n = 0
    for i in range(1, lengr + 1):  # Fortran: i=1,lengr（charloop DO）
        if n >= lengr:
            break
        # 不在pold位置范围内，直接复制gold字符
        if i < ibeg or i > iend:
            n += 1
            # Python字符串不可变，先转列表处理
            gnew_list = list(gnew)
            gnew_list[n - 1] = gold[i - 1]  # i是1-based，gold[i-1]是对应字符
            gnew = ''.join(gnew_list)
        # 在pold起始位置，替换为pnew
        if i == ibeg:
            for j in range(1, lnew + 1):  # Fortran: j=1,lnew
                n += 1
                if n - 1 >= lengr:
                    break
                gnew_list = list(gnew)
                gnew_list[n - 1] = pnew[j - 1]  # pnew[j-1]是对应字符
                gnew = ''.join(gnew_list)
    
    return gnew.strip()