import keyparameter
import rjtool

# ======================================================================
# PURPOSE:  write the branch starting at the position ia (alpha       
# position) and bounded to the position ig. One of the beta positions 
# is known (ib). The various ways of writing the branch are checked 
#  => A standard "formula" of the chain is given as output                                        
#                                                                     
# The structure of the subroutine is :                                
# 1- find longest trees starting from ia                              
# 2- write for each tree a copy of the corresponding formula          
# 3- find which formula has the highest priority                      
# ======================================================================
def getbrch(olobond, group, nca, rank, nring, ig, ia, ib, brch, ml):
    # 初始化输出参数（brch为列表传递，ml为列表传递）
    brch[0] = ' ' * len(brch[0])
    ml[0] = 0
    
    # 内部变量初始化
    np = 0
    tnp = 0
    maxpri = 0
    ptra = 0
    iadd = 0
    pre = 0
    nex = 0
    # 初始化brpath：维度与olobond一致（1-based）
    brpath_rows = len(olobond)
    brpath_cols = len(olobond[0]) if brpath_rows > 0 else nca + 2
    brpath = [[0 for _ in range(brpath_cols)] for _ in range(brpath_rows)]
    # 初始化brcopy和tbrcopy
    brcopy_len = len(brch[0])
    brcopy = [' ' * brcopy_len for _ in range(len(group) + 1)]  # 1-based
    tbrcopy = [' ' * brcopy_len for _ in range(len(group) + 1)]
    i = 0
    j = 0
    
    # 复制逻辑键矩阵（深拷贝）
    lobond = [[val for val in row] for row in olobond]
    
    # search the longest tree, starting from ia
    # -----------------------------------------
    lobond[ia][ig] = False
    lobond[ig][ia] = False
    treebr(lobond, ig, ia, ib, nca, brpath, ml, np)
    lobond[ia][ig] = True
    lobond[ig][ia] = True
    
    # write copy of longest branches
    # ------------------------------
    for i in range(1, np + 1):
        ptra = 1
        brcopy[i] = ' ' * brcopy_len
        for j in range(1, ml[0] + 1):
            iadd = brpath[i][j]
            if j == 1:
                pre = ig
            else:
                pre = brpath[i][j - 1]
            nex = brpath[i][j + 1] if (j + 1) < brpath_cols else 0
            mkbrcopy(lobond, group, nca, rank, iadd, pre, nex, ptra, brcopy[i])
    
    # If only 1 chain : write out and return
    # ---------------------------------------
    if np == 1:
        brch[0] = brcopy[1]
        return
    
    # if more than 1 chain : collapse identical copy
    # --------------------------------------------------------------
    for i in range(1, np + 1):
        for j in range(i + 1, np + 1):
            if brcopy[i] == brcopy[j]:
                brcopy[j] = ' ' * brcopy_len
    
    tnp = np
    np = 0
    for i in range(1, tnp + 1):
        if brcopy[i][0] != ' ':
            np += 1
            tbrcopy[np] = brcopy[i]
    
    # if only 1 chain remain then write out and return else search priority
    if np == 1:
        brch[0] = tbrcopy[1]
        return
    else:
        brpri(group, rank, tbrcopy, np, nring, maxpri)
        brch[0] = tbrcopy[maxpri[0]]
        return

# ======================================================================
# PURPOSE: THIS ROUTINE IS VERY SIMILAR TO MKCOPY. Main difference is 
# that mkcopy is called to write copies of the "full" molecule while 
# mkbrcopy is called to write copies of a given branch only.                           
# mkbrcopy makes a copy of a branch according to the longest tree. On 
# each call of mkcopy, the group ig is written in "copy". At the same 
# time and according to the bond-matrix, all groups which have a bond 
# to ig (except the bond of longest tree), are attached to it (in 
# parentheses). These attached must not have more than 1C 
# ======================================================================
def mkbrcopy(lobond, group, nca, rank, ig, pg, ng, ptr1, brcopy):
    ptr2 = 0
    i = 0
    j = 0
    ita = 0
    ia1 = 0
    ia2 = 0
    
    # ------------------------------------------------------
    # WRITE GROUP IG TO BRCOPY
    # ------------------------------------------------------
    group_ig = group[ig - 1] if (ig - 1) < len(group) else ''  # group可能是0-based列表
    group_ig_len = group_ig.find(' ')
    if group_ig_len == -1:
        group_ig_len = len(group_ig)
    ptr2 = ptr1 + group_ig_len - 1
    if ptr2 <= len(brcopy):
        brcopy = brcopy[:ptr1 - 1] + group_ig[:group_ig_len] + brcopy[ptr2:]
    ptr1 = ptr2 + 1
    
    # ------------------------------------------------------
    # LOOK FOR ALL ATTACHED GROUPS TO GROUP IG, NOT PG OR NG
    # ------------------------------------------------------
    ita = 0
    ia1 = 0
    ia2 = 0
    for i in range(1, nca + 1):
        if lobond[ig][i]:
            if (i != pg) and (i != ng):
                ita += 1
                if ita == 1:
                    ia1 = i
                if ita == 2:
                    ia2 = i
                
                # 检查附着基团是否只有一个键（到ig）
                for j in range(1, nca + 1):
                    if lobond[i][j]:
                        if j != ig:
                            print('--error--,in mkbrcopy. The branch has a group in beta ')
                            print('position of the closest C in the main (longest) branch')
                            raise SystemExit("in mkbrcopy")
    
    # ------------------------------------------------------
    # CHECK THE VARIOUS POSSIBLE CASES
    # ------------------------------------------------------
    # CASE 1 (no alpha group)
    if ita == 0:
        return brcopy
    
    # CASE 2 (1 alpha group)
    if ita == 1:
        group_ia1 = group[ia1 - 1] if (ia1 - 1) < len(group) else ''
        group_ia1_len = group_ia1.find(' ')
        if group_ia1_len == -1:
            group_ia1_len = len(group_ia1)
        ptr2 = ptr1 + group_ia1_len
        if ptr2 <= len(brcopy):
            brcopy = brcopy[:ptr1 - 1] + '(' + group_ia1[:group_ia1_len] + ')' + brcopy[ptr2:]
        ptr1 = ptr2 + 1
        return brcopy
    
    # CASE 3 (2 alpha groups)
    if ita == 2:
        if rank[ia2 - 1] < rank[ia1 - 1]:  # rank是0-based列表
            # 先写ia1
            group_ia1 = group[ia1 - 1] if (ia1 - 1) < len(group) else ''
            group_ia1_len = group_ia1.find(' ')
            if group_ia1_len == -1:
                group_ia1_len = len(group_ia1)
            ptr2 = ptr1 + group_ia1_len
            if ptr2 <= len(brcopy):
                brcopy = brcopy[:ptr1 - 1] + '(' + group_ia1[:group_ia1_len] + ')' + brcopy[ptr2:]
            ptr1 = ptr2 + 1
            
            # 再写ia2
            group_ia2 = group[ia2 - 1] if (ia2 - 1) < len(group) else ''
            group_ia2_len = group_ia2.find(' ')
            if group_ia2_len == -1:
                group_ia2_len = len(group_ia2)
            ptr2 = ptr1 + group_ia2_len
            if ptr2 <= len(brcopy):
                brcopy = brcopy[:ptr1 - 1] + '(' + group_ia2[:group_ia2_len] + ')' + brcopy[ptr2:]
            ptr1 = ptr2 + 1
            return brcopy
        else:
            # 先写ia2
            group_ia2 = group[ia2 - 1] if (ia2 - 1) < len(group) else ''
            group_ia2_len = group_ia2.find(' ')
            if group_ia2_len == -1:
                group_ia2_len = len(group_ia2)
            ptr2 = ptr1 + group_ia2_len
            if ptr2 <= len(brcopy):
                brcopy = brcopy[:ptr1 - 1] + '(' + group_ia2[:group_ia2_len] + ')' + brcopy[ptr2:]
            ptr1 = ptr2 + 1
            
            # 再写ia1
            group_ia1 = group[ia1 - 1] if (ia1 - 1) < len(group) else ''
            group_ia1_len = group_ia1.find(' ')
            if group_ia1_len == -1:
                group_ia1_len = len(group_ia1)
            ptr2 = ptr1 + group_ia1_len
            if ptr2 <= len(brcopy):
                brcopy = brcopy[:ptr1 - 1] + '(' + group_ia1[:group_ia1_len] + ')' + brcopy[ptr2:]
            ptr1 = ptr2 + 1
            return brcopy
    
    # more than 2 groups: error
    if ita > 2:
        print('--error--,in mkbrcopy. # of methyl group > 2')
        raise SystemExit("in mkbrcopy")

#=======================================================================
# PURPOSE: This routine is very similar to prioty. Main difference is 
# that prioty is called to find the standardized formula of the molecule
# while brpri is called to find which branch has the highest priority.
# The subroutine check which formula of the branch in "copy" has the 
# highest priority and return the corresponding index in "copy"                                                  
#
#              -SEE PRIOTY FOR ADDITIONAL COMMENT -                   
#=======================================================================
def brpri(ogroup, rank, copy_arr, ncp, nring, maxpri):
    # 初始化输出参数
    maxpri[0] = 0
    
    # 内部变量初始化
    group = [g.ljust(len(ogroup[0])) for g in ogroup]  # 复制group（保持长度）
    prigr = ' ' * len(ogroup[0])
    i = 0
    j = 0
    nelim = 0
    high = 0
    priort = [0 for _ in range(keyparameter.mxcp + 1)]  # 1-based
    ncg = 0
    k = 0
    gcntr = 0
    ps = 0
    ptr = 0
    rjg = [[0 for _ in range(2)] for _ in range(keyparameter.mxring + 1)]
    rjs = [[0 for _ in range(2)] for _ in range(keyparameter.mxring + 1)]
    lofind = False
    ngr = len(ogroup)
    
    # -----------------------------------------------------
    # INITIALIZE
    # -----------------------------------------------------
    for idx in range(len(priort)):
        priort[idx] = 0
    
    # -----------------------------------------------------
    # FIND PRIORITY AND ELIMINATE SUCCESSIVELY THE FORMULA
    # -----------------------------------------------------
    nelim = 0
    
    # loop over the groups, starting with the highest priority group
    if nring > 0:
        rjtool.rjgrm(nring, group, rjg)
    
    # 替换aloop标签循环为标准for循环
    for i in range(1, ngr + 1):
        # 找到rank为i的基团
        prigr = ' ' * len(prigr)
        for j in range(1, ngr + 1):
            if rank[j - 1] == i:  # rank是0-based
                prigr = group[j - 1]
                break
        ncg = prigr.find(' ')
        if ncg == -1:
            ncg = len(prigr)
        if ncg < 1:
            continue
        
        # count the number of groups identical to prigr
        gcntr = 0
        for j in range(1, ngr + 1):
            if group[j - 1] == prigr:
                gcntr += 1
        
        # loop over the copies of the branch and set priority
        for j in range(1, ncp + 1):
            priort[j] = 0
            if copy_arr[j].find(prigr[:ncg]) == -1:
                continue
            ps = 0
            
            # remove ring-join characters if present
            if nring > 0:
                rjtool.rjsrm(nring, copy_arr[j], rjs)
            
            # find the position of the group in the branch
            for k in range(1, gcntr + 1):
                lofind = False
                while True:
                    if lofind:
                        break
                    substr = copy_arr[j][ps:]
                    ptr = substr.find(prigr[:ncg])
                    if ptr == -1:
                        break
                    ps += ptr
                    pos = ps + ncg - 1  # 1-based对应copy_arr的索引（Python是0-based，需调整）
                    pos_py = pos - 1
                    
                    # 检查是否为完整基团
                    lofind = False
                    if pos_py + 1 < len(copy_arr[j]):
                        next_char = copy_arr[j][pos_py + 1]
                        if next_char == 'C':
                            lofind = True
                        elif (pos_py + 2 <= len(copy_arr[j])) and (copy_arr[j][pos_py + 1:pos_py + 3] == '(C'):
                            lofind = True
                        elif (pos_py + 2 <= len(copy_arr[j])) and (copy_arr[j][pos_py + 1:pos_py + 3] == ')C'):
                            lofind = True
                        elif (pos_py + 3 <= len(copy_arr[j])) and (copy_arr[j][pos_py + 1:pos_py + 4] == ')(C'):
                            lofind = True
                        elif next_char == 'c':
                            lofind = True
                        elif (pos_py + 2 <= len(copy_arr[j])) and (copy_arr[j][pos_py + 1:pos_py + 3] == '(c'):
                            lofind = True
                        elif (pos_py + 2 <= len(copy_arr[j])) and (copy_arr[j][pos_py + 1:pos_py + 3] == ')c'):
                            lofind = True
                        elif (pos_py + 3 <= len(copy_arr[j])) and (copy_arr[j][pos_py + 1:pos_py + 4] == ')(c'):
                            lofind = True
                        elif (pos_py + 3 <= len(copy_arr[j])) and (copy_arr[j][pos_py + 1:pos_py + 4] == ')-O'):
                            lofind = True
                        elif (pos_py + 2 <= len(copy_arr[j])) and (copy_arr[j][pos_py + 1:pos_py + 3] == '-O'):
                            lofind = True
                    if pos_py < len(copy_arr[j]) and copy_arr[j][pos_py] == ' ':
                        lofind = True
                    
                    if lofind:
                        priort[j] += priort[j] + ps  # 原逻辑的优先级计算
            # end grloop
            
            # replace ring-join characters if necessary
            if nring > 0:
                rjtool.rjsadd(nring, copy_arr[j], rjs)
        # end cploop
        
        # find the maximum priority
        high = 0
        for j in range(1, ncp + 1):
            if priort[j] > high:
                high = priort[j]
        
        # eliminate copy with priority < maximum
        for j in range(1, ncp + 1):
            if priort[j] < high:
                if copy_arr[j][0] != ' ':
                    copy_arr[j] = ' ' * len(copy_arr[j])
                    nelim += 1
                    if nelim == ncp - 1:
                        break  # 替代exit aloop
        # end eliminate
        if nelim == ncp - 1:
            break
    # end aloop
    
    # ---------------------------------------------------------
    # check that only one copy remains and return its ID number
    # ---------------------------------------------------------
    if nelim == ncp - 1:
        for i in range(1, ncp + 1):
            if copy_arr[i][0] != ' ':
                maxpri[0] = i
                break
    else:
        print('--error-- in brpri. More than one copy left:')
        for j in range(1, ncp + 1):
            print(copy_arr[j])
        raise SystemExit("in brpri")
    
    if nring > 0:
        rjtool.rjgadd(nring, group, rjg)

# ======================================================================
# PURPOSE : THIS ROUTINE IS VERY SIMILAR TO LNTREE. Main difference  
# is that lntree is called to find the longest tree in the molecule 
# while treebr is called to find the longest tree in a given branch of 
# the molecule. Treebr sets up the tree of the C-C bond starting at the  
# group given in the input (top) and evaluate the longest path (brpath), 
# its length (maxlng) and the total # of path having the maximum length.                
#
# SEE LNTREE FOR ADDITIONAL COMMENT                 
# ======================================================================
def treebr(lobond, ig, top, sec, nca, brpath, maxlng, npath):
    # 初始化输出参数
    maxlng[0] = 0
    npath[0] = 0
    
    # 内部变量初始化
    lobond_size = len(lobond) if len(lobond) > 0 else nca + 2
    left = [0 for _ in range(lobond_size + 1)]  # 1-based
    right = [0 for _ in range(lobond_size + 1)]
    center = [0 for _ in range(lobond_size + 1)]
    parent = [0 for _ in range(lobond_size + 1)]
    tlngth = [0 for _ in range(lobond_size + 1)]
    tpath = [[0 for _ in range(lobond_size + 4)] for _ in range(lobond_size + 1)]  # 扩展列数
    flag = [0 for _ in range(lobond_size + 1)]
    ptr = 0
    knt = 0
    nknt = 0
    nct1 = 0
    nct2 = 0
    iend = 0
    i = 0
    j = 0
    k = 0
    
    # 复制键矩阵
    tbond = [[val for val in row] for row in lobond]
    
    # -----------
    # initialize
    # -----------
    for idx in range(len(left)):
        left[idx] = 0
        right[idx] = 0
        center[idx] = 0
        parent[idx] = 0
        flag[idx] = 0
    for idx in range(len(tlngth)):
        tlngth[idx] = 0
    for row in range(len(brpath)):
        for col in range(len(brpath[row])):
            brpath[row][col] = 0
    for row in range(len(tpath)):
        for col in range(len(tpath[row])):
            tpath[row][col] = 0
    
    left[top] = sec
    parent[top] = ig
    parent[sec] = top
    tbond[top][sec] = False
    tbond[sec][top] = False
    
    # ---------------------------------------
    # get the relationships (parent/children)
    # ---------------------------------------
    for k in range(1, nca + 1):
        nknt = 0
        for i in range(1, nca + 1):
            knt = 0
            ptr = 0
            for j in range(1, nca + 1):
                if tbond[i][j]:
                    knt += 1
                    ptr = j
            nknt += knt
            
            # look for ith carbon with only one node, no parent
            if knt == 1:
                if parent[i] == 0 and i != top and i != ig:
                    parent[i] = ptr
                    tbond[i][ptr] = False
                    tbond[ptr][i] = False
                    if left[ptr] == 0:
                        left[ptr] = i
                    elif right[ptr] == 0:
                        right[ptr] = i
                    elif center[ptr] == 0:
                        center[ptr] = i
                    else:
                        # 输出键矩阵错误信息
                        print('--error-- in  treebr. No path, bonding pb')
                        for j_row in range(1, nca + 1):
                            print('Row', j_row, ':', end=' ')
                            for j_col in range(1, nca + 1):
                                print(lobond[j_row][j_col], end=' ')
                            print()
                        raise SystemExit("in treebr")
        # 无更多键则退出
        if nknt == 0:
            break
    
    # ---------------------------------------------
    # define all top-down paths starting at "top"
    # ---------------------------------------------
    nct1 = nca - 1
    nct2 = nca + 4
    
    # 替换brloop标签循环为标准for循环+条件判断
    for i in range(1, nct1 + 1):
        ptr = top
        tpath[i][1] = top
        loop_continue = False
        for j in range(2, nct2 + 1):
            if flag[ptr] == 0:
                if left[ptr] != 0:
                    ptr = left[ptr]
                    tpath[i][j] = ptr
                else:
                    flag[ptr] = 1
            elif flag[ptr] == 1:
                if right[ptr] != 0:
                    ptr = right[ptr]
                    tpath[i][j] = ptr
                else:
                    flag[ptr] = 2
            elif flag[ptr] == 2:
                if center[ptr] != 0:
                    ptr = center[ptr]
                    tpath[i][j] = ptr
                else:
                    flag[ptr] = 3
                    flag[parent[ptr]] += 1
                    loop_continue = True
                    break
            elif flag[ptr] == 3:
                flag[parent[ptr]] += 1
                loop_continue = True
                break
        if loop_continue:
            continue
    # end brloop
    
    # ---------------------
    # get the longest path 
    # ---------------------
    # 计算每条路径的长度，找到最大长度
    maxlng[0] = 0
    for i in range(1, nca + 1):
        tlngth[i] = 0
        for j in range(1, nca + 1):
            if tpath[i][j] != 0:
                tlngth[i] += 1
        if tlngth[i] > maxlng[0]:
            maxlng[0] = tlngth[i]
    
    # 找到所有最长路径
    npath[0] = 0
    for i in range(1, nca + 1):
        if tlngth[i] == maxlng[0]:
            npath[0] += 1
            for j in range(1, nca + 1):
                if j < len(brpath[npath[0]]):
                    brpath[npath[0]][j] = tpath[i][j]
    
    # 检查路径有效性
    if npath[0] == 0:
        print('--error-- in treebr. No path found')
        raise SystemExit("in treebr")
    if maxlng[0] < 2:
        print('--error-- in treebr. Length of the path shorter than expected')
        raise SystemExit("in treebr")
    if maxlng[0] > nca - 3:
        print('--error-- in treebr. Length of the path greater than expected')
        raise SystemExit("in treebr")