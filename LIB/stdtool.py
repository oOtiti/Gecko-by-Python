import keyparameter
import rjtool
import brchtool

# ================================================================
# PURPOSE : Set-up the tree of the C-C and C-O-C bond starting at    
#           the group in the input (top) and evaluate the longest    
#           path and its length.                                     
#    
#   The purpose of this routine is to evaluate the longest chain of  
#   the molecule. With the information in the bond-matrix, it is     
#   possible to evaluate the top-down trees of the molecule:         
#   e.g.                  CO(CH3)CH(CH(OH)(CH3))CHO                  
#                         1   2   3   4     5    6                   
#                   ___                                              
#                  | 1 |  parent of 2,3                              
#                   ---                                              
#                 /     \                                            
#             ___        ___                                         
#  child of 1| 2 |      | 3 | child of 1, parent of 4,6              
#             ---        ---                                         
#                      /     \                                       
#                   ___       ___                                    
#   child of 3 &   | 4 |     | 6 |  child of 3                       
#   parent of 5     ---       ---                                    
#                 /                                                  
#             ___                                                    
#            | 5 |  child of 4                                       
#             ---                                                    
#                                                                    
#   A child on the left hand is called "LEFT", on the right hand     
#   "RIGHT" and in the middle "CENTER".     
#   In CLNGTH the length of each longest tree is stored.  
#   Nevertheless the longest chain is the chain with most of the     
#   double bonds in it.                                              
#   First all relationships (parent - children) are evaluated. Then  
#   the top-down paths are defined so that first a LEFT, if not      
#   available, a RIGHT and then a CENTER is taken. As often as at    
#   least one child exists the specified path is followed further on.
#   In the last section the longest paths and the paths with the     
#   most double-bonds in it still remain and are given to the cal-   
#   ling routine.                                                    
#                                                                    
#********************************************************************
def lntree(bond, top, sec, nca, clngth, path):
    # 内部变量初始化
    bond_size = len(bond)
    left = [0] * (bond_size + 1)  # 1-based
    right = [0] * (bond_size + 1)
    center = [0] * (bond_size + 1)
    parent = [0] * (bond_size + 1)
    flag = [0] * (bond_size + 1)
    ptr = 0
    knt = 0
    nknt = 0
    nct1 = 0
    nct2 = 0
    maxlng = 0
    iend = 0
    i = 0
    j = 0
    k = 0
    
    # 初始化逻辑键矩阵
    tbond = [[False for _ in range(bond_size + 1)] for _ in range(bond_size + 1)]
    for i in range(1, bond_size + 1):
        for j in range(1, bond_size + 1):
            if bond[i-1][j-1] != 0:  # bond是0-based输入
                tbond[i][j] = True
    
    # 初始设置
    left[top] = sec
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
            
            # 寻找只有一个键且无父节点的基团
            if knt == 1:
                if parent[i] == 0 and i != top:
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
                        # 输出键矩阵错误
                        print('--error-- in lntree, ')
                        print('no path possible, error in bonding')
                        for j_row in range(1, nca + 1):
                            print('Row', j_row, ':', end=' ')
                            for j_col in range(1, nca + 1):
                                print(bond[j_row-1][j_col-1], end=' ')
                            print()
                        raise SystemExit("in lntree")
        # 无更多键则退出
        if nknt == 0:
            break
    
    # ---------------------------------------------
    # define all top-down pathes starting at "top"
    # ---------------------------------------------
    nct1 = nca - 1
    nct2 = nca + 4
    
    # 遍历所有可能路径
    for i in range(1, nct1 + 1):
        ptr = top
        path[top-1][i-1][0] = top  # path是3维0-based数组
        loop_continue = False
        for j in range(2, nct2 + 1):
            if flag[ptr] == 0:
                if left[ptr] != 0:
                    ptr = left[ptr]
                    path[top-1][i-1][j-1] = ptr
                else:
                    flag[ptr] = 1
            elif flag[ptr] == 1:
                if right[ptr] != 0:
                    ptr = right[ptr]
                    path[top-1][i-1][j-1] = ptr
                else:
                    flag[ptr] = 2
            elif flag[ptr] == 2:
                if center[ptr] != 0:
                    ptr = center[ptr]
                    path[top-1][i-1][j-1] = ptr
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
    
    # ---------------------
    # get the longest path 
    # ---------------------
    maxlng = 0
    for i in range(1, nca + 1):
        clngth[top-1][i-1] = 0  # clngth是2维0-based数组
        for j in range(1, nca + 1):
            if path[top-1][i-1][j-1] != 0:
                clngth[top-1][i-1] += 1
        
        if clngth[top-1][i-1] > maxlng:
            maxlng = clngth[top-1][i-1]
            iend = i - 1
            # 重置之前的路径长度为0
            for idx in range(iend):
                clngth[top-1][idx] = 0
        elif clngth[top-1][i-1] < maxlng:
            clngth[top-1][i-1] = 0

#=======================================================================
# PURPOSE : Check for misordered functionalities in parentheses next to 
# each other in the group given as input. The functionalities in the 
# group are sorted at the output For group priority see "pri" in module
# keyparameter where the priority invers TOP-DOWN.                                 
#                                                                     
#                       PTR21  PTR23   LNG3                           
#                         |      |____/___                            
#                C(..X1..)(..X2..)(..X3..)                            
#                 |      |        |      |                            
#               PTR11  PTR12    PTR32  PTR33                          
#=======================================================================
def ckgrppt(locat, group):
    # 初始化内部变量
    tgroup = ' ' * len(group)
    grp1 = ' ' * len(group)
    grp2 = ' ' * len(group)
    grp3 = ' ' * len(group)
    grp = ['   ' for _ in range(4)]  # 1-3索引使用
    ptr11 = 0
    ptr12 = 0
    ptr21 = 0
    ptr23 = 0
    ptr32 = 0
    ptr33 = 0
    cflg = 0
    i = 0
    p = 0
    lng3 = 0
    lengr = len(group)
    tri = False
    loswitch = False
    lopar = False
    
    # 设置初始指针（locat是1-based）
    ptr12 = locat
    ptr21 = locat + 1
    
    # ---------------------------------------
    # 寻找第一组括号的匹配
    # ---------------------------------------
    p = 1
    lopar = True
    for i in range(ptr12 - 1, 1, -1):
        if group[i-1:i] == ')':
            p += 1
        if group[i-1:i] == '(':
            p -= 1
        if p == 0:
            ptr11 = i
            lopar = False
            break
    if lopar:
        print('--error--, in ckgrppt. First group of parenthesis ')
        print(' mismatch for the group :', group.strip())
        raise SystemExit("in ckgrppt")
    
    # ---------------------------------------
    # 寻找第二组括号的匹配
    # ---------------------------------------
    p = 1
    lopar = True
    for i in range(ptr21 + 1, lengr + 1):
        if group[i-1:i] == '(':
            p += 1
        if group[i-1:i] == ')':
            p -= 1
        if p == 0:
            ptr23 = i
            lopar = False
            break
    if lopar:
        print('--error--, in ckgrppt. Second group of parenthesis ')
        print('mismatch for the group :', group.strip())
        raise SystemExit("in ckgrppt")
    
    # ---------------------------------------
    # 寻找第三组括号（如果存在）
    # ---------------------------------------
    ptr32 = ptr23 + 1
    lopar = True
    if ptr32 <= lengr and group[ptr32 - 1:ptr32] == '(':
        tri = True
        p = 1
        for i in range(ptr32 + 1, lengr + 1):
            if group[i-1:i] == ')':
                p -= 1
            if group[i-1:i] == '(':
                p += 1
            if p == 0:
                ptr33 = i
                lopar = False
                break
        if lopar:
            print('--error--, in ckgrppt. Third group of parenthesis ')
            print('mismatch for the group :', group.strip())
            raise SystemExit("in ckgrppt")
    elif ptr32 <= lengr and group[ptr32 - 1:ptr32] not in (' ', '.'):
        print('--error--, in ckgrppt. " " or "." is expected after')
        print('second group of parenthesis in :', group.strip())
        raise SystemExit("in ckgrppt")
    
    # -----------------
    # 定义所有基团
    # -----------------
    grp1 = group[ptr11 - 1:ptr12]
    grp[1] = group[ptr11:ptr11 + 3] if (ptr11 + 2) <= lengr else group[ptr11:]
    grp2 = group[ptr21 - 1:ptr23]
    grp[2] = group[ptr21:ptr21 + 3] if (ptr21 + 2) <= lengr else group[ptr21:]
    
    if tri:
        grp3 = group[ptr32 - 1:ptr33]
        grp[3] = group[ptr32:ptr32 + 3] if (ptr32 + 2) <= lengr else group[ptr32:]
        lng3 = ptr33 - ptr23
        
        # 检查第三组括号后是否为空格或点
        if ptr33 < lengr and group[ptr33:ptr33 + 1] != ' ':
            if ptr33 + 1 < lengr and group[ptr33:ptr33 + 2] != '. ':
                print('--error--, in ckgrppt. A " " is expected after')
                print('third group of parenthesis in :', group.strip())
                raise SystemExit("in ckgrppt")
    
    # 处理氟的特殊情况（仅1个字符）
    for i in range(1, 4):
        if grp[i][0] == 'F':
            grp[i] = 'F  '
    
    # --------------------------
    # 必要时交换基团顺序
    # --------------------------
    # 保存确定的前缀
    tgroup = group[:ptr11 - 1]
    
    # 交换前两组（如果顺序错误）
    if keyparameter.pri.find(grp[1]) < keyparameter.pri.find(grp[2]):
        loswitch = True
        ptr21_new = ptr11 + (ptr23 - ptr12)
        # 构建新字符串
        tgroup = tgroup.ljust(ptr11 - 1)
        tgroup += grp2
        tgroup = tgroup.ljust(ptr21_new - 1)
        tgroup += grp1
        if tri:
            tgroup = tgroup.ljust(ptr32 - 1)
            tgroup += group[ptr32 - 1:]
        # 截断到原长度
        if len(tgroup) > lengr:
            tgroup = tgroup[:lengr]
        group = tgroup
        ptr12 = ptr21_new - 1
    
    # 处理第三组（如果存在）
    if tri:
        tgroup = ' ' * lengr
        # 计算第三组的优先级
        cflg = 0
        for i in range(1, 3):
            if keyparameter.pri.find(grp[i]) < keyparameter.pri.find(grp[3]):
                cflg += 1
        
        if cflg == 1:
            # 第三组优先级介于前两组之间
            tgroup[:ptr12] = group[:ptr12]
            tgroup = tgroup[:ptr21] + grp3 + tgroup[ptr21 + len(grp3):]
            if loswitch:
                tgroup = tgroup[:ptr21 + lng3] + grp1 + tgroup[ptr21 + lng3 + len(grp1):]
            else:
                tgroup = tgroup[:ptr21 + lng3] + grp2 + tgroup[ptr21 + lng3 + len(grp2):]
            group = tgroup[:lengr]
        elif cflg == 2:
            # 第三组优先级最高
            tgroup[:ptr11 - 1] = group[:ptr11 - 1]
            tgroup = tgroup[:ptr11] + grp3 + tgroup[ptr11 + lng3:]
            tgroup = tgroup[:ptr11 + lng3] + group[ptr11:ptr23 + 1] + tgroup[ptr11 + lng3 + (ptr23 - ptr11 + 1):]
            group = tgroup[:lengr]

# ======================================================================
# PURPOSE: Check which formula in COPY has the highest priority,    
# comparing the position of different functional groups and return the
# "standardized" formula. 
# If the chemical is a radical, the formulas remain where the radical 
# group is at the end of the formula. If there are still more than one 
# writing left or the chemical is a non-radical molecule, according to 
# the group priorities, the formulas with the groups at the end remain, 
# respectively. This is done for all functional groups in the molecule 
# unless there are still more than one formulas left.                                         
# ======================================================================
def prioty(ogroup, rank, ocopy, ncp, nring, chem):
    # 初始化输出
    chem[0] = ' ' * len(chem[0])
    
    # 内部变量初始化
    i = 0
    j = 0
    ef = [0] * 6  # 1-5索引使用
    nelim = 0
    high = 0
    ncg = 0
    k = 0
    gcntr = 0
    eflag = 0
    ps = 0
    ptr = 0
    rjg = [[0 for _ in range(2)] for _ in range(keyparameter.mxring + 1)]
    rjs = [[0 for _ in range(2)] for _ in range(keyparameter.mxring + 1)]
    lofind = False
    temprjs = [[[0 for _ in range(2)] for __ in range(keyparameter.mxring + 1)] for ___ in range(len(ocopy) + 1)]
    minrjc = 0
    maxrjc = 0
    numgr = len(ogroup)
    prigr = ' ' * len(ogroup[0])
    group = [g.ljust(len(ogroup[0])) for g in ogroup]
    tempcopy = [' ' * len(chem[0]) for _ in range(len(ocopy) + 1)]
    copy = [' ' * len(chem[0]) for _ in range(len(ocopy) + 1)]
    
    # 复制输入
    for i in range(ncp):
        copy[i + 1] = ocopy[i]
    
    # --------------------------------------------------------
    # 检查自由基类型一致性
    # --------------------------------------------------------
    ef = [0] * 6
    for i in range(1, ncp + 1):
        if '.OO.' in copy[i]:
            ef[1] = 1
        elif 'CO(OO.)' in copy[i]:
            ef[2] = 1
        elif '(OO.)' in copy[i]:
            ef[3] = 1
        elif '(O.)' in copy[i]:
            ef[4] = 1
        else:
            ef[5] = 1
    eflag = sum(ef[1:6])
    if eflag > 1:
        print('--error--in prioty. Different radicals for molecule:')
        print(chem[0].strip())
        raise SystemExit("in prioty")
    
    # -----------------------------------------------------
    # 按优先级淘汰公式
    # -----------------------------------------------------
    nelim = 0
    if nring > 0:
        rjtool.rjgrm(nring, group, rjg)
    
    # 按基团优先级遍历
    for i in range(1, numgr + 1):
        # 找到当前优先级对应的基团
        prigr = ' ' * len(prigr)
        for j in range(1, numgr + 1):
            if rank[j - 1] == i:
                prigr = group[j - 1]
                break
        ncg = prigr.find(' ')
        if ncg == -1:
            ncg = len(prigr)
        if ncg < 1:
            continue
        
        # 统计该基团的数量
        gcntr = 0
        for j in range(1, numgr + 1):
            if group[j - 1] == prigr:
                gcntr += 1
        
        # 计算每个副本的优先级
        priort = [0] * (ncp + 1)
        for j in range(1, ncp + 1):
            priort[j] = 0
            ps = 0
            if copy[j][0] == ' ':
                continue
            
            # 移除环连接符
            if nring > 0:
                rjtool.rjsrm(nring, copy[j], rjs)
            
            # 找到每个基团的位置
            for k in range(1, gcntr + 1):
                lofind = False
                while True:
                    if lofind:
                        break
                    substr = copy[j][ps:]
                    ptr = substr.find(prigr[:ncg])
                    if ptr == -1:
                        print('--error--, from subroutine prioty')
                        print('Problem to detect ', prigr[:ncg])
                        print('in ', copy[j].strip())
                        raise SystemExit("in prioty")
                    ps += ptr
                    pos = ps + ncg - 1
                    pos_py = pos - 1  # 转换为0-based
                    
                    # 验证是否为完整基团
                    lofind = False
                    if pos_py + 1 < len(copy[j]):
                        next_char = copy[j][pos_py + 1]
                        if next_char == 'C':
                            lofind = True
                        elif (pos_py + 2 <= len(copy[j])) and (copy[j][pos_py + 1:pos_py + 3] == '(C'):
                            lofind = True
                        elif (pos_py + 2 <= len(copy[j])) and (copy[j][pos_py + 1:pos_py + 3] == ')C'):
                            lofind = True
                        elif (pos_py + 3 <= len(copy[j])) and (copy[j][pos_py + 1:pos_py + 4] == ')(C'):
                            lofind = True
                        elif next_char == 'c':
                            lofind = True
                        elif (pos_py + 2 <= len(copy[j])) and (copy[j][pos_py + 1:pos_py + 3] == '(c'):
                            lofind = True
                        elif (pos_py + 2 <= len(copy[j])) and (copy[j][pos_py + 1:pos_py + 3] == ')c'):
                            lofind = True
                        elif (pos_py + 3 <= len(copy[j])) and (copy[j][pos_py + 1:pos_py + 4] == ')(c'):
                            lofind = True
                        elif (pos_py + 2 <= len(copy[j])) and (copy[j][pos_py + 1:pos_py + 3] == '-O'):
                            lofind = True
                        elif (pos_py + 3 <= len(copy[j])) and (copy[j][pos_py + 1:pos_py + 4] == '(-O'):
                            lofind = True
                        elif (pos_py + 3 <= len(copy[j])) and (copy[j][pos_py + 1:pos_py + 4] == ')-O'):
                            lofind = True
                        elif (pos_py + 4 <= len(copy[j])) and (copy[j][pos_py + 1:pos_py + 5] == ')(-O'):
                            lofind = True
                    if pos_py < len(copy[j]) and copy[j][pos_py] == ' ':
                        lofind = True
                    
                    if lofind:
                        priort[j] += priort[j] + ps
            
            # 恢复环连接符
            if nring > 0:
                rjtool.rjsadd(nring, copy[j], rjs)
        
        # 找到最高优先级
        high = max(priort[1:ncp + 1])
        
        # 淘汰低优先级副本
        for j in range(1, ncp + 1):
            if priort[j] < high:
                if copy[j][0] != ' ':
                    copy[j] = ' ' * len(copy[j])
                    nelim += 1
                    if nelim == ncp - 1:
                        break
        if nelim == ncp - 1:
            break
    
    # -----------------------------------------------------
    # 处理环结构导致的剩余副本
    # -----------------------------------------------------
    if nelim != ncp - 1 and nring >= 2:
        # 保存环连接符位置
        for i in range(1, ncp + 1):
            tempcopy[i] = ' ' * len(tempcopy[i])
            for j in range(1, keyparameter.mxring + 1):
                for k in range(2):
                    temprjs[i][j][k] = 0
            if copy[i][0] != ' ':
                rjtool.rjsrm(nring, copy[i], rjs)
                tempcopy[i] = copy[i]
                for j in range(1, keyparameter.mxring + 1):
                    for k in range(2):
                        temprjs[i][j][k] = rjs[j - 1][k]  # rjs是0-based
                rjtool.rjsadd(nring, copy[i], rjs)
        
        # 验证副本除环连接符外是否一致
        for i in range(1, ncp):
            if tempcopy[i][0] != ' ' and tempcopy[i + 1][0] != ' ':
                if tempcopy[i] != tempcopy[i + 1]:
                    print('--error-- in prioty. Different formula in copy:')
                    print(copy[i].strip())
                    print(copy[i + 1].strip())
                    raise SystemExit("in prioty")
        
        # 根据环连接符位置选择副本
        for j in range(1, keyparameter.mxring + 1):
            # 找最小位置
            minrjc = len(chem[0])
            for i in range(1, ncp + 1):
                if copy[i][0] != ' ':
                    minrjc = min(temprjs[i][j][0], minrjc)
            # 淘汰非最小位置的副本
            for i in range(1, ncp + 1):
                if temprjs[i][j][0] != minrjc and copy[i][0] != ' ':
                    copy[i] = ' ' * len(copy[i])
                    nelim += 1
                    if nelim == ncp - 1:
                        break
            if nelim == ncp - 1:
                break
            # 找最大位置
            maxrjc = 0
            for i in range(1, ncp + 1):
                if copy[i][0] != ' ':
                    maxrjc = max(temprjs[i][j][1], maxrjc)
            # 淘汰非最大位置的副本
            for i in range(1, ncp + 1):
                if temprjs[i][j][1] != maxrjc and copy[i][0] != ' ':
                    copy[i] = ' ' * len(copy[i])
                    nelim += 1
                    if nelim == ncp - 1:
                        break
            if nelim == ncp - 1:
                break
    
    # ---------------------------------------------------------
    # 输出结果
    # ---------------------------------------------------------
    if nelim == ncp - 1:
        for i in range(1, ncp + 1):
            if copy[i][0] != ' ':
                chem[0] = copy[i]
                break
    else:
        print('--error-- in prioty. More than 1 copy left:')
        for j in range(1, ncp + 1):
            print(copy[j].strip())
        raise SystemExit("in prioty")
    
    if nring > 0:
        rjtool.rjgadd(nring, group, rjg)

# ===================================================================
# PURPOSE: make a copy of a molecule according to the longest tree.    
# On each call of mkcopy, the group ig is written in "copy". At the 
# same time and according to the bond-matrix, all groups having a bond 
# to ig (except on longest tree) are attached to it (with "()").      
#
# Ramifications are written according some priority rules. A max of 2  
# ramifications is expected and the formula must have the form :       
#        longest_tree-Cig(branch1)(branch2)-longest_tree               
# Since the branches may also contain ramification and functionalities,
# each branch must be written in a unique way. This is done by the        
# subroutine getbrch. If the carbon ig 'carry' 2 branches, then        
# branch1 and branch2 must be evaluated to know which branch has the   
# highest priority and must be written first.                          
#
# The structure of the subroutine is :                                 
#  1 - write the group ig and check if branching exist at pos. ig         
#  2 - if no branching => return                                       
#  3 - if only 1 branching => get and write the branch                 
#  4 - if 2 branching => get each branches, evaluate the priority of   
#       each branch and write them.    
#  It is expected that in most cases, the branch will only be 1C long. 
#  This case is tested first, since there is only 1 way to write the   
#  the branch.  
# ===================================================================
def mkcopy(lobond, group, nca, rank, nring, ig, pg, ng, ptr1, copy):
    # 内部变量初始化
    ptr2 = 0
    i = 0
    ialpha = 0
    ia1 = 0
    ia2 = 0
    ib1 = 0
    ib2 = 0
    ml1 = [0]
    ml2 = [0]
    maxpri = [0]
    tempbr = [' ' * len(copy) for _ in range(len(group) + 1)]
    brch1 = [' ' * len(copy)]
    brch2 = [' ' * len(copy)]
    
    # -----------------------------------------------------------
    # 写入当前基团ig
    # -----------------------------------------------------------
    group_ig = group[ig - 1]  # ig是1-based
    group_ig_len = group_ig.find(' ')
    if group_ig_len == -1:
        group_ig_len = len(group_ig)
    ptr2 = ptr1 + group_ig_len - 1
    if ptr2 < len(copy):
        copy = copy[:ptr1 - 1] + group_ig[:group_ig_len] + copy[ptr2:]
    ptr1 = ptr2 + 1
    
    # -----------------------------------------------------------
    # 寻找除pg和ng外的连接基团（alpha位置）
    # -----------------------------------------------------------
    ialpha = 0
    ia1 = 0
    ia2 = 0
    for i in range(1, nca + 1):
        if lobond[ig - 1][i - 1]:  # lobond是0-based
            if i != pg and i != ng:
                ialpha += 1
                if ialpha == 1:
                    ia1 = i
                if ialpha == 2:
                    ia2 = i
    
    # -----------------------------------------------------------
    # 无alpha基团，返回
    # -----------------------------------------------------------
    if ialpha == 0:
        return copy, ptr1
    
    # -----------------------------------------------------------
    # 一个alpha基团
    # -----------------------------------------------------------
    if ialpha == 1:
        # 检查beta位置是否有连接
        ib1 = 0
        for i in range(1, nca + 1):
            if lobond[ia1 - 1][i - 1]:
                if i != ig:
                    ib1 = i
                    break
        
        # 无beta位置，直接写入alpha基团
        if ib1 == 0:
            group_ia1 = group[ia1 - 1]
            group_ia1_len = group_ia1.find(' ')
            if group_ia1_len == -1:
                group_ia1_len = len(group_ia1)
            ptr2 = ptr1 + group_ia1_len
            if ptr2 < len(copy):
                copy = copy[:ptr1 - 1] + '(' + group_ia1[:group_ia1_len] + ')' + copy[ptr2:]
            ptr1 = ptr2 + 1
            return copy, ptr1
        else:
            # 获取分支并写入
            brchtool.getbrch(lobond, group, nca, rank, nring, ig, ia1, ib1, brch1, ml1)
            brch1_str = brch1[0].strip()
            brch1_len = len(brch1_str)
            ptr2 = ptr1 + brch1_len
            if ptr2 < len(copy):
                copy = copy[:ptr1 - 1] + '(' + brch1_str + ')' + copy[ptr2:]
            ptr1 = ptr2 + 1
            return copy, ptr1
    
    # -----------------------------------------------------------
    # 两个alpha基团
    # -----------------------------------------------------------
    if ialpha == 2:
        # 检查两个alpha基团的beta位置
        ib1 = 0
        for i in range(1, nca + 1):
            if lobond[ia1 - 1][i - 1]:
                if i != ig:
                    ib1 = i
                    break
        
        ib2 = 0
        for i in range(1, nca + 1):
            if lobond[ia2 - 1][i - 1]:
                if i != ig:
                    ib2 = i
                    break
        
        # 获取分支
        if ib1 != 0:
            brchtool.getbrch(lobond, group, nca, rank, nring, ig, ia1, ib1, brch1, ml1)
        if ib2 != 0:
            brchtool.getbrch(lobond, group, nca, rank, nring, ig, ia2, ib2, brch2, ml2)
        
        # Case 1: 两个都是甲基（无beta位置）
        if ib1 == 0 and ib2 == 0:
            if rank[ia2 - 1] < rank[ia1 - 1]:
                # 先写ia1
                group_ia1 = group[ia1 - 1]
                group_ia1_len = group_ia1.find(' ')
                if group_ia1_len == -1:
                    group_ia1_len = len(group_ia1)
                ptr2 = ptr1 + group_ia1_len
                if ptr2 < len(copy):
                    copy = copy[:ptr1 - 1] + '(' + group_ia1[:group_ia1_len] + ')' + copy[ptr2:]
                ptr1 = ptr2 + 1
                
                # 再写ia2
                group_ia2 = group[ia2 - 1]
                group_ia2_len = group_ia2.find(' ')
                if group_ia2_len == -1:
                    group_ia2_len = len(group_ia2)
                ptr2 = ptr1 + group_ia2_len
                if ptr2 < len(copy):
                    copy = copy[:ptr1 - 1] + '(' + group_ia2[:group_ia2_len] + ')' + copy[ptr2:]
                ptr1 = ptr2 + 1
            else:
                # 先写ia2
                group_ia2 = group[ia2 - 1]
                group_ia2_len = group_ia2.find(' ')
                if group_ia2_len == -1:
                    group_ia2_len = len(group_ia2)
                ptr2 = ptr1 + group_ia2_len
                if ptr2 < len(copy):
                    copy = copy[:ptr1 - 1] + '(' + group_ia2[:group_ia2_len] + ')' + copy[ptr2:]
                ptr1 = ptr2 + 1
                
                # 再写ia1
                group_ia1 = group[ia1 - 1]
                group_ia1_len = group_ia1.find(' ')
                if group_ia1_len == -1:
                    group_ia1_len = len(group_ia1)
                ptr2 = ptr1 + group_ia1_len
                if ptr2 < len(copy):
                    copy = copy[:ptr1 - 1] + '(' + group_ia1[:group_ia1_len] + ')' + copy[ptr2:]
                ptr1 = ptr2 + 1
            return copy, ptr1
        
        # Case 2: ia1是长链，ia2是甲基
        if ib1 != 0 and ib2 == 0:
            # 先写ia2（甲基）
            group_ia2 = group[ia2 - 1]
            group_ia2_len = group_ia2.find(' ')
            if group_ia2_len == -1:
                group_ia2_len = len(group_ia2)
            ptr2 = ptr1 + group_ia2_len
            if ptr2 < len(copy):
                copy = copy[:ptr1 - 1] + '(' + group_ia2[:group_ia2_len] + ')' + copy[ptr2:]
            ptr1 = ptr2 + 1
            
            # 再写ia1（长链）
            brch1_str = brch1[0].strip()
            brch1_len = len(brch1_str)
            ptr2 = ptr1 + brch1_len
            if ptr2 < len(copy):
                copy = copy[:ptr1 - 1] + '(' + brch1_str + ')' + copy[ptr2:]
            ptr1 = ptr2 + 1
            return copy, ptr1
        
        # Case 3: ia1是甲基，ia2是长链
        if ib1 == 0 and ib2 != 0:
            # 先写ia1（甲基）
            group_ia1 = group[ia1 - 1]
            group_ia1_len = group_ia1.find(' ')
            if group_ia1_len == -1:
                group_ia1_len = len(group_ia1)
            ptr2 = ptr1 + group_ia1_len
            if ptr2 < len(copy):
                copy = copy[:ptr1 - 1] + '(' + group_ia1[:group_ia1_len] + ')' + copy[ptr2:]
            ptr1 = ptr2 + 1
            
            # 再写ia2（长链）
            brch2_str = brch2[0].strip()
            brch2_len = len(brch2_str)
            ptr2 = ptr1 + brch2_len
            if ptr2 < len(copy):
                copy = copy[:ptr1 - 1] + '(' + brch2_str + ')' + copy[ptr2:]
            ptr1 = ptr2 + 1
            return copy, ptr1
        
        # Case 4: 两个都是长链
        if ib1 != 0 and ib2 != 0:
            if ml1[0] < ml2[0]:
                # 先写短链（ia1）
                brch1_str = brch1[0].strip()
                brch1_len = len(brch1_str)
                ptr2 = ptr1 + brch1_len
                if ptr2 < len(copy):
                    copy = copy[:ptr1 - 1] + '(' + brch1_str + ')' + copy[ptr2:]
                ptr1 = ptr2 + 1
                
                # 再写长链（ia2）
                brch2_str = brch2[0].strip()
                brch2_len = len(brch2_str)
                ptr2 = ptr1 + brch2_len
                if ptr2 < len(copy):
                    copy = copy[:ptr1 - 1] + '(' + brch2_str + ')' + copy[ptr2:]
                ptr1 = ptr2 + 1
                return copy, ptr1
            elif ml2[0] < ml1[0]:
                # 先写短链（ia2）
                brch2_str = brch2[0].strip()
                brch2_len = len(brch2_str)
                ptr2 = ptr1 + brch2_len
                if ptr2 < len(copy):
                    copy = copy[:ptr1 - 1] + '(' + brch2_str + ')' + copy[ptr2:]
                ptr1 = ptr2 + 1
                
                # 再写长链（ia1）
                brch1_str = brch1[0].strip()
                brch1_len = len(brch1_str)
                ptr2 = ptr1 + brch1_len
                if ptr2 < len(copy):
                    copy = copy[:ptr1 - 1] + '(' + brch1_str + ')' + copy[ptr2:]
                ptr1 = ptr2 + 1
                return copy, ptr1
            else:
                # 长度相同，比较优先级
                if brch1[0] == brch2[0]:
                    maxpri[0] = 1
                else:
                    tempbr[1] = brch1[0]
                    tempbr[2] = brch2[0]
                    for idx in range(3, len(tempbr)):
                        tempbr[idx] = ' '
                    brchtool.brpri(group, rank, tempbr, 2, nring, maxpri)
                
                if maxpri[0] == 1:
                    # 先写brch1
                    brch1_str = brch1[0].strip()
                    brch1_len = len(brch1_str)
                    ptr2 = ptr1 + brch1_len
                    if ptr2 < len(copy):
                        copy = copy[:ptr1 - 1] + '(' + brch1_str + ')' + copy[ptr2:]
                    ptr1 = ptr2 + 1
                    
                    # 再写brch2
                    brch2_str = brch2[0].strip()
                    brch2_len = len(brch2_str)
                    ptr2 = ptr1 + brch2_len
                    if ptr2 < len(copy):
                        copy = copy[:ptr1 - 1] + '(' + brch2_str + ')' + copy[ptr2:]
                    ptr1 = ptr2 + 1
                else:
                    # 先写brch2
                    brch2_str = brch2[0].strip()
                    brch2_len = len(brch2_str)
                    ptr2 = ptr1 + brch2_len
                    if ptr2 < len(copy):
                        copy = copy[:ptr1 - 1] + '(' + brch2_str + ')' + copy[ptr2:]
                    ptr1 = ptr2 + 1
                    
                    # 再写brch1
                    brch1_str = brch1[0].strip()
                    brch1_len = len(brch1_str)
                    ptr2 = ptr1 + brch1_len
                    if ptr2 < len(copy):
                        copy = copy[:ptr1 - 1] + '(' + brch1_str + ')' + copy[ptr2:]
                    ptr1 = ptr2 + 1
                return copy, ptr1
    
    # 超过两个alpha基团（错误）
    print('--error--, in mkcopy. More than 2 alpha groups')
    raise SystemExit("in mkcopy")

#========================================================
# PURPOSE :  Add "=" to formula with double bonds          
#========================================================
def dwrite(cpchem):
    # 初始化内部变量
    tempkc = ' ' * len(cpchem)
    nc = cpchem.find(' ')
    if nc == -1:
        nc = len(cpchem)
    else:
        nc -= 1
    idb = 0
    l = 0
    i = 0
    p = 0
    idp = 0
    ncd = 0
    j = 0
    ibflg = 0
    n = 0
    start = 0
    start2 = 0
    lengr = len(cpchem)
    
    # 遍历公式添加双键符号
    for i in range(1, nc + 1):
        # 跟踪括号嵌套
        if idp != 0:
            if cpchem[i-1:i] == '(':
                p += 1
            elif cpchem[i-1:i] == ')':
                p -= 1
            if p == 1:
                start2 = 1
        
        # 检测Cd基团
        if i + 1 <= nc and cpchem[i-1:i+1] == 'Cd':
            # 检查分支中是否有双键
            if ((i + 2 <= nc and cpchem[i+1:i+2] == '(') or (i + 3 <= nc and cpchem[i+2:i+3] == '(')) and idb == 0:
                ncd = 0
                ibflg = 0
                n = 0
                start = 0
                for j in range(i + 2, nc + 1):
                    if cpchem[j-1:j] == '(':
                        n += 1
                    elif cpchem[j-1:j] == ')':
                        n -= 1
                    if n == 1:
                        start = 1
                    if j + 1 <= nc and cpchem[j-1:j+1] == 'Cd' and ibflg == 0:
                        ncd += 1
                    if n == 0 and ibflg == 0 and start == 1:
                        ibflg = 1
                        if ncd == 2:
                            idp += 1
                            ncd = 0
                        else:
                            idb += 1
                            ncd = 0
            else:
                idb += 1
            
            # 添加等号
            if idp == 1 and p == 0 and start2 == 1:
                l += 1
                if l - 1 < len(tempkc):
                    tempkc = tempkc[:l-1] + '=' + tempkc[l:]
                idp = 0
                start2 = 0
                if idb == 1:
                    idb = 0
            elif idb > 1:
                l += 1
                if l - 1 < len(tempkc):
                    tempkc = tempkc[:l-1] + '=' + tempkc[l:]
                if i + 2 > nc or cpchem[i+1:i+3] != 'Cd':
                    idb = 0
        
        # 复制当前字符
        l += 1
        if l - 1 < len(tempkc):
            tempkc = tempkc[:l-1] + cpchem[i-1:i] + tempkc[l:]
    
    # 检查公式长度
    if l > lengr:
        print('--error-- in dwrite. Chemical formula is too long:')
        print(cpchem.strip())
        raise SystemExit("in dwrite")
    
    # 更新输入公式
    return tempkc

# ======================================================================
# PURPOSE: Write the reverse of chemical formula in COPY to output
#   The chain in input formula is written from right to left.     
# ======================================================================
def revers(copy, cc):
    # 初始化内部变量
    p = 0
    cc[0] = ' ' * len(cc[0])
    begin = 0  # 0-based
    iend = copy.find(' ')
    if iend == -1:
        iend = len(copy)
    else:
        iend -= 1
    lng = iend
    
    # 反转公式
    for i in range(lng, -1, -1):
        if copy[i:i+1] == ')':
            p += 1
        if copy[i:i+1] == '(':
            p -= 1
        if p == 0:
            # 检测C/c基团
            if copy[i:i+1] in ('C', 'c'):
                if i + 2 > len(copy) or copy[i:i+2] != 'Cl':
                    # 复制当前片段
                    fragment = copy[i:iend+1]
                    cc[0] = cc[0][:begin] + fragment + cc[0][begin + len(fragment):]
                    begin += len(fragment)
                    iend = i - 1
            # 检测-O-相关基团
            elif i + 3 <= len(copy) and copy[i:i+4] == '-O3-':
                cc[0] = cc[0][:begin] + '-O3-' + cc[0][begin + 4:]
                begin += 4
                iend = i - 1
            elif i + 3 <= len(copy) and copy[i:i+4] == '-O2-':
                cc[0] = cc[0][:begin] + '-O2-' + cc[0][begin + 4:]
                begin += 4
                iend = i - 1
            elif i + 3 <= len(copy) and copy[i:i+4] == '-O1-':
                cc[0] = cc[0][:begin] + '-O1-' + cc[0][begin + 4:]
                begin += 4
                iend = i - 1
            elif i + 2 <= len(copy) and copy[i:i+3] == '-O-':
                cc[0] = cc[0][:begin] + '-O-' + cc[0][begin + 3:]
                begin += 3
                iend = i - 1