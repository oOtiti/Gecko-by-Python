# 前置库已由用户提供，此处仅按需求导入对应模块内容
from keyparameter import mxcp
from mapping import gettrack
from reactool import rebond
from toolbox import erase_blank

#=======================================================================
# PURPOSE: fragments a species into 2 parts. The bond matrix provided
# as input contains the 2 fragments that must be converted as 2 new
# species provided as output (not standardized). It is thus assumed 
# that the bond breaking is already done in the calling program, by 
# setting the corresponding element of the bond matrix to zero
#=======================================================================
def fragm(bond, group):
    # 声明并初始化变量（保持原Fortran左闭右闭范围逻辑）
    ngr = sum(1 for g in group if g.strip() != '')  # 对应COUNT(group/=' ')
    len_group1 = len(group[0]) if group else 0
    
    # 初始化数组（维度完全匹配原Fortran）
    tgrp1 = [''] * len(group)
    tbnd1 = [[0 for _ in range(len(bond[0]))] for _ in range(len(bond))]
    trace1 = [0 for _ in range(len(bond))]
    
    tgrp2 = [''] * len(group)
    tbnd2 = [[0 for _ in range(len(bond[0]))] for _ in range(len(bond))]
    trace2 = [0 for _ in range(len(bond))]
    
    tgroup = [''] * len(group)
    tbond = [[0 for _ in range(len(bond[0]))] for _ in range(len(bond))]
    
    # 初始化track和trlen（对应原Fortran mxcp×SIZE(bond,1)维度）
    track = [[0 for _ in range(len(bond))] for _ in range(mxcp)]
    trlen = [0 for _ in range(mxcp)]
    ntr = 0
    
    # copy bond and groups（循环范围：1到len(group)左闭右闭）
    for i in range(1, len(group) + 1):
        tgroup[i - 1] = group[i - 1]
    # 循环范围：1到len(bond)左闭右闭，内层1到len(bond[0])左闭右闭
    for i in range(1, len(bond) + 1):
        for j in range(1, len(bond[0]) + 1):
            tbond[i - 1][j - 1] = bond[i - 1][j - 1]
    
    erase_blank(tbond, tgroup)
    
    # initialize tracer arrays
    chem1 = ' ' * (len(group) * len_group1)  # 适配rebond输出长度
    chem2 = ' ' * (len(group) * len_group1)
    # 循环范围：1到len(trace1)左闭右闭
    for i in range(1, len(trace1) + 1):
        trace1[i - 1] = 0
        trace2[i - 1] = 0
    for i in range(1, len(tgrp1) + 1):
        tgrp1[i - 1] = ' '
        tgrp2[i - 1] = ' '
    for i in range(1, len(tbnd1) + 1):
        for j in range(1, len(tbnd1[0]) + 1):
            tbnd1[i - 1][j - 1] = 0
            tbnd2[i - 1][j - 1] = 0
    
    # connect TRACE1 (if not connected, skip)
    trace1[0] = 1  # 原Fortran trace1(1)=1，列表索引转0-based
    gettrack(tbond, 1, ngr, ntr, track, trlen)  # 传入1（原Fortran 1-based起始节点）
    
    # 循环范围：1到ntr左闭右闭
    for i in range(1, ntr + 1):
        # 循环范围：1到trlen(i)左闭右闭
        for j in range(1, trlen[i - 1] + 1):
            node = track[i - 1][j - 1]  # track存储的是1-based节点
            if node != 0:
                trace1[node - 1] = 1  # 列表索引转0-based
    
    # construct TRACE2 from non-blank remainder（循环范围：1到ngr左闭右闭）
    for i in range(1, ngr + 1):
        if trace1[i - 1] == 0 and tgroup[i - 1].strip() != '':
            trace2[i - 1] = 1
    
    # write new group1, group2, bond1 and bond2 elements（循环范围：1到ngr左闭右闭）
    for i in range(1, ngr + 1):
        if trace1[i - 1] == 1:
            tgrp1[i - 1] = tgroup[i - 1]
            # 循环范围：1到ngr左闭右闭
            for j in range(1, ngr + 1):
                tbnd1[i - 1][j - 1] = tbond[i - 1][j - 1]
        
        if trace2[i - 1] == 1:
            tgrp2[i - 1] = tgroup[i - 1]
            # 循环范围：1到ngr左闭右闭
            for j in range(1, ngr + 1):
                tbnd2[i - 1][j - 1] = tbond[i - 1][j - 1]
    
    nring = 0
    rebond(tbnd1, tgrp1, chem1, nring)
    rebond(tbnd2, tgrp2, chem2, nring)
    
    return chem1, chem2