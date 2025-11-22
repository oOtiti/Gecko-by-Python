def fragm(bond, group, chem1, chem2):
    from keyparameter import mxcp, mxnode, mxlfo, mxlgr
    from mapping import gettrack
    from reactool import rebond

    ngr = 0
    for i in range(1, len(group) + 1):
        if group[i-1] != ' ':
            ngr += 1
    
    # copy bond and groups
    tgroup = [''] * len(group)
    for i in range(1, len(group) + 1):
        tgroup[i-1] = group[i-1]
    
    tbond = [[0 for _ in range(len(bond[0]))] for _ in range(len(bond))]
    for i in range(1, len(bond) + 1):
        for j in range(1, len(bond[0]) + 1):
            tbond[i-1][j-1] = bond[i-1][j-1]
    
    # erase blank nodes
    erase_blank(tbond, tgroup)
    
    # initialize tracer arrays
    chem1 = ''
    chem2 = ''
    trace1 = [0] * len(bond)
    trace2 = [0] * len(bond)
    tgrp1 = [''] * len(group)
    tgrp2 = [''] * len(group)
    tbnd1 = [[0 for _ in range(len(bond[0]))] for _ in range(len(bond))]
    tbnd2 = [[0 for _ in range(len(bond[0]))] for _ in range(len(bond))]
    
    # connect TRACE1
    trace1[0] = 1
    ntr = 0
    track = [[0 for _ in range(len(bond))] for _ in range(mxcp)]
    trlen = [0] * mxcp
    gettrack(tbond, 1, ngr, ntr, track, trlen)
    
    for i in range(1, ntr + 1):
        for j in range(1, trlen[i-1] + 1):
            trace1[track[i-1][j-1] - 1] = 1
    
    # construct TRACE2 from non-blank remainder
    for i in range(1, ngr + 1):
        if trace1[i-1] == 0 and tgroup[i-1] != ' ':
            trace2[i-1] = 1
    
    # write new group1, group2, bond1 and bond2 elements
    for i in range(1, ngr + 1):
        if trace1[i-1] == 1:
            tgrp1[i-1] = tgroup[i-1]
            for j in range(1, ngr + 1):
                tbnd1[i-1][j-1] = tbond[i-1][j-1]
        
        if trace2[i-1] == 1:
            tgrp2[i-1] = tgroup[i-1]
            for j in range(1, ngr + 1):
                tbnd2[i-1][j-1] = tbond[i-1][j-1]
    
    nring1 = 0
    nring2 = 0
    chem1, nring1 = rebond(tbnd1, tgrp1, chem1, nring1)
    chem2, nring2 = rebond(tbnd2, tgrp2, chem2, nring2)
    
    return chem1, chem2

def erase_blank(tbond, tgroup):
    ngr = len(tgroup)
    keep_indices = []
    
    # Find indices of non-blank groups (1-based)
    for i in range(1, ngr + 1):
        if tgroup[i-1] != ' ':
            keep_indices.append(i)
    
    # Create new arrays with same dimensions
    new_tbond = [[0 for _ in range(len(tbond[0]))] for _ in range(len(tbond))]
    new_tgroup = [''] * len(tgroup)
    
    # Copy non-blank groups and bonds to new arrays
    for i, idx_i in enumerate(keep_indices, 1):
        if i <= len(new_tgroup):
            new_tgroup[i-1] = tgroup[idx_i-1]
        for j, idx_j in enumerate(keep_indices, 1):
            if i <= len(new_tbond) and j <= len(new_tbond[0]):
                new_tbond[i-1][j-1] = tbond[idx_i-1][idx_j-1]
    
    # Fill remaining positions with blanks/zeros
    for i in range(len(keep_indices) + 1, len(new_tgroup) + 1):
        new_tgroup[i-1] = ' '
    for i in range(len(keep_indices) + 1, len(new_tbond) + 1):
        for j in range(1, len(new_tbond[0]) + 1):
            new_tbond[i-1][j-1] = 0
    for j in range(len(keep_indices) + 1, len(new_tbond[0]) + 1):
        for i in range(1, len(new_tbond) + 1):
            new_tbond[i-1][j-1] = 0
    
    # Update the input arrays
    for i in range(1, len(tgroup) + 1):
        tgroup[i-1] = new_tgroup[i-1]
    for i in range(1, len(tbond) + 1):
        for j in range(1, len(tbond[0]) + 1):
            tbond[i-1][j-1] = new_tbond[i-1][j-1]