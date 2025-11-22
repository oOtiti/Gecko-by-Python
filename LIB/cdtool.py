# 前置库已由用户提供，此处仅按需求导入对应模块内容
from keyparameter import mxcp, mxnode, mxlgr, mxring
from mapping import alkenetrack
from toolbox import stoperr
from rjtool import rjgrm
from stdgrbond import grbond
from reactool import swap
from normchem import stdchm
from reactool import rebond
#=======================================================================
# PURPOSE : check whether the C=C bonds are conjugated or not, including
# C=O bond (i.e. structure of type -C=C-C=O). The subroutine returns the
# "case" to which the species belongs. Six cases are considered :            
# CASE 1: regular Cd molecule : only >C=C< and >C=C-C=C< bonds in     
#         the molecule (i.e, without conjugated C=C-C=O)       
# CASE 2: is for structures containing the >C=C-C=O structure  
#         but no C=C-C=C (or C=C=O)      
# CASE 3: is for the -CO-C=C-C=C-C=O structure only (i.e. containing  
#         carbonyl at both sides of the conjugated C=C-C=C)     
# CASE 4: Two double bonds non conjugated (i.e. not C=C-C=C) but
#         containing at least one C=C-C=O                       
# CASE 5: is for the -CO-C=C-C=C< structure (i.e. containing    
#         carbonyl at only one side of the conjugated C=C-C=C)  
# CASE 6: -C=C=O e.g. ketene (only case for this group)         
# CASE 7: -C=C-O- => vinyl ether chemistry      
#=======================================================================
def cdcase2(chem, bond, group, rxnflg):
    # 声明并初始化输出变量（完全对应原Fortran维度，1-based逻辑保留）
    mxcd = len(bond[0])  # 对应SIZE(cdtrack,2)，保持原维度
    ncdtrack = 0
    # cdtrack保持原Fortran范围（左闭右闭），数组大小不变
    cdtrack = [[0 for _ in range(mxcd)] for _ in range(len(bond))]
    xcdconjug = [0 for _ in range(len(bond))]
    xcdsub = [0 for _ in range(len(group))]
    xcdeth = [[0 for _ in range(2)] for _ in range(len(group))]
    xcdcarbo = [[0 for _ in range(2)] for _ in range(len(group))]
    xcdcase = [0 for _ in range(len(bond))]
    
    cdtracklen = [0 for _ in range(len(bond))]
    cdbond = [[0 for _ in range(len(group))] for _ in range(len(group))]
    
    # get the number of groups
    ngr = sum(1 for g in group if g.strip() != '')
    
    # ------------------------------
    # Create and check the Cd tracks
    # ------------------------------
    alkenetrack(chem, bond, group, ngr, ncdtrack, cdtracklen, cdtrack)
    # 循环范围完全保留原Fortran左闭右闭：i从1到ncdtrack（含边界）
    for i in range(1, ncdtrack + 1):
        if cdtracklen[i - 1] == 4:  # 列表存储为0-based，索引转换不改变循环范围
            xcdconjug[i - 1] = 1
    
    # check for >C=CR-OH (enol: should have already tautomerised to ketone 
    # in subroutine alkcheck) and >C=CR-ONO2 (not available yet)
    if rxnflg != 0:
        # 循环范围：i从1到ncdtrack（左闭右闭），j从1到mxcd（左闭右闭）
        for i in range(1, ncdtrack + 1):
            for j in range(1, mxcd + 1):
                if cdtrack[i - 1][j - 1] != 0:
                    idx = cdtrack[i - 1][j - 1]  # 保留原Fortran 1-based索引
                    if '(ONO2)' in group[idx - 1]:  # 列表访问转0-based
                        mesg = ">C=C(ONO2)- not allowed"
                        stoperr('cdcase2 ', mesg, chem)
                    if '(OH)' in group[idx - 1]:
                        mesg = ">C=C(OH)- not allowed"
                        stoperr('cdcase2 ', mesg, chem)
    
    # ------------------------------
    # Fill the various matrix
    # ------------------------------
    # count # of nodes bonded to each Cd, except the C=C node (for which bond(i,j)=2)
    # and store -O- and CO nodes bonded to Cds
    # 循环范围：i从1到ncdtrack（左闭右闭），j从1到mxcd（左闭右闭）
    for i in range(1, ncdtrack + 1):
        for j in range(1, mxcd + 1):
            if cdtrack[i - 1][j - 1] != 0:
                icd = cdtrack[i - 1][j - 1]  # current Cd node（保留1-based）
                # 循环范围：k从1到ngr（左闭右闭）
                for k in range(1, ngr + 1):
                    if bond[icd - 1][k - 1] == 1:  # 列表访问转0-based
                        xcdsub[icd - 1] += 1
                        if group[k - 1] == 'CHO' or group[k - 1][:2] == 'CO':
                            if xcdcarbo[icd - 1][0] == 0:
                                xcdcarbo[icd - 1][0] = k  # 保留1-based存储
                            elif xcdcarbo[icd - 1][1] == 0:
                                xcdcarbo[icd - 1][1] = k
                            else:
                                mesg = "more than 2 carbonyls bonded to a Cd"
                                stoperr('cdcase2 ', mesg, chem)
                    
                    if bond[icd - 1][k - 1] == 3:
                        xcdsub[icd - 1] += 1
                        if xcdeth[icd - 1][0] == 0:
                            xcdeth[icd - 1][0] = k  # 保留1-based存储
                        elif xcdeth[icd - 1][1] == 0:
                            xcdeth[icd - 1][1] = k
                        else:
                            mesg = "more than 2 ethers bonded to a Cd"
                            stoperr('cdcase2 ', mesg, chem)
    
    # -------------------------------
    # set "cd case" for each Cd track 
    # -------------------------------
    # 循环范围：i从1到ncdtrack（左闭右闭）
    for i in range(1, ncdtrack + 1):
        # check ketene >C=C=O (case 6)
        ketene_flag = False
        # 循环范围：j从1到mxcd（左闭右闭）
        for j in range(1, mxcd + 1):
            if cdtrack[i - 1][j - 1] == 0:
                break
            icd = cdtrack[i - 1][j - 1]
            if group[icd - 1] == 'CdO ':
                xcdcase[i - 1] = 6
                ketene_flag = True
                break
        if ketene_flag:
            continue
        
        # check vinylic structure C=C-OR (case 7)
        vinyl_flag = False
        # 循环范围：j从1到mxcd，步长2（左闭右闭，对应原Fortran j=1,mxcd,2）
        for j in range(1, mxcd + 1, 2):
            icd1st = cdtrack[i - 1][j - 1]
            if icd1st == 0:
                break
            # j+1不超过mxcd（保持原范围约束）
            if j + 1 > mxcd:
                icd2nd = 0
            else:
                icd2nd = cdtrack[i - 1][j]
            
            if xcdeth[icd1st - 1][0] != 0:
                xcdcase[i - 1] = 7
                vinyl_flag = True
            if icd2nd != 0 and xcdeth[icd2nd - 1][0] != 0:
                xcdcase[i - 1] = 7
                vinyl_flag = True
            
            if (xcdeth[icd1st - 1][0] != 0) and (icd2nd != 0 and xcdeth[icd2nd - 1][0] != 0):
                mesg = "ether found both sides of a Cd=Cd"
                stoperr('cdcase2 ', mesg, chem)
            
            if vinyl_flag:
                break
        if vinyl_flag:
            continue
        
        # check conjugated C=C-C=C tracks (case 1, 3 or 5)
        if xcdconjug[i - 1] != 0:
            nco = 0
            # 循环范围：j从1到4，步长3（左闭右闭，对应原Fortran j=1,4,3）
            for j in range(1, 5, 3):
                if j > mxcd:
                    break
                icd = cdtrack[i - 1][j - 1]
                if icd == 0:
                    raise Exception('--error-- in cdcase, no Cd unexpected')
                if xcdcarbo[icd - 1][0] != 0:
                    nco += 1
            if nco > 1:
                xcdcase[i - 1] = 3
                continue
            elif nco == 1:
                xcdcase[i - 1] = 5
                continue
            else:
                xcdcase[i - 1] = 1
                continue
        
        # check simple C=C tracks (case 1 or 2)
        else:
            nco = 0
            # 循环范围：j从1到2（左闭右闭）
            for j in range(1, 3):
                if j > mxcd:
                    break
                icd = cdtrack[i - 1][j - 1]
                if icd == 0:
                    raise Exception('--error-- in cdcase, no Cd unexpected')
                if xcdcarbo[icd - 1][0] != 0:
                    nco += 1
            if nco != 0:
                xcdcase[i - 1] = 2
                continue
            else:
                xcdcase[i - 1] = 1
                continue
    
    return ncdtrack, cdtrack, xcdconjug, xcdsub, xcdeth, xcdcarbo, xcdcase

#=======================================================================
# PURPOSE : switch enol structure to the corresponding keto structure.
# NOTE: the routine is similar to alkcheck below, but limited to enol
# only.
#=======================================================================
def switchenol(pchem):
    # 声明并初始化变量（保留原范围逻辑）
    loswitch = False
    # bond维度：mxnode×mxnode（左闭右闭，对应原Fortran范围）
    bond = [[0 for _ in range(mxnode)] for _ in range(mxnode)]
    dbflg = 0
    nring = 0
    group = [''] * mxnode
    tgroup = [''] * mxnode
    pold = ''
    pnew = ''
    tempkc = ' ' * len(pchem)
    # rjg维度：mxring×2（左闭右闭）
    rjg = [[0 for _ in range(2)] for _ in range(mxring)]
    
    if '.' in pchem:
        return pchem, loswitch
    
    grbond(pchem, group, bond, dbflg, nring)
    rjgrm(nring, group, rjg)
    
    # 循环范围：i从1到mxnode（左闭右闭）
    for i in range(1, mxnode + 1):
        if 'Cd' in group[i - 1]:
            if '(OH)' in group[i - 1]:
                pold = '(OH)'
                pnew = 'O'
                swap(group[i - 1], pold, tgroup[i - 1], pnew)
                group[i - 1] = tgroup[i - 1]
                
                pold = 'Cd'
                pnew = 'C'
                swap(group[i - 1], pold, tgroup[i - 1], pnew)
                group[i - 1] = tgroup[i - 1]
                
                # 循环范围：j从1到mxnode（左闭右闭）
                for j in range(1, mxnode + 1):
                    if bond[i - 1][j - 1] == 2:
                        bond[i - 1][j - 1] = 1
                        bond[j - 1][i - 1] = 1
                        
                        if group[j - 1][:4] == 'CdH2':
                            pold = 'CdH2'
                            pnew = 'CH3'
                        elif group[j - 1][:3] == 'CdH':
                            pold = 'CdH'
                            pnew = 'CH2'
                        elif group[j - 1][:2] == 'Cd':
                            pold = 'Cd'
                            pnew = 'CH'
                        
                        swap(group[j - 1], pold, tgroup[j - 1], pnew)
                        group[j - 1] = tgroup[j - 1]
                        rebond(bond, group, tempkc, nring)
                        pchem = tempkc
                        stdchm(pchem)
                
                loswitch = True
    
    return pchem, loswitch

#=======================================================================
#=======================================================================
# PURPOSE : This subroutine fragments the substitued alkenes such as
# >Cd=Cd(OH)- or >Cd=Cd(OOH)- or >Cd=Cd(ONO2)- which come from Norrish
# II (alkenes photolysis), or from fragmentation after oxidation of
# conjugated alkenes.
# We consider that these alkenes are energy-rich and decompose to :
#  >Cd=Cd(OH)-    -> >CH-CO- 
#  >Cd=Cd(OOH)-   -> >C(.)-CO- + OH.
#  >Cd=Cd(ONO2)-  -> >C(.)-CO- + NO2
# Each Cd can't support more than one group like OH, OOH or ONO2
# If there is more than one group on the double bond, treat first the
# -OH next the -OOH and at last the -ONO2
# Note: the routine is similar the switchenol routine above. This 
# alkcheck is called after adding the enol to the dictionary and 
# often act as a "patch" for treating enol structure. Must progressively
# be avoid with the progressive development of gecko.
#=======================================================================
def alkcheck(pchem):
    # 声明并初始化变量（保留原范围逻辑）
    coprod = ' '
    acom = ' '
    # bond维度：mxnode×mxnode（左闭右闭）
    bond = [[0 for _ in range(mxnode)] for _ in range(mxnode)]
    dbflg = 0
    nring = 0
    group = [''] * mxnode
    tgroup = [''] * mxnode
    pold = ''
    pnew = ''
    tempkc = ' ' * len(pchem)
    # rjg维度：mxring×2（左闭右闭）
    rjg = [[0 for _ in range(2)] for _ in range(mxring)]
    
    if '.' in pchem:
        return pchem, coprod, acom
    
    grbond(pchem, group, bond, dbflg, nring)
    rjgrm(nring, group, rjg)
    
    # 循环范围：i从1到mxnode（左闭右闭）
    for i in range(1, mxnode + 1):
        if 'Cd' in group[i - 1]:
            if '(OH)' in group[i - 1]:
                pold = '(OH)'
                pnew = 'O'
                swap(group[i - 1], pold, tgroup[i - 1], pnew)
                group[i - 1] = tgroup[i - 1]
                
                pold = 'Cd'
                pnew = 'C'
                swap(group[i - 1], pold, tgroup[i - 1], pnew)
                group[i - 1] = tgroup[i - 1]
                
                # 循环范围：j从1到mxnode（左闭右闭）
                for j in range(1, mxnode + 1):
                    if bond[i - 1][j - 1] == 2:
                        bond[i - 1][j - 1] = 1
                        bond[j - 1][i - 1] = 1
                        
                        if group[j - 1][:4] == 'CdH2':
                            pold = 'CdH2'
                            pnew = 'CH3'
                        elif group[j - 1][:3] == 'CdH':
                            pold = 'CdH'
                            pnew = 'CH2'
                        elif group[j - 1][:2] == 'Cd':
                            pold = 'Cd'
                            pnew = 'CH'
                        
                        swap(group[j - 1], pold, tgroup[j - 1], pnew)
                        group[j - 1] = tgroup[j - 1]
                        rebond(bond, group, tempkc, nring)
                        pchem = tempkc
                        stdchm(pchem)
                
                acom = 'KETOENOL '
                break
            
            # 注释部分保持原Fortran逻辑
            # elif '(OOH)' in group[i - 1]:
            #     pold = '(OOH)'
            #     pnew = '(O.)'
            #     swap(group[i - 1], pold, tgroup[i - 1], pnew)
            #     group[i - 1] = tgroup[i - 1]
            #     rebond(bond, group, tempkc, nring)
            #     pchem = tempkc
            #     stdchm(pchem)
            #     coprod = 'HO'
            #     acom = 'xxxxxxxx'
            #     break
            
            elif '(ONO2)' in group[i - 1]:
                pold = '(ONO2)'
                pnew = '(O.)'
                swap(group[i - 1], pold, tgroup[i - 1], pnew)
                group[i - 1] = tgroup[i - 1]
                rebond(bond, group, tempkc, nring)
                pchem = tempkc
                stdchm(pchem)
                coprod = 'NO2'
                acom = 'KETOENENIT '
                break
    
    return pchem, coprod, acom