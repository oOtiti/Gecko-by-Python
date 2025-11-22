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
from keyparameter import mxcp
from mapping import alkenetrack
from toolbox import stoperr

def cdcase2(chem, bond, group, rxnflg, ncdtrack, cdtrack, xcdconjug, xcdsub, xcdeth, xcdcarbo, xcdcase):
    # 初始化变量
    ncdtrack[0] = 0  # 用列表包裹实现输出传递
    # 初始化二维数组cdtrack为0
    for i in range(len(cdtrack)):
        for j in range(len(cdtrack[i])):
            cdtrack[i][j] = 0
    # 初始化一维数组xcdconjug为0
    for i in range(len(xcdconjug)):
        xcdconjug[i] = 0
    # 初始化一维数组xcdsub为0
    for i in range(len(xcdsub)):
        xcdsub[i] = 0
    # 初始化二维数组xcdeth为0
    for i in range(len(xcdeth)):
        for j in range(len(xcdeth[i])):
            xcdeth[i][j] = 0
    # 初始化二维数组xcdcarbo为0
    for i in range(len(xcdcarbo)):
        for j in range(len(xcdcarbo[i])):
            xcdcarbo[i][j] = 0
    # 初始化一维数组xcdcase为0
    for i in range(len(xcdcase)):
        xcdcase[i] = 0
    
    mxcd = len(cdtrack[0])  # SIZE(cdtrack,2)
    # 初始化二维数组cdbond为0
    cdbond = [[0 for _ in range(len(group))] for _ in range(len(group))]
    
    # get the number of groups
    ngr = sum(1 for g in group if g.strip() != '')
    
    # ------------------------------
    # Create and check the Cd tracks
    # ------------------------------
    # 声明cdtracklen数组，长度为SIZE(cdtrack,1)
    cdtracklen = [0 for _ in range(len(cdtrack))]
    alkenetrack(chem, bond, group, ngr, ncdtrack, cdtracklen, cdtrack)
    for i in range(1, ncdtrack[0] + 1):
        if cdtracklen[i-1] == 4:
            xcdconjug[i-1] = 1
    
    # check for >C=CR-OH (enol: should have already tautomerised to ketone 
    # in subroutine alkcheck) and >C=CR-ONO2 (not available yet)
    if rxnflg != 0:
        for i in range(1, ncdtrack[0] + 1):
            for j in range(1, mxcd + 1):
                if cdtrack[i-1][j-1] != 0:
                    grp_idx = cdtrack[i-1][j-1] - 1  # 转换为0-based索引
                    if '(ONO2)' in group[grp_idx]:
                        mesg = ">C=C(ONO2)- not allowed"
                        stoperr('cdcase2 ', mesg, chem)
                    if '(OH)' in group[grp_idx]:
                        mesg = ">C=C(OH)- not allowed"
                        stoperr('cdcase2 ', mesg, chem)
    
    # ------------------------------
    # Fill the various matrix
    # ------------------------------
    
    # count # of nodes bonded to each Cd, except the C=C node (for which bond(i,j)=2)
    # and store -O- and CO nodes bonded to Cds
    for i in range(1, ncdtrack[0] + 1):
        for j in range(1, mxcd + 1):
            if cdtrack[i-1][j-1] != 0:
                icd = cdtrack[i-1][j-1]  # current Cd node is icd (1-based)
                for k in range(1, ngr + 1):  # check Cd neighbours (1-based)
                    if bond[icd-1][k-1] == 1:  # bond是0-based数组
                        xcdsub[icd-1] += 1  # count sub on Cd grp. (xcdsub是0-based)
                        grp_k = group[k-1]  # group是0-based
                        if grp_k == 'CHO' or grp_k[:2] == 'CO':
                            if xcdcarbo[icd-1][0] == 0:
                                xcdcarbo[icd-1][0] = k  # store CO node (1-based)
                            elif xcdcarbo[icd-1][1] == 0:
                                xcdcarbo[icd-1][1] = k  # store CO node (1-based)
                            else:
                                mesg = "more than 2 carbonyls bonded to a Cd"
                                stoperr('cdcase2 ', mesg, chem)
                    
                    if bond[icd-1][k-1] == 3:
                        xcdsub[icd-1] += 1  # count sub on Cd grp. (xcdsub是0-based)
                        if xcdeth[icd-1][0] == 0:
                            xcdeth[icd-1][0] = k  # store -O- node (1-based)
                        elif xcdeth[icd-1][1] == 0:
                            xcdeth[icd-1][1] = k  # store -O- node (1-based)
                        else:
                            mesg = "more than 2 ethers bonded to a Cd"
                            stoperr('cdcase2 ', mesg, chem)
    
    # -------------------------------
    # set "cd case" for each Cd track 
    # -------------------------------
    trloop = True
    for i in range(1, ncdtrack[0] + 1):
        # check ketene >C=C=O (case 6)
        case6_flag = False
        for j in range(1, mxcd + 1):
            icd = cdtrack[i-1][j-1]
            if icd == 0:
                break
            grp_idx = icd - 1  # 转换为0-based
            if group[grp_idx] == 'CdO ':
                xcdcase[i-1] = 6
                case6_flag = True
                break
        if case6_flag:
            continue
        
        # check vinylic structure C=C-OR (case 7)
        case7_flag = False
        for j in range(1, mxcd + 1, 2):
            icd1st = cdtrack[i-1][j-1]
            if icd1st == 0:
                break
            # j+1不超过mxcd（因为mxcd是偶数？按原逻辑）
            if j + 1 <= mxcd:
                icd2nd = cdtrack[i-1][j]  # j是1-based，j-1是当前索引，j是下一个
            else:
                icd2nd = 0
            
            if xcdeth[icd1st-1][0] != 0:
                xcdcase[i-1] = 7
                case7_flag = True
            if icd2nd != 0 and xcdeth[icd2nd-1][0] != 0:
                xcdcase[i-1] = 7
                case7_flag = True
            
            if (icd1st != 0 and xcdeth[icd1st-1][0] != 0) and (icd2nd != 0 and xcdeth[icd2nd-1][0] != 0):
                mesg = "ether found both sides of a Cd=Cd"
                stoperr('cdcase2 ', mesg, chem)
            
            if case7_flag:
                break
        if case7_flag:
            continue
        
        # check conjugated C=C-C=C tracks (case 1, 3 or 5)
        if xcdconjug[i-1] != 0:
            nco = 0
            # check terminal Cd only for -CO- (j=1和j=4，1-based)
            for j in [1, 4]:
                if j > mxcd:
                    break
                icd = cdtrack[i-1][j-1]
                if icd == 0:
                    raise RuntimeError('--error-- in cdcase, no Cd unexpected')
                if xcdcarbo[icd-1][0] != 0:
                    nco += 1
            if nco > 1:
                xcdcase[i-1] = 3
            elif nco == 1:
                xcdcase[i-1] = 5
            else:
                xcdcase[i-1] = 1
            continue
        
        # check simple C=C tracks (case 1 or 2)
        else:
            nco = 0
            for j in range(1, 3):  # j=1,2 (1-based)
                if j > mxcd:
                    break
                icd = cdtrack[i-1][j-1]
                if icd == 0:
                    raise RuntimeError('--error-- in cdcase, no Cd unexpected')
                if xcdcarbo[icd-1][0] != 0:
                    nco += 1
            if nco != 0:
                xcdcase[i-1] = 2
            else:
                xcdcase[i-1] = 1
            continue

#=======================================================================
# PURPOSE : switch enol structure to the corresponding keto structure.
# NOTE: the routine is similar to alkcheck below, but limited to enol
# only.
#=======================================================================
from keyparameter import mxnode, mxlgr, mxring
from rjtool import rjgrm
from stdgrbond import grbond
from reactool import swap, rebond
from normchem import stdchm

def switchenol(pchem, loswitch):
    # 初始化输出变量
    loswitch[0] = False  # 用列表包裹实现输出传递
    
    if '.' in pchem:
        return
    
    # 初始化变量
    bond = [[0 for _ in range(mxnode)] for _ in range(mxnode)]
    dbflg = 0
    nring = 0
    group = [''] * mxnode
    pold = ''
    pnew = ''
    tgroup = [''] * mxnode
    tempkc = [''] * len(pchem)
    tempkc = ''.join(tempkc)  # 初始化空字符串
    rjg = [[0 for _ in range(2)] for _ in range(mxring)]  # ring-join group pairs
    
    grbond(pchem, group, bond, dbflg, nring)
    rjgrm(nring, group, rjg)
    
    # grloop循环
    for i in range(1, mxnode + 1):  # 1-based到mxnode
        grp_idx = i - 1  # 转换为0-based
        if 'Cd' in group[grp_idx]:
            if '(OH)' in group[grp_idx]:
                pold = '(OH)'
                pnew = 'O'
                swap(group[grp_idx], pold, tgroup[grp_idx], pnew)
                group[grp_idx] = tgroup[grp_idx]
                
                pold = 'Cd'
                pnew = 'C'
                swap(group[grp_idx], pold, tgroup[grp_idx], pnew)
                group[grp_idx] = tgroup[grp_idx]
                
                for j in range(1, mxnode + 1):  # 1-based到mxnode
                    if bond[i-1][j-1] == 2:  # bond是0-based
                        bond[i-1][j-1] = 1
                        bond[j-1][i-1] = 1
                        
                        grp_j = group[j-1]
                        if grp_j[:4] == 'CdH2':
                            pold = 'CdH2'
                            pnew = 'CH3'
                        elif grp_j[:3] == 'CdH':
                            pold = 'CdH'
                            pnew = 'CH2'
                        elif grp_j[:2] == 'Cd':
                            pold = 'Cd'
                            pnew = 'CH'
                        else:
                            # 无匹配基团，不执行swap
                            continue
                        
                        swap(group[j-1], pold, tgroup[j-1], pnew)
                        group[j-1] = tgroup[j-1]
                        rebond(bond, group, tempkc, nring)
                        # 更新pchem（输入输出）
                        pchem[0] = tempkc  # 用列表包裹pchem实现输出传递
                        stdchm(pchem[0])
                
                loswitch[0] = True

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
def alkcheck(pchem, coprod, acom):
    # 初始化输出变量
    coprod[0] = ' '  # 用列表包裹实现输出传递
    acom[0] = ' '    # 用列表包裹实现输出传递
    
    if '.' in pchem[0]:
        return
    
    # 初始化变量
    bond = [[0 for _ in range(mxnode)] for _ in range(mxnode)]
    dbflg = 0
    nring = 0
    group = [''] * mxnode
    pold = ''
    pnew = ''
    tgroup = [''] * mxnode
    tempkc = [''] * len(pchem[0])
    tempkc = ''.join(tempkc)  # 初始化空字符串
    rjg = [[0 for _ in range(2)] for _ in range(mxring)]  # ring-join group pairs
    
    grbond(pchem[0], group, bond, dbflg, nring)
    rjgrm(nring, group, rjg)
    
    # grloop循环
    for i in range(1, mxnode + 1):  # 1-based到mxnode
        grp_idx = i - 1  # 转换为0-based
        if 'Cd' in group[grp_idx]:
            if '(OH)' in group[grp_idx]:
                pold = '(OH)'
                pnew = 'O'
                swap(group[grp_idx], pold, tgroup[grp_idx], pnew)
                group[grp_idx] = tgroup[grp_idx]
                
                pold = 'Cd'
                pnew = 'C'
                swap(group[grp_idx], pold, tgroup[grp_idx], pnew)
                group[grp_idx] = tgroup[grp_idx]
                
                for j in range(1, mxnode + 1):  # 1-based到mxnode
                    if bond[i-1][j-1] == 2:  # bond是0-based
                        bond[i-1][j-1] = 1
                        bond[j-1][i-1] = 1
                        
                        grp_j = group[j-1]
                        if grp_j[:4] == 'CdH2':
                            pold = 'CdH2'
                            pnew = 'CH3'
                        elif grp_j[:3] == 'CdH':
                            pold = 'CdH'
                            pnew = 'CH2'
                        elif grp_j[:2] == 'Cd':
                            pold = 'Cd'
                            pnew = 'CH'
                        else:
                            # 无匹配基团，不执行swap
                            continue
                        
                        swap(group[j-1], pold, tgroup[j-1], pnew)
                        group[j-1] = tgroup[j-1]
                        rebond(bond, group, tempkc, nring)
                        # 更新pchem（输入输出）
                        pchem[0] = tempkc
                        stdchm(pchem[0])
                
                acom[0] = 'KETOENOL '
                break  # EXIT grloop
            
            # ELSE IF (INDEX(group(i),'(OOH)')/=0) THEN
            #     pold='(OOH)'  ;  pnew='(O.)'
            #     CALL swap(group(i),pold,tgroup(i),pnew)
            #     group(i)=tgroup(i)
            #     CALL rebond(bond,group,tempkc,nring)
            #     pchem=tempkc
            #     CALL stdchm (pchem)
            #     coprod='HO'
            #     acom='xxxxxxxx'
            #     EXIT grloop
            
            elif '(ONO2)' in group[grp_idx]:
                pold = '(ONO2)'
                pnew = '(O.)'
                swap(group[grp_idx], pold, tgroup[grp_idx], pnew)
                group[grp_idx] = tgroup[grp_idx]
                
                rebond(bond, group, tempkc, nring)
                pchem[0] = tempkc
                stdchm(pchem[0])
                
                coprod[0] = 'NO2'
                acom[0] = 'KETOENENIT'
                break  # EXIT grloop