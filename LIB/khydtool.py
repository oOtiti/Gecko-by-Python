import math
# 前置库已由用户提供，此处仅按需求导入对应模块内容
from database import nkhydb, khydb_chem, khydb_dat
from searching import srh5
from keyparameter import mxnode, mxlgr, mxring, mxcp
from rjtool import rjgrm
from stdgrbond import grbond
from reactool import rebond
from mapping import abcde_map, chemmap
from reactool import swap
from normchem import stdchm
from toolbox import stoperr

# ======================================================================
# PURPOSE: get all hydrates that can be made from a given formula and 
# related hydration constants. 
#
# For details about the Gromhe SAR and Khydration, see 
# T. Raventos-Duran et al., Atmospheric Chemistry and Physics, 7643-7654, 2010.
# C. Mouchel-Vallon et al., Geoscientific Model Development, 1339-1362, 2017 
# https://doi.org/10.5194/gmd-10-1339-2017, 
# 
# OUTPUT SUMMARY:
# - nwa: number of water molecule that can be added to chem. (e.g. if 
#     the molecule bear 3 ketones, the nwa is return as 3).
# - nhyd(i): number of distinct hydrate that can be made for i molecule
#     of water added.
# - chemhyd(j,i): formula of the jth hydrate bearing ith added water 
#     molecule
# - yhyd(j,i): hydration constant (with respect to the fully non-hydrated 
#     molec) of the jth hydrate bearing ith added water molecules (log scale)
# - khydstar: hydration constant taking all possible hydrate into
#     account (apparent constant K_star). NOT A LOG SCALE !
# ======================================================================
def get_hydrate(chem, nwa, nhyd, chemhyd, yhyd, khydstar):
    from database import nkhydb, khydb_chem, khydb_dat
    from searching import srh5
    # 声明并初始化局部变量（结合Fortran声明与初始化）
    i = 0
    j = 0
    t1chemhyd = [""] * len(yhyd[0])  # 维度匹配yhyd第2维
    t1yhyd = [0.0 for _ in range(len(yhyd[0]))]
    t1nhyd = 0

    chemdat = [""] * len(nhyd)  # 维度匹配nhyd
    tempkc = ""
    yield_ = 0.0  # 避免与Python关键字冲突
    numy = 0.0
    sumy = 0.0
    ydat = [0.0 for _ in range(len(nhyd))]
    ndat = 0
    jp = 0
    hlev = 0

    # table of functions - Index of functionalities
    # --------------------------
    #  1= -OH    ;  2= -NO2     ;  3= -ONO2   ; 4= -OOH    ;  5= -F     
    #  6= -Cl    ;  7=  -Br     ;  8= -I      ; 9= -CHO    ; 10= -CO-   
    # 11= -COOH  ; 12= -CO(OOH) ; 13= -PAN    ; 14= -O-    ; 15= R-COO-R
    # 16 = HCO-O-R; 17= -CO(F) ; 18= -CO(Cl) ;  19= -CO(Br) ; 20= -CO(I)

    # 初始化输出变量（严格匹配Fortran逻辑）
    nwa = 0
    khydstar = 0.0
    for idx in range(len(nhyd)):
        nhyd[idx] = 0
    for idx1 in range(len(chemhyd)):
        for idx2 in range(len(chemhyd[idx1])):
            chemhyd[idx1][idx2] = ' '
    for idx1 in range(len(yhyd)):
        for idx2 in range(len(yhyd[idx1])):
            yhyd[idx1][idx2] = 0.0
      
    # -----------
    # monohydrate
    # -----------
    t1nhyd = 0
    for idx in range(len(t1chemhyd)):
        t1chemhyd[idx] = ' '
    for idx in range(len(t1yhyd)):
        t1yhyd[idx] = 0.0

    yield_ = 0.0  # parent compound
    khydration(chem, yield_, ndat, chemdat, ydat)
    if ndat == 0:
        return  # 无水合物时直接返回（通过引用修改输出参数）
    t1nhyd = ndat
    for idx in range(1, ndat + 1):
        t1chemhyd[idx - 1] = chemdat[idx - 1]
        t1yhyd[idx - 1] = ydat[idx - 1]
      
    # collapse identical formula
    for i in range(1, t1nhyd):
        numy = 1.0
        sumy = t1yhyd[i - 1]
        for j in range(i + 1, t1nhyd + 1):
            if t1chemhyd[i - 1] == t1chemhyd[j - 1]:
                numy += 1.0
                sumy += t1yhyd[j - 1]
                t1chemhyd[j - 1] = ' '
                t1yhyd[j - 1] = 0.0
        if numy > 1.0:
            t1yhyd[i - 1] = sumy / float(numy)
      
    # write the table for monohydrate
    for i in range(1, t1nhyd + 1):
        if t1chemhyd[i - 1][:1] != ' ':
            nhyd[0] += 1
            chemhyd[0][nhyd[0] - 1] = t1chemhyd[i - 1]
            yhyd[0][nhyd[0] - 1] = t1yhyd[i - 1]

    # -----------
    # multi-hydrate
    # -----------
    multiloop = True
    hlev = 2
    while multiloop and hlev <= len(nhyd):
        if nhyd[hlev - 2] == 0:  # all hydrate were found (hlev-1对应索引hlev-2)
            nwa = hlev - 1
            multiloop = False
            break

        t1nhyd = 0
        for idx in range(len(t1chemhyd)):
            t1chemhyd[idx] = ' '
        for idx in range(len(t1yhyd)):
            t1yhyd[idx] = 0.0
        
        for jp in range(1, nhyd[hlev - 2] + 1):
            yield_ = yhyd[hlev - 2][jp - 1]  # parent compound (hlev-1对应索引hlev-2)
            tempkc = chemhyd[hlev - 2][jp - 1]
            khydration(tempkc, yield_, ndat, chemdat, ydat)
            for i in range(1, ndat + 1):
                t1chemhyd[t1nhyd + i - 1] = chemdat[i - 1]
                t1yhyd[t1nhyd + i - 1] = ydat[i - 1]
            t1nhyd += ndat

        # collapse identical formula
        for i in range(1, t1nhyd):
            numy = 1.0
            sumy = t1yhyd[i - 1]
            for j in range(i + 1, t1nhyd + 1):
                if t1chemhyd[i - 1] == t1chemhyd[j - 1]:
                    numy += 1.0
                    sumy += t1yhyd[j - 1]
                    t1chemhyd[j - 1] = ' '
                    t1yhyd[j - 1] = 0.0
            if numy > 1.0:
                t1yhyd[i - 1] = sumy / float(numy)
      
        # write the table for monohydrate
        for i in range(1, t1nhyd + 1):
            if t1chemhyd[i - 1][:1] != ' ':
                nhyd[hlev - 1] += 1
                chemhyd[hlev - 1][nhyd[hlev - 1] - 1] = t1chemhyd[i - 1]
                yhyd[hlev - 1][nhyd[hlev - 1] - 1] = t1yhyd[i - 1]

        hlev += 1

    # Adjust nwa if exit did not occur in the formula loop
    if hlev >= len(nhyd):
        nwa = len(nhyd)

    # compute Kstar
    for i in range(1, nwa + 1):
        for j in range(1, nhyd[i - 1] + 1):
            khydstar += 10 ** (yhyd[i - 1][j - 1])

    #-------------------------------------------------------------
    # Check if hydration constant is already known in the database
    #-------------------------------------------------------------
    i = srh5(chem, khydb_chem, nkhydb)
    if i > 0:
        khydstar = 10 ** khydb_dat[i - 1][0]

# ======================================================================
# PURPOSE: return the list of hydrate (1 water molecule added) that
# can be made from the formula given as input
#
# OUTPUT SUMMARY:
# - ndat : number of distinct molecule (hydrate) that can be made from
#     chem (1 water molecule added only !)
# - chemdat(i) : formula of the i'th hydrate
# - ydat(i) : hydration constant of the i'th hydrate (with respect
#             to the fully non hydrated molecule !)
#
# INTERNAL
# - nabcde(k): number of distinct pathways that end up at a position k
#     relative to top (e.g. nabcde(3) gives the number of distinct 
#     pathways finishing in a beta position relative to top.
# - tabcde(k,i,j): give the pathways (nodes j), for the track number i
#     to reach the position k (k=2 is beta position ...). For example, 
#     tabcde(4,1,j) give the first track to reach a gamma position (node
#     given by tabcde(4,1,4), using the track given by tabcde(4,1,*).
# - mapfun(a,b,c): provide the number of function of type 'c' at position
#     (node) 'a'. index 'b' if for node type with 1=aliphatic, 2=cd and
#     3=aromatic. For example, the molecule CH2(OH)CdH=CdHCHO should 
#     have non zero values at the positions : mapfun(1,1,1)=1 and 
#     mapfun(4,2,9)=1
# - funflg(a): get the number of functional group at node a. For the 
#     example above, non-zero values are found at position 1 and 4, 
#     where it is set to 1.
# - nodetype : table of character for type node:
#      'y' = carbonyl      'r' = aromatic        'o'= -O- node
#      'd' = Cd            'n' = others (i.e. normal)
#
# table of functions - Index of functionalities (in mapfun)
# -----------------------------------------------------------------
#  1= -OH    ;  2= -NO2     ;  3= -ONO2   ; 4= -OOH    ;  5= -F
#  6= -Cl    ;  7=  -Br     ;  8= -I      ; 9= -CHO    ; 10= -CO-
# 11= -COOH  ; 12= -CO(OOH) ; 13= -PAN    ; 14= -O-    ; 15= R-COO-R
# 16 = HCO-O-R; 17= -CO(F) ; 18= -CO(Cl) ;  19= -CO(Br) ; 20= -CO(I)
# This list is used in alisig which provide sigma for a given group
# ======================================================================
def khydration(chem, yield_val, ndat, chemdat, ydat):
    # grbond and rjgrm variables
    bond = [[0 for _ in range(mxnode)] for _ in range(mxnode)]
    group = [' ' * mxlgr for _ in range(mxnode)]
    nring = 0
    dbflg = 0
    rjg = [[0 for _ in range(2)] for _ in range(mxring)]  # ring-join group pairs

    # local
    tchemdat = [' ' for _ in range(len(chemdat))]
    tydat = [0.0 for _ in range(len(ydat))]
    tndat = 0
    node = 0
    tempkc = ' ' * len(chem)
    tgroup = [' ' * mxlgr for _ in range(mxnode)]
    pold = ''
    pnew = ''
    tbond = [[0 for _ in range(mxnode)] for _ in range(mxnode)]
    sigma = 0.0
    hsigma = 0.0
    khyd = 0.0
    sumy = 0.0
    numy = 0

    # chemmap variables 
    nodetype = [' ' for _ in range(mxnode)]
    alifun = [0.0 for _ in range(21)]
    cdfun = [0.0 for _ in range(21)]
    arofun = [0.0 for _ in range(21)]
    mapfun = [[[0.0 for _ in range(21)] for __ in range(3)] for ___ in range(mxnode)]
    funflg = [0 for _ in range(mxnode)]
    tabester = [[0 for _ in range(2)] for _ in range(4)]  # 1= -O- side, 2= CO side
    ngrp = 0
    nfcd = 0
    nfcr = 0
    ierr = 0

    # abcde_map variables 
    mxdeep = 7  # sigma effect scroll up to position 7 (para group in aromatic)
    nabcde = [0 for _ in range(mxdeep)]
    tabcde = [[[0 for _ in range(mxnode)] for __ in range(mxcp)] for ___ in range(mxdeep)]

    progname = 'khydration '
    mesg = ''

    nfcr = 0
    nfcd = 0
    ndat = 0
    tndat = 0
    # 初始化数组（左闭右闭范围）
    for i in range(1, len(chemdat) + 1):
        chemdat[i - 1] = ' '
        tchemdat[i - 1] = ' '
    for i in range(1, len(ydat) + 1):
        ydat[i - 1] = 0.0
        tydat[i - 1] = 0.0

    # build the group and bond matrix for chem      
    grbond(chem, group, bond, dbflg, nring)
    for i in range(1, mxnode + 1):
        bond[i - 1][i - 1] = 0
    rjgrm(nring, group, rjg)  # rm ring index and get ring closing nodes

    # make copy 
    for i in range(1, mxnode + 1):
        tgroup[i - 1] = group[i - 1]
    for i in range(1, mxnode + 1):
        for j in range(1, mxnode + 1):
            tbond[i - 1][j - 1] = bond[i - 1][j - 1]

    # count the number of group in chem and get the map of the functional 
    # group (see index number related the organic function)
    node = 0
    for i in range(1, mxnode + 1):
        if group[i - 1][0] != ' ':
            node += 1
    chemmap(chem, node, group, bond, ngrp, nodetype,
            alifun, cdfun, arofun, mapfun, funflg,
            tabester, nfcd, nfcr, ierr)
    if ierr != 0:
        mesg = "--error--, after a chemmap call "
        stoperr(progname, mesg, chem)

    # ---------------------------------
    # hydration of aliphatic aldehyde 
    # ---------------------------------
    if (alifun[8] != 0) or (cdfun[8] != 0):  # alifun(9)对应索引8，cdfun(9)对应索引8

        # find aldehyde position（原alidloop循环，Python用for+continue实现）
        for i in range(1, node + 1):
            # 检查当前节点是否为脂肪族/烯烃醛基（mapfun(i,1,9)或mapfun(i,2,9)）
            if not ((mapfun[i - 1][0][8] == 1.0) or (mapfun[i - 1][1][8] == 1.0)):
                continue
            
            sigma = 0.0
            abcde_map(bond, i, node, nabcde, tabcde)  # get the neighboors.

            # check that the aldehyde is not an hidden anhydride (-CO-O-CO-) 
            is_anhydride = False
            for k in range(1, nabcde[1] + 1):  # nabcde(2)对应索引1
                if nodetype[tabcde[1][k - 1][1] - 1] == 'o':
                    is_anhydride = True
                    break
            if is_anhydride:
                continue  # 对应原CYCLE alidloop，跳过当前节点

            # Remove the path that would lead to count the ester twice. 
            kill_ester(tabester, nabcde, tabcde)

            # get taft sigma (look for groups up to delta position)
            sigma = get_tsig(bond, group, 5, nabcde, tabcde, mapfun,
                             funflg, nodetype, sigma)

            # compute the hydration constant for aliphatic aldehyde 
            khyd = 0.0818 + 1.2663 * sigma

            # If the species is a 'Cd' aldehyde, then add the correction applied
            # for aromatic (no data available)
            if mapfun[i - 1][1][8] == 1.0:
                khyd = khyd - 1.5817

            # make the hydrate, rebuild, rename and reset tgroup & tbond
            pold = 'CHO '
            pnew = 'CH(OH)(OH) '
            swap(group[i - 1], pold, tgroup[i - 1], pnew)
            rebond(tbond, tgroup, tempkc, nring)
            stdchm(tempkc)
            # reset tgroup & tbond
            for j in range(1, mxnode + 1):
                tgroup[j - 1] = group[j - 1]
            for j in range(1, mxnode + 1):
                for k in range(1, mxnode + 1):
                    tbond[j - 1][k - 1] = bond[j - 1][k - 1]

            # store the data 
            tndat += 1
            if tndat > len(ydat):
                mesg = "# exceed max hydrate "
                stoperr(progname, mesg, chem)
            tydat[tndat - 1] = khyd + yield_val
            tchemdat[tndat - 1] = tempkc

    # ---------------------------------
    # hydration of aliphatic ketone 
    # ---------------------------------
    if alifun[9] != 0:  # alifun(10)对应索引9

        # find ketone position（原alikloop循环，Python用for循环+continue实现跳转）
        for i in range(1, node + 1):
            # 检查当前节点是否为脂肪族酮基（mapfun(i,1,10)）
            if mapfun[i - 1][0][9] != 1.0:
                continue  # 不是酮基，直接跳过当前节点
            
            sigma = 0.0
            # get the neighboors.
            abcde_map(bond, i, node, nabcde, tabcde)

            # check that the ketone is not an hidden anhydride (-CO-O-CO-) 
            is_anhydride = False
            for k in range(1, nabcde[1] + 1):  # nabcde(2)对应索引1
                # tabcde(2,k,2)对应tabcde[1][k-1][1]，检查相邻节点是否为-O-
                if nodetype[tabcde[1][k - 1][1] - 1] == 'o':
                    is_anhydride = True
                    break
            if is_anhydride:
                continue  # 是隐藏酸酐，跳过当前节点后续处理（对应原CYCLE alikloop）

            # Remove the path that would lead to count the ester twice. 
            # Performed by "killing" the "second" node of the functionality.
            kill_ester(tabester, nabcde, tabcde)

            # get taft sigma (look for groups up to delta position)
            sigma = get_tsig(bond, group, 5, nabcde, tabcde, mapfun,
                             funflg, nodetype, sigma)

            # compute the hydration constant for aliphatic ketone 
            khyd = 0.0818 + 1.2663 * sigma - 2.5039

            # make the hydrate, rebuild, rename and reset tgroup & tbond
            pold = 'CO '
            pnew = 'C(OH)(OH) '
            swap(group[i - 1], pold, tgroup[i - 1], pnew)
            rebond(tbond, tgroup, tempkc, nring)
            stdchm(tempkc)
            # reset tgroup & tbond
            for j in range(1, mxnode + 1):
                tgroup[j - 1] = group[j - 1]
            for j in range(1, mxnode + 1):
                for k in range(1, mxnode + 1):
                    tbond[j - 1][k - 1] = bond[j - 1][k - 1]

            # store the data 
            tndat += 1
            if tndat > len(ydat):
                mesg = "# exceed max hydrate "
                stoperr(progname, mesg, chem)
            tydat[tndat - 1] = khyd + yield_val
            tchemdat[tndat - 1] = tempkc

    # ---------------------------------
    # hydration of aromatic aldehyde 
    # ---------------------------------
    if arofun[8] != 0:  # arofun(9)对应索引8

        # find aldehyde position
        for i in range(1, node + 1):
            if mapfun[i - 1][2][8] == 1.0:  # mapfun(i,3,9)对应mapfun[i-1][2][8]
                sigma = 0.0

                abcde_map(bond, i, node, nabcde, tabcde)  # get the neighboors

                # Remove the path that would lead to count the ester twice. 
                # Performed by "killing" the "second" node of the functionality.
                kill_ester(tabester, nabcde, tabcde)

                # get hammet sigma (look for groups up to delta position)
                hsigma = get_hsig(nabcde, tabcde, mapfun, funflg, nodetype, hsigma)

                # compute the hydration constant for alihpatic aldehyde 
                khyd = 0.0818 + 0.4969 * hsigma - 1.5817

                # make the hydrate, rebuild, rename and reset tgroup & tbond
                pold = 'CHO '
                pnew = 'CH(OH)(OH) '
                swap(group[i - 1], pold, tgroup[i - 1], pnew)
                rebond(tbond, tgroup, tempkc, nring)
                stdchm(tempkc)
                # reset tgroup & tbond
                for j in range(1, mxnode + 1):
                    tgroup[j - 1] = group[j - 1]
                for j in range(1, mxnode + 1):
                    for k in range(1, mxnode + 1):
                        tbond[j - 1][k - 1] = bond[j - 1][k - 1]

                # store the data 
                tndat += 1
                if tndat > len(ydat):
                    mesg = "# exceed max hydrate "
                    stoperr(progname, mesg, chem)
                tydat[tndat - 1] = khyd + yield_val
                tchemdat[tndat - 1] = tempkc

    # ----------------------------------------------
    # collapse identical product 
    # ----------------------------------------------
    # note : must add khyd, not the log of it
    for i in range(1, tndat + 1):
        numy = 1
        sumy = 10 ** (tydat[i - 1])
        for j in range(i + 1, tndat + 1):
            if tchemdat[i - 1] == tchemdat[j - 1]:
                numy += 1
                sumy += 10 ** (tydat[j - 1])
                tchemdat[j - 1] = ' '
        if numy > 1:
            tydat[i - 1] = math.log10(sumy)

    for i in range(1, tndat + 1):
        if tchemdat[i - 1][0] != ' ':
            ndat += 1
            chemdat[ndat - 1] = tchemdat[i - 1]
            ydat[ndat - 1] = tydat[i - 1]

    return ndat, chemdat, ydat

# ======================================================================
# Purpose: Return the taft sigma value at for the positions given by the 
# of the functional group given in "tabcde". Note : the position for 
# which sigma is computed is given by tabcde(1,1,1).
#
# - ndeep: deepest position for which sigma must be considered (ndeep=5 
#   is delta position)
# - nabcde(k): number of distinct pathways that end up at a position k
#     relative to top (e.g. nabcde(3) gives the number of distinct 
#     pathways finishing in a beta position relative to top.
# - tabcde(k,i,j): give the pathways (nodes j), for the track number i
#     to reach the position k (k=2 is beta position ...). For example, 
#     tabcde(4,1,j) give the first track to reach a gamma position (node
#     given by tabcde(4,1,4), using the track given by tabcde(4,1,*).
# - mapfun(a,b,c): provide the number of function of type 'c' at position
#     (node) 'a'. index 'b' if for node type with 1=aliphatic, 2=cd and
#     3=aromatic. For example, the molecule CH2(OH)CdH=CdHCHO should 
#     have non zero values at the positions : mapfun(1,1,1)=1 and 
#     mapfun(4,2,9)=1
# - funflg(a): get the number of functional group at node a. For the 
#     example above, non-zero values are found at position 1 and 4, 
#     where it is set to 1.
# - nodetype: table of character for type node:
#      'y' = carbonyl      'r' = aromatic        'o'= -O- node
#      'd' = Cd            'n' = others (i.e. normal)
# - sigma: is the taft sigma due to the neighboors relative to the
#     position being considered (and given by tabcde(1,1,1).
#
# table of functions - Index of functionalities (in mapfun)
# -----------------------------------------------------------------
#  1= -OH    ;  2= -NO2     ;  3= -ONO2   ; 4= -OOH    ;  5= -F     
#  6= -Cl    ;  7=  -Br     ;  8= -I      ; 9= -CHO    ; 10= -CO-   
# 11= -COOH  ; 12= -CO(OOH) ; 13= -PAN    ; 14= -O-    ; 15= R-COO-R
# 16 = HCO-O-R; 17= -CO(F) ; 18= -CO(Cl) ;  19= -CO(Br) ; 20= -CO(I)
# 21 = -CO(O-)
# This list is used in alisig which provide sigma for a given group
# ======================================================================
def get_tsig(bond, group, ndeep, nabcde, tabcde, mapfun, funflg, nodetype, sigma):
    # 声明并初始化局部变量
    idp = 0
    k = 0
    ncd = 0
    locn = 0
    l = 0
    i1 = 0
    i2 = 0
    ifun = 0
    dist = 0.0
    tsig = 0.0
    aroflg = 0

    # Tafta sigma values for aliphatic compounds. sigma for ester depends
    # whether ester is connected from the -CO- side or the -O- side. 
    # Sigma=2.56 for -O- side and 2.00 for CO side.
    # CMV 16/06/14 : add sigma taft for CO(O-), n°21
    alisig = [0.62, 1.47, 1.38, 0.62, 1.10, 0.94, 1.00, 1.00, 2.15, 1.81,
              2.08, 2.08, 2.00, 1.81, 2.56, 2.90, 2.44, 2.37, 2.37, 2.37,
              -1.06]
    sigester = 2.00  # sigma for ester, CO side

    if ndeep > len(nabcde):
        print(f'--error-- in get_tsig, ndeep>SIZE(nabcde)={len(nabcde)}')
        raise Exception("in get_tsig")
    
    sigma = 0.0
    # scroll the neighboors up to ndeep position 
    for idp in range(2, ndeep + 1):  # idp从2到ndeep（左闭右闭）
        for k in range(1, nabcde[idp - 1] + 1):  # nabcde(idp)对应索引idp-1
            ncd = 0
            locn = tabcde[idp - 1][k - 1][idp - 1]  # tabcde(idp,k,idp)对应tabcde[idp-1][k-1][idp-1]
            if locn == 0:
                continue
            if funflg[locn - 1] != 0:  # funflg(locn)对应funflg[locn-1]

                # start the loop to see if Cd are in between
                for l in range(1, idp):  # l从1到idp-1（左闭右闭）
                    i1 = tabcde[idp - 1][k - 1][l - 1]  # tabcde(idp,k,l)对应tabcde[idp-1][k-1][l-1]
                    if i1 == 0:
                        continue
                    i2 = tabcde[idp - 1][k - 1][l]  # tabcde(idp,k,l+1)对应tabcde[idp-1][k-1][l]
                    if i2 == 0:
                        continue
                    # bond(i1,i2)对应bond[i1-1][i2-1]
                    if bond[i1 - 1][i2 - 1] == 2:
                        ncd += 1
                    # nodetype(i1)对应nodetype[i1-1]，nodetype(i2)对应nodetype[i2-1]
                    if (nodetype[i1 - 1] == 'd') and (nodetype[i2 - 1] == 'd'):
                        ncd += 1
                dist = float(idp) - float(ncd) - 2.0
                # group(locn)对应group[locn-1]
                if group[locn - 1][:3] == '-O-':
                    dist = dist - 1.0  # decrease the distance
   
                aroflg = 0  # check for aromatic nodes between groups
                for l in range(2, idp + 1):  # l从2到idp（左闭右闭）
                    # tabcde(idp,k,l)对应tabcde[idp-1][k-1][l-1]
                    if tabcde[idp - 1][k - 1][l - 1] == 0:
                        continue
                    node_val = tabcde[idp - 1][k - 1][l - 1]
                    if nodetype[node_val - 1] == 'r':
                        aroflg = 1

                # find the function on node locn (local node)
                for ifun in range(1, 22):  # ifun从1到21（左闭右闭）
                    # mapfun(locn,1,ifun)对应mapfun[locn-1][0][ifun-1]
                    if mapfun[locn - 1][0][ifun - 1] != 0:
                        tsig = mapfun[locn - 1][0][ifun - 1] * alisig[ifun - 1]
                        if ifun == 15:  # check ester
                            # group(locn)对应group[locn-1]
                            if group[locn - 1][:2] == 'CO':
                                tsig = sigester  # CO instead of -O-
                        sigma += (tsig * (0.4 ** dist))

                    # mapfun(locn,2,ifun)对应mapfun[locn-1][1][ifun-1]
                    if mapfun[locn - 1][1][ifun - 1] != 0:  # =CdC(X)CHO be must counted
                        tsig = mapfun[locn - 1][1][ifun - 1] * alisig[ifun - 1]
                        if ifun == 15:  # check ester
                            if group[locn - 1][:2] == 'CO':
                                tsig = sigester  # CO instead of -O-
                        sigma += (tsig * (0.4 ** dist))

                    # mapfun(locn,3,ifun)对应mapfun[locn-1][2][ifun-1]
                    if mapfun[locn - 1][2][ifun - 1] != 0:  # Phi-COCHO must be counted
                        if aroflg == 0:
                            tsig = mapfun[locn - 1][2][ifun - 1] * alisig[ifun - 1]
                            if ifun == 15:  # check ester
                                if group[locn - 1][:2] == 'CO':
                                    tsig = sigester  # CO instead of -O-
                            sigma += (tsig * (0.4 ** dist))

    return sigma

# ======================================================================
# Purpose: Return the hammet sigma (see header for get_tsig above for)
# group numbers. arorto, arometa, aropara give the hammet sigma for 
# each group in o,m,p.
# ======================================================================
def get_hsig(nabcde, tabcde, mapfun, funflg, nodetype, hsigma):
    # 声明并初始化局部变量
    idp = 0
    k = 0
    locn = 0
    l = 0
    ifun = 0
    funindex = 0
    dist = 0.0

    # CMV - 16/06/14. Add hammet values for -CO(O-)
    # benzoic as a reference for substituents in ortho
    arometa = [0.13, 0.74, 0.55, 0.00, 0.34, 0.37, 0.39, 0.35, 0.36, 0.36,
               0.35, 0.00, 0.00, 0.11, 0.32, 0.00, 0.55, 0.53, 0.53, 0.53,
               0.09]
    aropara = [-0.38, 0.78, 0.70, 0.00, 0.06, 0.24, 0.22, 0.21, 0.44, 0.47,
               0.44, 0.00, 0.00, -0.28, 0.39, 0.00, 0.70, 0.69, 0.69, 0.69,
               -0.05]
    arorto = [1.22, 1.99, 0.00, 0.00, 0.93, 1.28, 1.35, 1.34, 0.72, 0.07,
              0.95, 0.00, 0.00, 0.12, 0.63, 0.00, 0.00, 0.00, 0.00, 0.00,
              -0.91]
    
    hsigma = 0.0
    funindex = 0
    # scroll the neighboors up to 7 position (i.e. a carbonyl on para)
    for idp in range(2, 8):  # idp从2到7（左闭右闭）
        for k in range(1, nabcde[idp - 1] + 1):  # nabcde(idp)对应索引idp-1

            # exit if already taken into account by another pathway (it happens at 
            # every para position)
            if k >= 2:
                duplicate = False
                for l in range(1, k):  # l从1到k-1（左闭右闭）
                    # tabcde(idp,k,idp)对应tabcde[idp-1][k-1][idp-1]
                    # tabcde(idp,l,idp)对应tabcde[idp-1][l-1][idp-1]
                    if tabcde[idp - 1][k - 1][idp - 1] == tabcde[idp - 1][l - 1][idp - 1]:
                        duplicate = True
                        break
                if duplicate:
                    continue

            locn = tabcde[idp - 1][k - 1][idp - 1]  # tabcde(idp,k,idp)对应tabcde[idp-1][k-1][idp-1]
            if funflg[locn - 1] != 0:  # funflg(locn)对应funflg[locn-1]

                # find the function on node locn (local node)
                funindex = 0
                for ifun in range(1, 22):  # ifun从1到21（左闭右闭）
                    # CMV - 20/06/16. Look for function on aromatic AND aliphatic carbons
                    # mapfun(locn,3,ifun)对应mapfun[locn-1][2][ifun-1]
                    # mapfun(locn,1,ifun)对应mapfun[locn-1][0][ifun-1]
                    if (mapfun[locn - 1][2][ifun - 1] != 0) or (mapfun[locn - 1][0][ifun - 1] != 0):
                        funindex = ifun
            
                # CMV - 20/06/16. funindex should not be 0 at this stage
                if funindex == 0:
                    print("error, funindex = 0 in get_hsig")
                    raise Exception("in get_hsig")

                # find the distance (ortho, meta, para) on node locn (local node)
                dist = 0
                for l in range(1, idp + 1):  # l从1到idp（左闭右闭）
                    # tabcde(idp,k,l)对应tabcde[idp-1][k-1][l-1]
                    node_val = tabcde[idp - 1][k - 1][l - 1]
                    if node_val == 0:
                        continue
                    # nodetype(node_val)对应nodetype[node_val-1]
                    if nodetype[node_val - 1] == 'r':
                        dist += 1

                if dist == 2:
                    hsigma += arorto[funindex - 1]
                elif dist == 3:
                    hsigma += arometa[funindex - 1]
                elif dist == 4:
                    hsigma += aropara[funindex - 1]
                # tabcde(idp,k,2)对应tabcde[idp-1][k-1][1]
                if tabcde[idp - 1][k - 1][1] == 0:
                    hsigma = 0.0  # kill second node on ester

    return hsigma

# ======================================================================
# Purpose: Remove the path that would lead to count the ester twice. 
# Performed by "killing" the "second" node of the functionality.
#
# tabester: provide the position of ester "couple" (i.e. the O and CO 
# nodes. For example, the molecule CH3CO-O-CH2-O-COCH3 has the following
# values: tabester(1,1)=3,tabester(1,2)=2 
#         tabester(2,1)=5,tabester(2,2)=6
# nabcde(k): number of distinct pathways that end up at a position k 
# relative to top (e.g. nabcde(3) gives the number of distinct pathways
# finishing in a beta position relative to top.                   
# tabcde(k,i,j) : give the pathways (node j), for the track number i  
# to reach the position k (k=2 is beta position ...). For example, 
# tabcde(4,1,j) give the first track to reach a gamma position (node 
# given by tabcde(4,1,4), using the track given by tabcde(4,1,*).
# ======================================================================
def kill_ester(tabester, nabcde, tabcde):
    # 声明并初始化局部变量
    k = 0
    ipos = 0
    it = 0
    j = 0

    for k in range(1, len(tabester) + 1):  # k从1到SIZE(tabester,1)（左闭右闭）
        # tabester(k,1)对应tabester[k-1][0]
        if tabester[k - 1][0] != 0:

            for ipos in range(2, len(nabcde) + 1):  # ipos从2到SIZE(nabcde)（左闭右闭）
                # nabcde(ipos)对应nabcde[ipos-1]
                for it in range(1, nabcde[ipos - 1] + 1):  # it从1到nabcde(ipos)（左闭右闭）
                    for j in range(1, ipos):  # j从1到ipos-1（左闭右闭）

                        # tabcde(ipos,it,j)对应tabcde[ipos-1][it-1][j-1]
                        # tabcde(ipos,it,j+1)对应tabcde[ipos-1][it-1][j]
                        # tabester(k,1)对应tabester[k-1][0]，tabester(k,2)对应tabester[k-1][1]
                        if tabcde[ipos - 1][it - 1][j - 1] == tabester[k - 1][0]:
                            if tabcde[ipos - 1][it - 1][j] == tabester[k - 1][1]:
                                tabcde[ipos - 1][it - 1][j] = 0  # set the node to 0

                        if tabcde[ipos - 1][it - 1][j - 1] == tabester[k - 1][1]:
                            if tabcde[ipos - 1][it - 1][j] == tabester[k - 1][0]:
                                tabcde[ipos - 1][it - 1][j] = 0  # set the node to 0

    return tabcde