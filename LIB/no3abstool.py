def rabsno3(tbond, tgroup, ig, arrhc):
    # PURPOSE : Find rate constants for H-atom abstraction by NO3 based on
    # Kerdouci et al., chemphyschem, 3909, 2010 and Atmos. Env., 363, 2014.

    # 保持原Fortran循环范围（1-based遍历），仅对数组访问时-1转换索引
    arrhc[:] = 0.0

    ngr = sum(1 for g in tgroup if g.strip() != '')

    nca = 0
    saturated = True

    # 原循环范围：i=1到ngr
    for i in range(1, ngr + 1):
        # 数组访问转换为0-based
        if tgroup[i-1][:1] == 'C':
            nca += 1
        if tgroup[i-1][:2] == 'Cd':
            saturated = False

    # 检查ig值有效性（原逻辑：ig>ngr报错）
    if ig > ngr:
        print('--error--, in rabsno3, => ig is greater than ngr')
        raise SystemExit("in rabsno3")

    # FIND K(0) VALUE
    # ------------------
    # Values are from Kerdouci et al., 2011
    # ig为1-based，数组访问时-1
    if tgroup[ig-1][:3] == 'CH3':
        arrhc[0] = 1.00e-18
        arrhc[1] = 0.0
        arrhc[2] = 0.0
    elif tgroup[ig-1][:3] == 'CH2':
        arrhc[0] = 2.56e-17
        arrhc[1] = 0.0
        arrhc[2] = 0.0
    elif tgroup[ig-1][:3] == 'CHO':
        # 原循环范围：i=1到ngr
        for i in range(1, ngr + 1):
            if tbond[ig-1, i-1] == 3:
                arrhc[:] = 0.0
                return
        arrhc[0] = 2.415e-15
        arrhc[1] = 0.0
        arrhc[2] = 0.0
    elif tgroup[ig-1][:2] == 'CH':
        arrhc[0] = 1.05e-16
        arrhc[1] = 0.0
        arrhc[2] = 0.0
    else:
        print(f'--error-- in rabsno3, no NO3 reaction for: {tgroup[ig-1].strip()}')
        raise SystemExit("in rabsno3")

    # ---------------------------------------
    # FIND MULTIPLIERS BASED ON SUBSTITUENTS
    # ---------------------------------------

    # on same carbon:
    mult = 1.0
    if '(OH)' in tgroup[ig-1]:
        mult = 18.0
        arrhc[0] *= mult

    if 'CHO' in tgroup[ig-1] and nca > 3 and saturated:
        # 使用math.exp替代numpy.exp，避免额外库依赖
        import math
        mult = -14.5 + 21.4 * (1 - math.exp(-0.43 * nca))
        arrhc[0] *= mult

    # on alpha carbons:
    nether = 0
    # 原循环范围：i=1到ngr
    for i in range(1, ngr + 1):
        mult = 1.0
        # 数组访问转换为0-based
        if tbond[ig-1, i-1] != 0:
            # simple alkyl:
            if tgroup[i-1][:3] == 'CH3':
                mult = 1.0
                arrhc[0] *= mult

            if tgroup[i-1][:4] == 'CH2 ':
                mult = 1.02
                # 原循环范围：j=1到ngr
                for j in range(1, ngr + 1):
                    if tbond[i-1, j-1] == 3 and j != ig:
                        mult = 1.0
                        break
                arrhc[0] *= mult

            if tgroup[i-1][:4] == 'CH2(':
                mult = 1.02
                arrhc[0] *= mult

            if tgroup[i-1][:3] == 'CH ':
                # 原循环范围：j=1到ngr
                for j in range(1, ngr + 1):
                    if tbond[i-1, j-1] == 3 and j != ig:
                        mult = 10.0
                        break
                    else:
                        mult = 1.61
                arrhc[0] *= mult

            if tgroup[i-1][:3] == 'CH(':
                # 原循环范围：j=1到ngr
                for j in range(1, ngr + 1):
                    if tbond[i-1, j-1] == 3 and j != ig:
                        mult = 10.0
                        break
                    else:
                        mult = 1.61
                arrhc[0] *= mult

            if tgroup[i-1][:2] == 'C(':
                # 原循环范围：j=1到ngr
                for j in range(1, ngr + 1):
                    if tbond[i-1, j-1] == 3 and j != ig:
                        mult = 48.0
                        break
                    else:
                        mult = 2.03
                arrhc[0] *= mult

            if tgroup[i-1][:2] == 'C ':
                # 原循环范围：j=1到ngr
                for j in range(1, ngr + 1):
                    if tbond[i-1, j-1] == 3 and j != ig:
                        mult = 48.0
                        break
                    else:
                        mult = 2.03
                arrhc[0] *= mult

            # overwrite for carbonyls, alcohols and double bonds values
            if tgroup[i-1][:2] == 'CO':
                mult = 0.64
                arrhc[0] *= mult

            if tgroup[i-1][:3] == 'CHO':
                mult = 143.0
                arrhc[0] *= mult

            if tgroup[i-1][:2] == 'Cd':
                mult = 1.0
                arrhc[0] *= mult

            # Ethers and esters
            # For acetal, consider only one -O- influence
            if '-O-' in tgroup[i-1]:
                # 原循环范围：j=1到ngr
                for j in range(1, ngr + 1):
                    if tbond[i-1, j-1] == 3 and j != ig:
                        # ether in alpha of ig
                        if tgroup[j-1][:3] == 'CH3':
                            mult = 130.0
                            arrhc[0] *= mult
                        elif tgroup[j-1][:4] == 'CH2 ':
                            mult = 58.0
                            arrhc[0] *= mult
                        elif tgroup[j-1][:4] == 'CH2(':
                            mult = 58.0
                            arrhc[0] *= mult
                        elif tgroup[j-1][:3] == 'CH ':
                            mult = 23.0
                            arrhc[0] *= mult
                        elif tgroup[j-1][:3] == 'CH(':
                            mult = 23.0
                            arrhc[0] *= mult
                        elif tgroup[j-1][:2] == 'C(':
                            mult = 495.0
                            arrhc[0] *= mult
                        elif tgroup[j-1][:2] == 'C ':
                            mult = 495.0
                            arrhc[0] *= mult