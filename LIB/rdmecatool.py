# 前置库已由用户提供，此处仅按需求导入对应模块内容
from keyparameter import dirgecko, tfu1
from normchem import stdchm
from sortstring import sort_string
from searching import srh5
from references import nreac_info, code
from database import (
    nkwoh, nkwoh_pd, kwoh_rct, kwoh_pd, kwoh_copd, kwoh_yld, kwoh_com,
    nkwo3, nkwo3_pd, kwo3_rct, kwo3_pd, kwo3_copd, kwo3_yld, kwo3_com,
    nkwno3, nkwno3_pd, kwno3_rct, kwno3_pd, kwno3_copd, kwno3_yld, kwno3_com,
    nkwro2, kwro2_arrh, kwro2_stoi, kwro2_rct, kwro2_prd, kwro2_com,
    nkwrco3, kwrco3_arrh, kwrco3_stoi, kwrco3_rct, kwrco3_prd, kwrco3_com,
    nkwro, kwro_arrh, kwro_stoi, kwro_rct, kwro_prd, kwro_com,
    nkwcri, kwcri_arrh, kwcri_stoi, kwcri_rct, kwcri_prd, kwcri_com,
    njdat, jlabel, jchem, jprod, coprodj, nj40, jlab40, j40, j_com
)

# ======================================================================
# PURPOSE: read known mechanism and store the data in the 
# database module
# ======================================================================
def rdmeca():
    # 声明并初始化变量（对应Fortran声明，结合初始化）
    filename = ''

    # VOC+OH 
    print('  ...reading known VOC+OH reactions ...')
    filename = f"{dirgecko.strip()}/DATA/oh_prod.dat"
    rdvocmchdb(filename, nkwoh, nkwoh_pd, kwoh_rct, kwoh_pd, kwoh_copd, kwoh_yld, kwoh_com)

    # VOC+NO3 
    print('  ...reading known VOC+NO3 reactions ...')
    filename = f"{dirgecko.strip()}/DATA/no3_prod.dat"
    rdvocmchdb(filename, nkwno3, nkwno3_pd, kwno3_rct, kwno3_pd, kwno3_copd, kwno3_yld, kwno3_com)

    # VOC+O3 
    print('  ...reading known VOC+O3 reactions ...')
    filename = f"{dirgecko.strip()}/DATA/o3_prod.dat"
    rdvocmchdb(filename, nkwo3, nkwo3_pd, kwo3_rct, kwo3_pd, kwo3_copd, kwo3_yld, kwo3_com)

    # RO2 chemistry
    print('  ...reading known RO2 reactions ...')
    filename = f"{dirgecko.strip()}/DATA/ro2.dat"
    rdradmchdb(filename, nkwro2, kwro2_rct, kwro2_prd, kwro2_arrh, kwro2_stoi, kwro2_com)

    # RCOO2 chemistry
    print('  ...reading known RO2 reactions ...')
    filename = f"{dirgecko.strip()}/DATA/rcoo2.dat"
    rdradmchdb(filename, nkwrco3, kwrco3_rct, kwrco3_prd, kwrco3_arrh, kwrco3_stoi, kwrco3_com)

    # RO chemistry
    print('  ...reading known RO reactions ...')
    filename = f"{dirgecko.strip()}/DATA/ro.dat"
    rdradmchdb(filename, nkwro, kwro_rct, kwro_prd, kwro_arrh, kwro_stoi, kwro_com)

    # Criegee chemistry
    print('  ...reading known criegee reactions ...')
    filename = f"{dirgecko.strip()}/DATA/criegee.dat"
    rdradmchdb(filename, nkwcri, kwcri_rct, kwcri_prd, kwcri_arrh, kwcri_stoi, kwcri_com)

    # read known species for photolysis reactions
    rdhvdb()

# ======================================================================
# PURPOSE: Read mechanism data in the file provided as input. The
# routine is in particular called to read known VOC+OH, VOC+NO3, VOC+O3
# reactions.
# ======================================================================
def rdvocmchdb(filename, nkwdat, nkwpd, kwrct, kwpd, kwcopd, kwyld, kwcom):
    # 声明并初始化变量（对应Fortran声明，结合初始化）
    lenlin = 200  # length (max) of the line to be read
    line = ' ' * lenlin
    i = 0
    j = 0
    k = 0
    n1 = 0
    n2 = 0
    ndat = 0
    ncop = 0
    istart = 0
    iend = 0
    ierr = 0
    ilin = 0
    ncpbeg = 0
    ncpend = 0
    maxprd = 0
    maxcoprd = 0
    maxcom = 0
    maxent = 0

    # temporary copies of the table (before sorting)
    tpnkwpd = [0 for _ in range(len(nkwpd))]
    tpkwyld = [[0.0 for _ in range(len(kwyld[0]))] for __ in range(len(kwyld))]
    tpkwrct = [' ' * len(kwrct[0]) for _ in range(len(kwrct))]
    tpkwpd = [[' ' * len(kwpd[0][0]) for _ in range(len(kwpd[0]))] for __ in range(len(kwpd))]
    tpkwcopd = [[[' ' * len(kwcopd[0][0][0]) for _ in range(len(kwcopd[0][0]))] for __ in range(len(kwcopd[0]))] for ___ in range(len(kwcopd))]
    tpkwcom = [[' ' * len(kwcom[0][0]) for _ in range(len(kwcom[0]))] for __ in range(len(kwcom))]

    # 初始化变量和数组（nkwdat是标量，直接赋值）
    nkwdat = 0
    for idx in range(len(nkwpd)):
        nkwpd[idx] = 0
    for idx in range(len(kwrct)):
        kwrct[idx] = ' '
    for idx1 in range(len(kwpd)):
        for idx2 in range(len(kwpd[0])):
            kwpd[idx1][idx2] = ' '
    for idx1 in range(len(kwyld)):
        for idx2 in range(len(kwyld[0])):
            kwyld[idx1][idx2] = 0.0
    for idx1 in range(len(kwcopd)):
        for idx2 in range(len(kwcopd[0])):
            for idx3 in range(len(kwcopd[0][0])):
                kwcopd[idx1][idx2][idx3] = ' '
    for idx1 in range(len(kwcom)):
        for idx2 in range(len(kwcom[0])):
            kwcom[idx1][idx2] = ' '

    for idx in range(len(tpkwrct)):
        tpkwrct[idx] = ' '
    for idx1 in range(len(tpkwpd)):
        for idx2 in range(len(tpkwpd[0])):
            tpkwpd[idx1][idx2] = ' '
    for idx1 in range(len(tpkwyld)):
        for idx2 in range(len(tpkwyld[0])):
            tpkwyld[idx1][idx2] = 0.0
    for idx1 in range(len(tpkwcopd)):
        for idx2 in range(len(tpkwcopd[0])):
            for idx3 in range(len(tpkwcopd[0][0])):
                tpkwcopd[idx1][idx2][idx3] = ' '
    for idx1 in range(len(tpkwcom)):
        for idx2 in range(len(tpkwcom[0])):
            tpkwcom[idx1][idx2] = ' '

    maxent = len(nkwpd)
    maxprd = len(kwcopd[0])
    maxcoprd = len(kwcopd[0][0])
    maxcom = len(kwcom[0])

    # 打开文件
    try:
        f = open(filename, 'r')
    except Exception as e:
        print(f'--error--, failed to open file: {filename.strip()}')
        raise Exception("in rdvocmchdb") from e

    # read the data
    ilin = 0
    while True:
        ilin += 1
        try:
            line = f.readline()
            if not line:
                raise EOFError("End of file reached without 'END' keyword")
        except EOFError as e:
            print(f'--error--, while reading file: {filename.strip()}')
            print(f'at line number: {ilin}')
            print('keyword "END" missing ?')
            raise Exception("in rdvocmchdb") from e
        except Exception as e:
            print(f'--error--, while reading file: {filename.strip()}')
            print(f'at line number: {ilin}')
            raise Exception("in rdvocmchdb") from e

        # 处理行内容
        line = line.rstrip('\n\r').ljust(lenlin)

        if line[:3] == 'END':
            break
        if line[0] == '!':
            continue

        # check that the line is correctly formatted
        n1 = line.find(':')
        if n1 == -1:
            print(f'--error--, while reading file: {filename.strip()}')
            print(f'in rdvocmchdb, ":" not found at line:', ilin)
            raise Exception("in rdvocmchdb")

        ncpbeg = line.find(';')
        if ncpbeg == -1:
            print(f'--error--, while reading file: {filename.strip()}')
            print(f'in rdvocmchdb, ";" not found at line:', ilin)
            raise Exception("in rdvocmchdb")

        nkwdat += 1
        if nkwdat > maxent:
            print(f'--error--, while reading file: {filename.strip()}')
            print(f'tnumber of reactions is greater than mxkwr:', maxent)
            raise Exception("in rdvocmchdb")

        # read reactant(s) and the number of products to be read
        try:
            tpkwrct[nkwdat - 1] = line[:n1].strip().ljust(len(kwrct[0]))
        except Exception as e:
            print(f'--error--, while reading file: {filename.strip()}')
            print(f'while reading reactant at line number:', ilin)
            print(line.strip())
            raise Exception("in rdvocmchdb") from e

        try:
            ndat = int(line[n1+1:ncpbeg].strip())
        except Exception as e:
            print(f'--error--, while reading file: {filename.strip()}')
            print(f'while reading # of data at line number:', ilin)
            print(line.strip())
            raise Exception("in rdvocmchdb") from e

        if ndat > maxprd:
            print(f'--error--, while reading file: {filename.strip()}')
            print(f'at line number: {ilin}')
            print(f' number of products is greater than maxprd:', maxprd)
            raise Exception("in rdvocmchdb")

        # read the comment code
        for i in range(1, maxcom + 1):
            ncpend = line[ncpbeg+1:].find(';')
            if ncpend == -1:
                if line[ncpbeg+1:].strip() != "":
                    try:
                        tpkwcom[nkwdat - 1][i - 1] = line[ncpbeg+1:].strip().ljust(len(kwcom[0][0]))
                    except Exception as e:
                        print(f'--error--, while reading file: {filename.strip()}')
                        print(f'while reading reference code at line:', ilin)
                        print(line.strip())
                        raise Exception("in rdvocmchdb") from e
                break
            else:
                ncpend = ncpbeg + ncpend + 1  # count ncpend from the 1st char in line
                if line[ncpbeg+1:ncpend-1].strip() != "":
                    try:
                        tpkwcom[nkwdat - 1][i - 1] = line[ncpbeg+1:ncpend-1].strip().ljust(len(kwcom[0][0]))
                    except Exception as e:
                        print(f'--error--, while reading file: {filename.strip()}')
                        print(f'while reading reference code at line:', ilin)
                        print(line.strip())
                        raise Exception("in rdvocmchdb") from e
                ncpbeg = ncpend
                continue

        stdchm(tpkwrct[nkwdat - 1])
        tpnkwpd[nkwdat - 1] = ndat

        # loop over the products
        for j in range(1, ndat + 1):
            ilin += 1
            try:
                line = f.readline()
                if not line:
                    raise EOFError("Data missing")
            except EOFError as e:
                print(f'--error--, while reading file: {filename.strip()}')
                print(f'at line number: {ilin}')
                print('data missing ??? ')
                raise Exception("in rdvocmchdb") from e
            except Exception as e:
                print(f'--error--, while reading file: {filename.strip()}')
                print(f'at line number: {ilin}')
                raise Exception("in rdvocmchdb") from e

            line = line.rstrip('\n\r').ljust(lenlin)
            n1 = line.find(':')
            if n1 == -1:
                print(f'--error--, while reading file: {filename.strip()}')
                print(f'--error-- in rdvocmchdb, ":" not found at line:', ilin)
                raise Exception("in rdvocmchdb")

            # count the number of coproducts
            ncop = 0
            for k in range(n1 + 1, lenlin):
                if line[k:k+1] == '+':
                    ncop += 1
            if ncop > maxcoprd:
                print(f'--error--, while reading file: {filename.strip()}')
                print(f'at line number: {ilin}')
                print(f' number of coproducts is greater than maxcoprd:', maxcoprd)
                raise Exception("in rdvocmchdb")

            # read the yield and the 'main' product
            try:
                tpkwyld[nkwdat - 1][j - 1] = float(line[:n1].strip())
            except Exception as e:
                print(f'--error--, while reading file: {filename.strip()}')
                print(f'while reading yield at line number: {ilin}')
                print(line.strip())
                raise Exception("in rdvocmchdb") from e

            if ncop == 0:  # only 1 product
                try:
                    tpkwpd[nkwdat - 1][j - 1] = line[n1+2:].strip().ljust(len(kwpd[0][0]))
                except Exception as e:
                    print(f'--error--, while reading file: {filename.strip()}')
                    print(f'while reading product at line number: {ilin}')
                    print(line.strip())
                    raise Exception("in rdvocmchdb") from e
            else:  # read the product and coproducts
                n2 = line.find('+')
                try:
                    tpkwpd[nkwdat - 1][j - 1] = line[n1+2:n2-1].strip().ljust(len(kwpd[0][0]))
                except Exception as e:
                    print(f'--error--, while reading file: {filename.strip()}')
                    print(f'while reading product at line number: {ilin}')
                    print(line.strip())
                    raise Exception("in rdvocmchdb") from e
                istart = n2
                for k in range(1, ncop + 1):
                    iend = line[istart+1:].find('+')
                    if iend == -1:
                        iend = lenlin
                    else:
                        iend = istart + iend + 1
                    try:
                        tpkwcopd[nkwdat - 1][j - 1][k - 1] = line[istart+2:iend-1].strip().ljust(len(kwcopd[0][0][0]))
                    except Exception as e:
                        print(f'--error--, while reading file: {filename.strip()}')
                        print(f'while reading co-products at line number: {ilin}')
                        print(line.strip())
                        raise Exception("in rdvocmchdb") from e
                    istart = iend

            stdchm(tpkwpd[nkwdat - 1][j - 1])

            if ncop > 0:
                for k in range(1, ncop + 1):
                    if tpkwcopd[nkwdat - 1][j - 1][k - 1][0] == ' ':
                        print(f'--error--, while reading file: {filename.strip()}')
                        print(f'expected co-product not found in line #: {ilin}')
                        print(line.strip())
                        raise Exception("in rdvocmchdb")

    f.close()

    # sort the formula
    for idx in range(len(kwrct)):
        kwrct[idx] = tpkwrct[idx]
    sort_string(kwrct[:nkwdat])  # 对应kwrct(1:nkwdat)

    # check for duplicate 
    ierr = 0
    for i in range(1, nkwdat):  # i=1到nkwdat-1（左闭右闭）
        if kwrct[i - 1] == kwrct[i]:
            print(f'--error--, while reading file: {filename.strip()}')
            print(f'Following species identified 2 times: {kwrct[i - 1].strip()}')
            ierr = 1
    if ierr != 0:
        raise Exception("in rdvocmchdb")

    # sort the tables according to formula
    for i in range(1, nkwdat + 1):  # i=1到nkwdat（左闭右闭）
        ilin = srh5(tpkwrct[i - 1], kwrct, nkwdat)
        if ilin <= 0:
            print(f'--error--, while sorting species in: {filename.strip()}')
            print(f'species "lost" after sorting the list: {tpkwrct[i - 1].strip()}')
            raise Exception("in rdvocmchdb")
        nkwpd[ilin - 1] = tpnkwpd[i - 1]
        for j in range(len(kwyld[0])):
            kwyld[ilin - 1][j] = tpkwyld[i - 1][j]
        for j in range(len(kwpd[0])):
            kwpd[ilin - 1][j] = tpkwpd[i - 1][j]
        for j in range(len(kwcopd[0])):
            for k in range(len(kwcopd[0][0])):
                kwcopd[ilin - 1][j][k] = tpkwcopd[i - 1][j][k]
        for j in range(len(kwcom[0])):
            kwcom[ilin - 1][j] = tpkwcom[i - 1][j]

    # check the availability of the codes
    ierr = 0
    for i in range(1, nkwdat + 1):
        for j in range(1, maxcom + 1):
            if kwcom[i - 1][j - 1][0] != ' ':
                n1 = srh5(kwcom[i - 1][j - 1], code, nreac_info)
                if n1 <= 0:
                    print(f'--error--, while reading file {filename.strip()}')
                    print(f'unknow comment: {kwcom[i - 1][j - 1].strip()}')
                    ierr = 1
    if ierr != 0:
        raise Exception("in rdvocmchdb")

# ======================================================================
# PURPOSE: Read mechanism data in the file provided as input. The
# routine is in particular called to read know RO2, RCO3 and RO 
# chemistry.
#
# Data are returned to be next stored the database module:  
# - nkr : # of reaction given as input in the file
# - krct(:,2) : formula of the 2 reactants for reaction i
# - kprd(:,3) : formula the main products (up to 3) 
# - arrh(:,3) : arrhenius para. for the reactions         
# - kcost(:,3): stoi. coef. for the main products
# - kcom(:3)  : comment's code for the reactions
# ======================================================================
def rdradmchdb(filename, nkr, krct, kprd, arrh, kcost, kcom):
    # 声明并初始化变量（对应Fortran声明，结合初始化）
    i = 0
    j = 0
    n1 = 0
    n2 = 0
    n3 = 0
    n4 = 0
    n5 = 0
    n6 = 0
    n7 = 0
    n8 = 0
    n9 = 0
    n10 = 0
    nlast = 0
    ierr = 0
    loerr = 0
    ilin = 0
    ncpbeg = 0
    ncpend = 0
    lenlin = 300  # length (max) of the line to be read
    line = ' ' * lenlin
    maxent = 0
    maxcom = 0
    nkr = 0
    for idx1 in range(len(arrh)):
        for idx2 in range(len(arrh[0])):
            arrh[idx1][idx2] = 0.0
    for idx1 in range(len(kcost)):
        for idx2 in range(len(kcost[0])):
            kcost[idx1][idx2] = 0.0
    for idx1 in range(len(krct)):
        for idx2 in range(len(krct[0])):
            krct[idx1][idx2] = ' '
    for idx1 in range(len(kprd)):
        for idx2 in range(len(kprd[0])):
            kprd[idx1][idx2] = ' '
    for idx1 in range(len(kcom)):
        for idx2 in range(len(kcom[0])):
            kcom[idx1][idx2] = ' '

    maxent = len(arrh)
    maxcom = len(kcom[0])

    # 打开文件
    try:
        f = open(filename, 'r')
    except Exception as e:
        print(f'--error--, failed to open file: {filename.strip()}')
        raise Exception("in rdradmchdb") from e

    # read data
    ilin = 0
    while True:
        ilin += 1
        try:
            line = f.readline()
            if not line:
                raise EOFError("End of file reached without 'END' keyword")
        except EOFError as e:
            print(f'--error--, while reading file: {filename.strip()}')
            print(f'at line number: {ilin}')
            print('keyword "END" missing ?')
            raise Exception("in rdradmchdb") from e
        except Exception as e:
            print(f'--error--, while reading file: {filename.strip()}')
            print(f'at line number: {ilin}')
            raise Exception("in rdradmchdb") from e

        # 处理行内容
        line = line.rstrip('\n\r').ljust(lenlin)

        if line[:3] == 'END':
            break
        if line[0] == '!':
            continue

        # 计算|的位置（Fortran INDEX是1-based，Python find是0-based，需调整）
        n1 = line.find('|')
        if n1 == -1:
            n2 = -1
        else:
            n2 = line[n1+1:].find('|')
            n2 = n2 + n1 + 1 if n2 != -1 else -1
        if n2 == -1:
            n3 = -1
        else:
            n3 = line[n2+1:].find('|')
            n3 = n3 + n2 + 1 if n3 != -1 else -1
        if n3 == -1:
            n4 = -1
        else:
            n4 = line[n3+1:].find('|')
            n4 = n4 + n3 + 1 if n4 != -1 else -1
        if n4 == -1:
            n5 = -1
        else:
            n5 = line[n4+1:].find('|')
            n5 = n5 + n4 + 1 if n5 != -1 else -1
        if n5 == -1:
            n6 = -1
        else:
            n6 = line[n5+1:].find('|')
            n6 = n6 + n5 + 1 if n6 != -1 else -1
        if n6 == -1:
            n7 = -1
        else:
            n7 = line[n6+1:].find('|')
            n7 = n7 + n6 + 1 if n7 != -1 else -1
        if n7 == -1:
            n8 = -1
        else:
            n8 = line[n7+1:].find('|')
            n8 = n8 + n7 + 1 if n8 != -1 else -1
        if n8 == -1:
            n9 = -1
        else:
            n9 = line[n8+1:].find('|')
            n9 = n9 + n8 + 1 if n9 != -1 else -1
        if n9 == -1:
            nlast = -1
        else:
            nlast = line[n9+1:].find('|')
        n10 = n9 + nlast + 1 if (n9 != -1 and nlast != -1) else -1

        # check that the line is correctly formatted
        if nlast == -1:
            print(f'--error--, while reading: {filename.strip()}')
            print(f' missing "|" at line: {ilin}')
            print(line.strip())
            raise Exception("in rdradmchdb")

        nkr += 1
        if nkr >= maxent:
            print('--error--, while reading ro.dat')
            print('number of reaction greater than mkr')
            raise Exception("in rdradmchdb")

        # search for ';' (reference code)
        ncpbeg = line.find(';')
        if ncpbeg == -1:
            print(f'--error--, while reading file: {filename.strip()}')
            print(f'in rdvocmchdb, ";" not found at line:', ilin)
            raise Exception("in rdradmchdb")

        # read data
        loerr = 0
        try:
            krct[nkr - 1][0] = line[:n1].strip().ljust(len(krct[0][0]))
        except Exception as e:
            loerr += 1
        try:
            krct[nkr - 1][1] = line[n1+1:n2].strip().ljust(len(krct[0][0]))
        except Exception as e:
            loerr += 1
        try:
            kcost[nkr - 1][0] = float(line[n2+1:n3].strip())
        except Exception as e:
            loerr += 1
        try:
            kprd[nkr - 1][0] = line[n3+1:n4].strip().ljust(len(kprd[0][0]))
        except Exception as e:
            loerr += 1
        try:
            kcost[nkr - 1][1] = float(line[n4+1:n5].strip())
        except Exception as e:
            loerr += 1
        try:
            kprd[nkr - 1][1] = line[n5+1:n6].strip().ljust(len(kprd[0][0]))
        except Exception as e:
            loerr += 1
        try:
            kcost[nkr - 1][2] = float(line[n6+1:n7].strip())
        except Exception as e:
            loerr += 1
        try:
            kprd[nkr - 1][2] = line[n7+1:n8].strip().ljust(len(kprd[0][0]))
        except Exception as e:
            loerr += 1
        try:
            kcost[nkr - 1][3] = float(line[n8+1:n9].strip())
        except Exception as e:
            loerr += 1
        try:
            kprd[nkr - 1][3] = line[n9+1:n10].strip().ljust(len(kprd[0][0]))
        except Exception as e:
            loerr += 1
        try:
            arrh_parts = line[n10+1:].strip().split()
            for j in range(1, 4):
                arrh[nkr - 1][j - 1] = float(arrh_parts[j - 1])
        except Exception as e:
            loerr += 1

        if loerr != 0:
            print(f'--error--, while reading: {filename.strip()}')
            print(f'An error identified while reading line #: {ilin}')
            print(line.strip())
            raise Exception("in rdradmchdb")

        # 处理'- '为空字符串
        if krct[nkr - 1][1][:2] == '- ':
            krct[nkr - 1][1] = ' ' * len(krct[0][0])
        if kprd[nkr - 1][1][:2] == '- ':
            kprd[nkr - 1][1] = ' ' * len(kprd[0][0])
        if kprd[nkr - 1][2][:2] == '- ':
            kprd[nkr - 1][2] = ' ' * len(kprd[0][0])
        if kprd[nkr - 1][3][:2] == '- ':
            kprd[nkr - 1][3] = ' ' * len(kprd[0][0])

        # check the formula
        stdchm(krct[nkr - 1][0])
        stdchm(krct[nkr - 1][1])
        stdchm(kprd[nkr - 1][0])
        stdchm(kprd[nkr - 1][1])
        stdchm(kprd[nkr - 1][2])
        stdchm(kprd[nkr - 1][3])

        # read the comment code
        for i in range(1, maxcom + 1):
            ncpend = line[ncpbeg+1:].find(';')
            if ncpend == -1:
                if line[ncpbeg+1:].strip() != "":
                    try:
                        kcom[nkr - 1][i - 1] = line[ncpbeg+1:].strip().ljust(len(kcom[0][0]))
                    except Exception as e:
                        print(f'--error--, while reading file: {filename.strip()}')
                        print(f'while reading reference code at line:', ilin)
                        print(line.strip())
                        raise Exception("in rdradmchdb") from e
                break
            else:
                ncpend = ncpbeg + ncpend + 1  # count ncpend from the 1st char in line
                if line[ncpbeg+1:ncpend-1].strip() != "":
                    try:
                        kcom[nkr - 1][i - 1] = line[ncpbeg+1:ncpend-1].strip().ljust(len(kcom[0][0]))
                    except Exception as e:
                        print(f'--error--, while reading file: {filename.strip()}')
                        print(f'while reading reference code at line:', ilin)
                        print(line.strip())
                        raise Exception("in rdradmchdb") from e
                ncpbeg = ncpend
                continue

    f.close()

    # check the availability of the codes
    ierr = 0
    for i in range(1, nkr + 1):
        for j in range(1, maxcom + 1):
            if kcom[i - 1][j - 1][0] != ' ':
                n1 = srh5(kcom[i - 1][j - 1], code, nreac_info)
                if n1 <= 0:
                    print(f'--error--, while reading file {filename.strip()}')
                    print(f'unknown comment: {kcom[i - 1][j - 1].strip()}')
                    ierr = 1
    if ierr != 0:
        raise Exception("in rdradmchdb")

# ======================================================================
# PURPOSE: (1) read the data in photo.dat file (i.e the photolytic 
# data for the species used as "reference") and (2) read typical J 
# values (at a zenithal angle of 40).
#
# Data are stored in the database module:                    
# - njdat : # of photolytic reaction.    
# - jlabel(:) : ID # of photolytic reaction 
# - jchem(:)  : formula of the species being photolyzed  
# - jprod(:,2): the two main photodissociation fragments
# - coprodj(:): additional inorganic coproduct
# - j_com(:,:): comment's code for the photolysis reactions
# - nj40      : # of J40 data in database
# - jlab40(i) : label for which J4O is provided    
# - j40(i)    : J40 values for various labels      
# ======================================================================
def rdhvdb():
    # 声明并初始化变量（对应Fortran声明，结合初始化）
    i = 0
    j = 0
    n1 = 0
    n2 = 0
    n3 = 0
    n4 = 0
    ilin = 0
    ierr = 0
    loerr = 0
    idat = 0
    jdat = 0
    maxj = 0
    maxj40 = 0
    maxcom = 0
    ncpend = 0
    ncpbeg = 0
    lenlin = 300  # length (max) of the line to be read
    line = ' ' * lenlin

    # 初始化数组（njdat、nj40是标量，直接赋值）
    global njdat, nj40  # 声明使用全局标量变量
    njdat = 0
    nj40 = 0
    for idx in range(len(jchem)):
        jchem[idx] = ' '
    for idx in range(len(jlabel)):
        jlabel[idx] = 0
    for idx1 in range(len(jprod)):
        for idx2 in range(len(jprod[0])):
            jprod[idx1][idx2] = ' '
    for idx in range(len(coprodj)):
        coprodj[idx] = ' '
    for idx1 in range(len(j_com)):
        for idx2 in range(len(j_com[0])):
            j_com[idx1][idx2] = ' '
    for idx in range(len(j40)):
        j40[idx] = 0.0
    for idx in range(len(jlab40)):
        jlab40[idx] = 0

    jdat = 0
    idat = 0
    maxj = len(jlabel)
    maxj40 = len(jlab40)
    maxcom = len(j_com[0])

    # read the known reactions and their associated label
    # ---------------------------------------------------
    try:
        f = open(f"{dirgecko.strip()}/DATA/photo.dat", 'r')
    except Exception as e:
        print(f'--error--, failed to open file: photo.dat')
        raise Exception("in rdhvdb") from e

    # read photolysis data
    ilin = 0
    while True:
        ilin += 1
        try:
            line = f.readline()
            if not line:
                raise EOFError("End of file reached without 'END' keyword")
        except EOFError as e:
            print('--error--, while reading file: photo.dat')
            print(f'at line number: {ilin}')
            print('keyword "END" missing ?')
            raise Exception("in rdhvdb") from e
        except Exception as e:
            print('--error--, while reading file: photo.dat')
            print(f'at line number: {ilin}')
            raise Exception("in rdhvdb") from e

        # 处理行内容
        line = line.rstrip('\n\r').ljust(lenlin)

        if line[:3] == 'END':
            break
        if line[0] == '!':
            continue

        # 计算|的位置（Fortran INDEX是1-based，Python find是0-based，需调整）
        n1 = line.find('|')
        if n1 == -1:
            n2 = -1
        else:
            n2 = line[n1+1:].find('|')
            n2 = n2 + n1 + 1 if n2 != -1 else -1
        if n2 == -1:
            n3 = -1
        else:
            n3 = line[n2+1:].find('|')
            n3 = n3 + n2 + 1 if n3 != -1 else -1
        if n3 == -1:
            n4 = -1
        else:
            n4 = line[n3+1:].find('|')
            n4 = n3 + n4 + 1 if n4 != -1 else -1

        # check that the line is correctly formatted
        if n4 == -1:
            print('--error--, while reading photo.dat')
            print(f' missing "|" at line: {ilin}')
            print(line.strip())
            raise Exception("in rdhvdb")

        # search for ';' (reference code)
        ncpbeg = line.find(';')
        if ncpbeg == -1:
            print('--error--, while reading file: phot.dat')
            print(f'in rdvocmchdb, ";" not found at line:', ilin)
            raise Exception("in rdhvdb")

        jdat += 1
        if jdat >= maxj:
            print('--error--, while reading photo.dat')
            print('number of reaction greater than mkr')
            raise Exception("in rdhvdb")

        loerr = 0
        try:
            jchem[jdat - 1] = line[:n1].strip().ljust(len(jchem[0]))
        except Exception as e:
            loerr += 1
        try:
            jlabel[jdat - 1] = int(line[n1+1:n2].strip())
        except Exception as e:
            loerr += 1
        try:
            jprod[jdat - 1][0] = line[n2+1:n3].strip().ljust(len(jprod[0][0]))
        except Exception as e:
            loerr += 1
        try:
            jprod[jdat - 1][1] = line[n3+1:n4].strip().ljust(len(jprod[0][0]))
        except Exception as e:
            loerr += 1
        try:
            coprodj[jdat - 1] = line[n4+1:].strip().ljust(len(coprodj[0]))
        except Exception as e:
            loerr += 1

        if loerr != 0:
            print('--error--, while reading photo.dat ')
            print(f'An error identified while reading line #: {ilin}')
            print(line.strip())
            raise Exception("in rdhvdb")

        # check the formula (if C>1 only)
        stdchm(jchem[jdat - 1])
        stdchm(jprod[jdat - 1][0])
        stdchm(jprod[jdat - 1][1])

        if coprodj[jdat - 1][:2] == '- ':
            coprodj[jdat - 1] = ' ' * len(coprodj[0])

        # read the comment code
        for i in range(1, maxcom + 1):
            ncpend = line[ncpbeg+1:].find(';')
            if ncpend == -1:
                if line[ncpbeg+1:].strip() != "":
                    try:
                        j_com[jdat - 1][i - 1] = line[ncpbeg+1:].strip().ljust(len(j_com[0][0]))
                    except Exception as e:
                        print('--error--, while reading file: photo.dat')
                        print(f'while reading reference code at line:', ilin)
                        print(line.strip())
                        raise Exception("in rdhvdb") from e
                break
            else:
                ncpend = ncpbeg + ncpend + 1  # count ncpend from the 1st char in line
                if line[ncpbeg+1:ncpend-1].strip() != "":
                    try:
                        j_com[jdat - 1][i - 1] = line[ncpbeg+1:ncpend-1].strip().ljust(len(j_com[0][0]))
                    except Exception as e:
                        print('--error--, while reading file: photo.dat')
                        print(f'while reading reference code at line:', ilin)
                        print(line.strip())
                        raise Exception("in rdhvdb") from e
                ncpbeg = ncpend
                continue

    f.close()
    njdat = jdat  # 标量直接赋值，无索引

    # check the availability of the codes
    ierr = 0
    for i in range(1, jdat + 1):
        for j in range(1, maxcom + 1):
            if j_com[i - 1][j - 1][0] != ' ':
                n1 = srh5(j_com[i - 1][j - 1], code, nreac_info)
                if n1 <= 0:
                    print('--error--, while reading file: photo.dat')
                    print(f'unknown comment: {j_com[i - 1][j - 1].strip()}')
                    ierr = 1
    if ierr != 0:
        raise Exception("in rdradmchdb")

    # Read the photolysis rates for a 40° zenith angle
    # ---------------------------------------------------
    try:
        f = open(f"{dirgecko.strip()}/DATA/j40.dat", 'r')
    except Exception as e:
        print(f'--error--, failed to open file: j40.dat')
        raise Exception("in rdhvdb") from e

    ilin = 0
    while True:
        ilin += 1
        try:
            line = f.readline()
            if not line:
                raise EOFError("End of file reached without 'END' keyword")
        except EOFError as e:
            print('--error--, while reading file j40.dat')
            print(f'at line number: {ilin}')
            print('keyword "END" missing ?')
            raise Exception("in rdhvdb") from e
        except Exception as e:
            print('--error--, while reading file j40.dat')
            print(f'at line number: {ilin}')
            raise Exception("in rdhvdb") from e

        # 处理行内容
        line = line.rstrip('\n\r')

        if line[:3] == 'END':
            break
        if line[0] == '!':
            continue

        idat += 1
        if idat >= maxj40:
            print('--error--,  while reading j40.dat')
            print(' number of data is greater than maxj40')
            raise Exception("in rdhvdb")

        try:
            parts = line.strip().split()
            jlab40[idat - 1] = int(parts[0])
            j40[idat - 1] = float(parts[1])
        except Exception as e:
            print('--error--, while reading j40.dat')
            print(f'at line number: {ilin}')
            raise Exception("in rdhvdb") from e

    f.close()
    nj40 = idat  # 标量直接赋值，无索引