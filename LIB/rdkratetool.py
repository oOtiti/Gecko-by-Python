# 前置库已由用户提供，此处仅按需求导入对应模块内容
from keyparameter import mxlfo, dirgecko, tfu1
from database import (
    mxkdb,
    nkohdb, kohdb_chem, kohdb_298, kohdb_arr, kohdb_com,
    nko3db, ko3db_chem, ko3db_298, ko3db_arr, ko3db_com,
    nkno3db, kno3db_chem, kno3db_298, kno3db_arr, kno3db_com
)
from sortstring import sort_string
from searching import srh5
from normchem import stdchm
from references import nreac_info, code

# ======================================================================
# Purpose: 
# ======================================================================
def rdkoxfiles():
    # 初始化变量（对应Fortran声明，保持初值一致）
    filename = ''
    
    print('  ...reading OH rate constants')
    filename = f"{dirgecko.strip()}/DATA/koh_rate.dat"
    rdratedb(filename, nkohdb, kohdb_chem, kohdb_298, kohdb_arr, kohdb_com)

    print('  ...reading O3 rate constants')
    filename = f"{dirgecko.strip()}/DATA/ko3_rate.dat"
    rdratedb(filename, nko3db, ko3db_chem, ko3db_298, ko3db_arr, ko3db_com)

    print('  ...reading NO3 rate constants')
    filename = f"{dirgecko.strip()}/DATA/kno3_rate.dat"
    rdratedb(filename, nkno3db, kno3db_chem, kno3db_298, kno3db_arr, kno3db_com)

# ======================================================================
# Purpose: read a database file for a given set of reactions (e.g. 
# VOC+OH, VOC+O3, ...)
# ======================================================================
def rdratedb(filename, ndat, oxchem, k298, arrh, oxcom):
    # 声明并初始化变量（对应Fortran声明，结合初始化）
    line = ''  # string in the db
    i = 0
    j = 0
    n1 = 0
    n2 = 0
    ierr = 0
    ilin = 0
    ipc = [0] * 7  # expected # of ";" char in string
    maxdat = 0  # max # of data allowed
    tpcom = ' ' * len(oxcom[0][0])  # a temporary comment
    # temporary copy of the output tables (need to sort tables)
    tpoxchem = [' ' * len(oxchem[0]) for _ in range(len(oxchem))]
    tpk298 = [0.0 for _ in range(len(arrh))]
    tparrh = [[0.0 for _ in range(len(arrh[0]))] for __ in range(len(arrh))]
    tpoxcom = [[' ' * len(oxcom[0][0]) for _ in range(len(oxcom[0]))] for __ in range(len(oxcom))]
    
    # initialize
    maxdat = len(arrh)
    ndat = 0
    for idx in range(len(k298)):
        k298[idx] = 0.0
    for idx1 in range(len(arrh)):
        for idx2 in range(len(arrh[0])):
            arrh[idx1][idx2] = 0.0
    for idx in range(len(oxchem)):
        oxchem[idx] = ' '
    for idx1 in range(len(oxcom)):
        for idx2 in range(len(oxcom[0])):
            oxcom[idx1][idx2] = ' '
    
    # 复制到临时数组
    for idx in range(len(tpoxchem)):
        tpoxchem[idx] = oxchem[idx]
    for idx in range(len(tpk298)):
        tpk298[idx] = k298[idx]
    for idx1 in range(len(tparrh)):
        for idx2 in range(len(tparrh[0])):
            tparrh[idx1][idx2] = arrh[idx1][idx2]
    for idx1 in range(len(tpoxcom)):
        for idx2 in range(len(tpoxcom[0])):
            tpoxcom[idx1][idx2] = oxcom[idx1][idx2]
    
    # open the file
    try:
        f = open(filename, 'r')
    except Exception as e:
        print(f'--error--, failed to open file: {filename}')
        raise Exception("in rdratedb") from e
    
    # read the file
    ilin = 0
    while True:
        ilin += 1
        try:
            line = f.readline()
            if not line:  # 文件结束
                raise EOFError("End of file reached without 'END' keyword")
        except EOFError as e:
            print(f'--error--, while reading file: {filename.strip()}')
            print(f'at line number: {ilin}')
            print('keyword "END" missing ?')
            raise Exception("in rdratedb") from e
        except Exception as e:
            print(f'--error--, while reading file: {filename.strip()}')
            print(f'at line number: {ilin}')
            raise Exception("in rdratedb") from e
        
        # 处理行内容（去除换行符）
        line = line.rstrip('\n\r')
        
        if line.startswith('!'):
            continue
        if line.startswith('END'):
            break
        
        # remove anything after "!" (not necessary as 1st character)
        n1 = line.find("!")
        if n1 > 0:
            line = line[:n1]
        
        # get the position of the ";" separating the field
        n1 = 0
        for i in range(1, len(line.strip()) + 1):
            if line[i-1:i] == ';':
                n1 += 1
                ipc[n1-1] = i  # Python列表0-based，ipc(1)对应ipc[0]
            if n1 == 7:
                break
        if n1 != 7:
            print(f'--error--, while reading file: {filename.strip()}')
            print(f'at line number: {ilin}')
            print(f'missing ";" separator. n1= {n1}')
            raise Exception("in rdratedb")
        if ';' in line[ipc[6]:]:  # ipc(7)对应ipc[6]
            print(f'--error--, while reading file: {filename.strip()}')
            print(f'at line number: {ilin}')
            print('too many ";" separator. n1 > 7')
            raise Exception("in rdratedb")
        
        # read the species
        ndat += 1
        if ndat > maxdat:
            print(f'--error--, while reading file: {filename.strip()}')
            print(f'too many species (check table size). maxdat= {maxdat}')
            raise Exception("in rdratedb")
        try:
            # line(1:ipc(1)-1)对应line[0:ipc[0]-1]
            tpoxchem[ndat-1] = line[:ipc[0]-1].strip()
        except Exception as e:
            print(f'--error--, while reading file: {filename.strip()}')
            print(f'at line number: {ilin}')
            print(f'while reading species: {line[:ipc[0]-1]}')
            raise Exception("in rdratedb") from e
        
        # read the rate constant @ 298K
        try:
            # line(ipc(1)+1:ipc(2)-1)对应line[ipc[0]:ipc[1]-1]
            tpk298[ndat-1] = float(line[ipc[0]:ipc[1]-1].strip())
        except Exception as e:
            print(f'--error--, while reading file: {filename.strip()}')
            print(f'while reading k298 at: {line[ipc[0]:ipc[1]-1].strip()}')
            print(f'at line number: {ilin}')
            raise Exception("in rdratedb") from e
        
        # read arrhenius coefficients
        for i in range(1, 4):  # i=1,2,3
            try:
                # line(ipc(i+1)+1:ipc(i+2)-1)对应line[ipc[i]:ipc[i+1]-1]
                tparrh[ndat-1][i-1] = float(line[ipc[i]:ipc[i+1]-1].strip())
            except Exception as e:
                print(f'--error--, while reading file: {filename.strip()}')
                print(f'while reading rate coef. at: {line[ipc[i]:ipc[i+1]-1].strip()}')
                print(f'at line number: {ilin}')
                raise Exception("in rdratedb") from e
        
        # read the comments
        for i in range(1, 4):  # i=1,2,3
            # line(ipc(i+4)+1:)对应line[ipc[i+3]:]
            tpcom = line[ipc[i+3]:].strip().ljust(len(tpcom))
            n2 = tpcom.find(';')
            if n2 != -1:
                tpcom = tpcom[:n2].ljust(len(tpcom))  # remove string not part of comment
            if tpcom[0] == ' ':
                continue  # nothing to read
            try:
                tpoxcom[ndat-1][i-1] = tpcom.strip()
            except Exception as e:
                print(f'--error--, while reading file: {filename.strip()}')
                print(f'while reading comment (code): {tpcom.strip()}')
                print(f'at line number: {ilin}')
                print(f'{line.strip()}')
                raise Exception("in rdratedb") from e
    
    f.close()

    # check and normalize the formula
    for i in range(1, ndat + 1):
        stdchm(tpoxchem[i-1])
    
    # sort the formula
    for i in range(len(oxchem)):
        oxchem[i] = tpoxchem[i]
    sort_string(oxchem[:ndat])  # 传递前ndat个元素（对应oxchem(1:ndat)）

    # check for duplicate 
    ierr = 0
    for i in range(1, ndat):  # i=1到ndat-1（左闭右闭）
        if oxchem[i-1] == oxchem[i]:
            print(f'--error--, while reading file: {filename.strip()}')
            print(f'Following species identified 2 times: {oxchem[i-1].strip()}')
            ierr = 1
    if ierr != 0:
        raise Exception("in rdratedb")
    
    # sort the tables according to formula
    for i in range(1, ndat + 1):
        ilin = srh5(tpoxchem[i-1], oxchem, ndat)
        if ilin <= 0:
            print(f'--error--, while sorting species in: {filename.strip()}')
            print(f'species "lost" after sorting the list: {tpoxchem[i-1].strip()}')
            raise Exception("in rdratedb")
        k298[ilin-1] = tpk298[i-1]
        for j in range(len(arrh[0])):
            arrh[ilin-1][j] = tparrh[i-1][j]
        for j in range(len(oxcom[0])):
            oxcom[ilin-1][j] = tpoxcom[i-1][j]
    
    # check the availability of the codes
    ierr = 0
    for i in range(1, ndat + 1):
        for j in range(1, len(oxcom[0]) + 1):
            if oxcom[i-1][j-1].strip() != '':
                n1 = srh5(oxcom[i-1][j-1], code, nreac_info)
                if n1 <= 0:
                    print(f'--error--, while reading file {filename.strip()}')
                    print(f'unknow comment: {oxcom[i-1][j-1].strip()}')
                    ierr = 1
    if ierr != 0:
        raise Exception("in rdratedb")