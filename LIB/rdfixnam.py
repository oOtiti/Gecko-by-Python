# 前置库已由用户提供，此处仅按需求导入对应模块内容
from keyparameter import mxldi, mxlco, mxlfo, tfu1, dirgecko
from normchem import stdchm
from database import mxfn, nfn, namfn, chemfn
from searching import srh5
from sortstring import sort_string

#*****************************************************************
# PURPOSE: Read the file prescibing (forcing) the name of some
# particular species. Needed to systematically set the same name
# for a given species.
#                                                                 
# Data are store in the fixed name (fn) database:                                                         
#   - nfn       : total nb. of species having a fixed name   
#   - namfn(i)  : table of the fixed name (6 character)       
#   - chemfn(i) : formula corresponding the species  withfixed name                               
#*****************************************************************
def rdfixnam():
    # 声明并初始化变量（对应Fortran声明，结合初始化）
    i = 0
    ilin = 0
    ierr = 0
    filename = ''
    line = ' ' * mxldi
    tempnam = [' ' * mxlco for _ in range(mxfn)]
    tempchem = [' ' * mxlfo for _ in range(mxfn)]

    # 初始化变量和数组（左闭右闭范围）
    global nfn  # 声明使用全局变量nfn
    nfn = 0
    for idx in range(len(namfn)):
        namfn[idx] = ' '
    for idx in range(len(chemfn)):
        chemfn[idx] = ' '
    for idx in range(len(tempnam)):
        tempnam[idx] = ' '
    for idx in range(len(tempchem)):
        tempchem[idx] = ' '

    filename = f"{dirgecko.strip()}/DATA/fixedname.dat"
    try:
        f = open(filename, 'r')
    except Exception as e:
        print(f'--error-- in rdfixnam: failed to open file {filename.strip()}')
        raise Exception("in rdfixnam, while opening file") from e

    # 对应原rdloop循环
    while True:
        try:
            line = f.readline()
            if not line:  # 文件结束
                raise EOFError("End of file reached without 'END' keyword")
        except EOFError as e:
            print('--error-- in rdfixnam. Missing keyword END ?')
            raise Exception("in rdfixnam, while reading inputs") from e
        except Exception as e:
            print('--error-- in rdfixnam. Error reading file')
            raise Exception("in rdfixnam, while reading inputs") from e
        
        # 处理行内容（去除换行符，补全到mxldi长度）
        line = line.rstrip('\n\r').ljust(mxldi)
        
        if line[0] == '!':
            continue
        if line[:3] == 'END':
            break

        nfn += 1
        if nfn >= mxfn:
            print('--error--, in rdfixnam: nfn > mxfn')
            raise Exception("in rdfixnam")
        
        # 读取固定名称和对应化学式（按空格分割，匹配原READ逻辑）
        try:
            # 按空格分割，第一个字段为tempnam（固定长度），剩余为tempchem
            parts = line.strip().split(maxsplit=1)
            if len(parts) < 2:
                raise ValueError("Insufficient fields")
            # 确保tempnam长度符合mxlco，超出截断，不足补空格
            tempnam[nfn - 1] = parts[0][:mxlco].ljust(mxlco)
            # 确保tempchem长度符合mxlfo，超出截断，不足补空格
            tempchem[nfn - 1] = parts[1][:mxlfo].ljust(mxlfo)
        except Exception as e:
            print(f'--error-- in rdfixnam, while reading data: {line.strip()}')
            raise Exception("in rdfixnam, while reading inputs") from e
        
        stdchm(tempchem[nfn - 1])

    f.close()

    # sort the formula
    for idx in range(len(chemfn)):
        chemfn[idx] = tempchem[idx]
    sort_string(chemfn[:nfn])  # 对应chemfn(1:nfn)

    # check for duplicate 
    ierr = 0
    if nfn > 1:
        for i in range(1, nfn):  # i=1到nfn-1（左闭右闭）
            if chemfn[i - 1] == chemfn[i]:
                print(f'--error--, while reading file: {filename.strip()}')
                print(f'Following species identified 2 times: {chemfn[i - 1].strip()}')
                ierr = 1
    if ierr != 0:
        raise Exception("in rdfixnam")
    
    # sort the tables according to formula
    for i in range(1, nfn + 1):  # i=1到nfn（左闭右闭）
        ilin = srh5(tempchem[i - 1], chemfn, nfn)
        if ilin <= 0:
            print(f'--error--, while sorting species in: {filename.strip()}')
            print(f'species "lost" after sorting the list: {tempchem[i - 1].strip()}')
            raise Exception("in rdfixnam")
        namfn[ilin - 1] = tempnam[i - 1]