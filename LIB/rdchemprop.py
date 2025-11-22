from database import nhlcdb, hlcdb_chem, hlcdb_dat, hlcdb_com, nkhydb, khydb_chem, khydb_dat, khydb_com
from keyparameter import dirgecko, tfu1
from sortstring import sort_string
from searching import srh5
from normchem import stdchm
from references import nreac_info, code


def rdhenry():
    nhlcdb = 0
    for i in range(len(hlcdb_chem)):
        hlcdb_chem[i] = ''
    for i in range(len(hlcdb_dat)):
        for j in range(len(hlcdb_dat[0])):
            hlcdb_dat[i][j] = 0.0
    for i in range(len(hlcdb_com)):
        for j in range(len(hlcdb_com[0])):
            hlcdb_com[i][j] = ''
    filename = dirgecko.strip() + 'DATA/henry_const.dat'
    rddatabase(filename, 3, nhlcdb, hlcdb_chem, hlcdb_dat, hlcdb_com)


def rdhydrat():
    nkhydb = 0
    for i in range(len(khydb_chem)):
        khydb_chem[i] = ''
    for i in range(len(khydb_dat)):
        for j in range(len(khydb_dat[0])):
            khydb_dat[i][j] = 0.0
    for i in range(len(khydb_com)):
        for j in range(len(khydb_com[0])):
            khydb_com[i][j] = ''
    filename = dirgecko.strip() + 'DATA/hydrat_const.dat'
    rddatabase(filename, 1, nkhydb, khydb_chem, khydb_dat, khydb_com)


def rddatabase(filename, ncol, ndat, chemframe, dataframe, comframe):
    maxdat = len(dataframe)
    ndat = 0
    for i in range(len(dataframe)):
        for j in range(len(dataframe[0])):
            dataframe[i][j] = 0.0
    for i in range(len(chemframe)):
        chemframe[i] = ''
    for i in range(len(comframe)):
        for j in range(len(comframe[0])):
            comframe[i][j] = ''
    tpchemframe = chemframe.copy()
    tpdataframe = [row[:] for row in dataframe]
    tpcomframe = [row[:] for row in comframe]
    nsep = ncol + 3
    
    tfu1 = open(filename, 'r')
    
    ilin = 0
    while True:
        ilin = ilin + 1
        line = tfu1.readline()
        if not line:
            print('--error--, while reading file: ', filename.strip())
            print('at line number: ', ilin)
            print('keyword "END" missing ?')
            raise Exception("in rddatabase")
        if line[0] == '!':
            continue
        if line[0:3] == 'END':
            break
        if '!' in line:
            line = line[0:line.index('!')]
        
        n1 = 0
        ipc = [0] * (nsep + 1)
        for i in range(len(line.strip())):
            if line[i] == ';':
                n1 = n1 + 1
                ipc[n1] = i
            if n1 == nsep:
                break
        if n1 != nsep:
            print('--error--, while reading file: ', filename.strip())
            print('at line number: ', ilin)
            print('missing ";" separator. n1= ', n1)
            raise Exception("in rddatabase")
        if ';' in line[ipc[nsep]:]:
            print('--error--, while reading file: ', filename.strip())
            print('at line number: ', ilin)
            print('too many ";" separator. n1 > ', nsep)
            raise Exception("in rddatabase")
        
        ndat = ndat + 1
        if ndat > maxdat:
            print('--error--, while reading file: ', filename.strip())
            print('too many species (check table size). maxdat= ', maxdat)
            raise Exception("in rddatabase")
        tpchemframe[ndat-1] = line[0:ipc[1]].strip()
        
        for i in range(1, ncol+1):
            try:
                tpdataframe[ndat-1][i-1] = float(line[ipc[i]+1:ipc[i+1]])
            except:
                print('--error--, while reading file: ', filename.strip())
                print('while reading rate coef. at: ', line[ipc[i]+1:ipc[i+1]])
                print('at line number: ', ilin)
                raise Exception("in rddatabase")
        
        for i in range(1, 4):
            tpcom = line[ipc[i+ncol]+1:].strip()
            n2 = tpcom.find(';')
            if n2 != -1:
                tpcom = tpcom[0:n2]
            if len(tpcom) == 0 or tpcom[0] == ' ':
                continue
            try:
                tpcomframe[ndat-1][i-1] = tpcom
            except:
                print('--error--, while reading file: ', filename.strip())
                print('while reading comment (code): ', tpcom)
                print('at line number: ', ilin)
                print(line.strip())
                raise Exception("in rddatabase")
                
    tfu1.close()
    
    for i in range(1, ndat+1):
        stdchm(tpchemframe[i-1])
    
    chemframe[0:ndat] = tpchemframe[0:ndat]
    sort_string(chemframe[0:ndat])
    
    ierr = 0
    for i in range(1, ndat):
        if chemframe[i-1] == chemframe[i]:
            print('--error--, while reading file: ', filename.strip())
            print('Following species identified 2 times: ', chemframe[i-1].strip())
            ierr = 1
    if ierr != 0:
        raise Exception("in rddatabase")
    
    for i in range(1, ndat+1):
        ilin = srh5(tpchemframe[i-1], chemframe, ndat)
        if ilin <= 0:
            print('--error--, while sorting species in: ', filename.strip())
            print('species "lost" after sorting the list: ', tpchemframe[i-1].strip())
            raise Exception("in rddatabase")
        for j in range(len(dataframe[0])):
            dataframe[ilin-1][j] = tpdataframe[i-1][j]
        for j in range(len(comframe[0])):
            comframe[ilin-1][j] = tpcomframe[i-1][j]
    
    ierr = 0
    for i in range(1, ndat+1):
        for j in range(1, len(comframe[0])+1):
            if len(comframe[i-1][j-1]) > 0 and comframe[i-1][j-1][0] != ' ':
                n1 = srh5(comframe[i-1][j-1], code, nreac_info)
                if n1 <= 0:
                    print('--error--, while reading file ', filename.strip())
                    print('unknown comment: ', comframe[i-1][j-1], " ", n1)
                    ierr = 1
    if ierr != 0:
        raise Exception("in rddatabase")