import sys
from keyparameter import tfu1, dirgecko
from references import mxreac_info, mxlcod, mxlreac_info, nreac_info, code, reac_info
from sortstring import sort_string
from searching import srh5

def rdreac_info():
    global nreac_info, code, reac_info
    
    filename = dirgecko + 'DATA/references.dat'
    try:
        tfu1_file = open(filename, 'r')
    except IOError:
        print('--error--, while trying to open: ', filename)
        sys.exit("in rdreac_info")
    
    nreac_info = 0
    tpcode = [''] * mxreac_info
    tpreac_info = [''] * mxreac_info
    
    ilin = 0
    while True:
        line = tfu1_file.readline()
        if not line:
            print('--error--, while reading file: ', filename)
            print('at line number: ', ilin)
            print('keyword "END" missing ?')
            sys.exit("in rdreac_info")
            
        ilin += 1
        line = line.strip()
        if not line or line[0] == '!':
            continue
            
        if line[0:3] == 'END':
            break
            
        n1 = line.find(':')
        if n1 <= 0:
            print('--error--, while reading file: ', filename)
            print('at line number: ', ilin)
            print('Expect ":" separator in line: ', line)
            sys.exit("in rdreac_info")
            
        leftline = line[0:n1].strip()
        rightline = line[n1+1:].strip()
        
        lcode = len(leftline)
        if lcode > mxlcod:
            print('--error--, while reading file: ', filename)
            print('at line number: ', ilin)
            print('Length of the code exceed mxlcod in: ', line)
            sys.exit("in rdreac_info")
            
        cblank = leftline.find(' ')
        if cblank != -1:
            print('--error--, while reading file: ', filename)
            print('at line number: ', ilin)
            print('Unexpected " " char in code: ', line)
            sys.exit("in rdreac_info")
            
        nreac_info += 1
        if nreac_info > mxreac_info:
            print('--error--, while reading file: ', filename)
            print('too many reac_infos (check table size). mxreac_info= ', mxreac_info)
            sys.exit("in rdreac_info")
            
        tpcode[nreac_info-1] = leftline.ljust(mxlcod)[:mxlcod]
        
        if len(rightline) > mxlreac_info:
            print('--error--, while reading file: ', filename)
            print('at line number: ', ilin)
            print('Length of the code exceed mxlreac_info in: ', line)
            sys.exit("in rdreac_info")
            
        tpreac_info[nreac_info-1] = rightline.ljust(mxlreac_info)[:mxlreac_info]
        
    tfu1_file.close()
    
    code = tpcode.copy()
    # Fortran did: CALL sort_string(code(1:nreac_info))
    # i.e. sort only the first nreac_info entries and leave the rest unchanged.
    sub_code = code[:nreac_info]
    sort_string(sub_code)
    code[:nreac_info] = sub_code
    
    ierr = 0
    for i in range(nreac_info-1):
        if code[i] == code[i+1]:
            print('--error--, while reading file: ', filename)
            print('Following code identified 2 times: ', code[i].strip())
            ierr = 1
            
    if ierr != 0:
        sys.exit("in rdreac_info")
        
    for i in range(nreac_info):
        ilin = srh5(tpcode[i], code, nreac_info)
        if ilin <= 0:
            print('--error--, while sorting reac_info in: ', filename)
            print('reac_info "lost" after sorting the list: ', tpcode[i].strip())
            sys.exit("in rdreac_info")
            
        reac_info[ilin-1] = tpreac_info[i]

def fullref(shortcode):
    from references import mxreac_info, code, reac_info, mxlreac_info
    
    fullref_str = ' ' * mxlreac_info
    for i in range(mxreac_info):
        if code[i].strip() == shortcode.strip():
            fullref_str = reac_info[i]
            break
            
    return fullref_str