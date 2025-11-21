import keyparameter
import rjtool
import stdgrbond
import toolbox

# ======================================================================
# Return the # of carbon atoms in a molecule                    
# ======================================================================
def cnum(chem):
    ic = 0
    nc = chem.find(' ')
    if nc == -1:
        nc = len(chem)
    else:
        nc += 1
    for i in range(1, nc + 1):
        idx = i - 1
        if idx < len(chem) and chem[idx] == 'C':
            if idx + 1 < len(chem) and chem[idx:idx+2] != 'Cl':
                ic += 1
        if idx < len(chem) and chem[idx] == 'c':
            ic += 1
    return ic

# ======================================================================
# Return the # of C-O-C (ether) in a molecule.
# (e.g. to estimate the # of of nodes in chem  
# ======================================================================
def onum(chem):
    io = 0
    nc = chem.find(' ')
    if nc == -1:
        nc = len(chem)
    else:
        nc -= 1
    for i in range(1, nc + 1):
        idx = i - 1
        if idx + 1 < len(chem) and chem[idx:idx+2] == '-O':
            io += 1
    return io

# ======================================================================
# PURPOSE: compute the molecular weight of the species provided as input
# ======================================================================
def molweight(chem):
    ica, iha, ina, ioa, ira, isa, ifl, ibr, icl = [0], [0], [0], [0], [0], [0], [0], [0], [0]
    getatoms(chem, ica, iha, ina, ioa, ira, isa, ifl, ibr, icl)
    
    weight = (ica[0] * 12.01 + iha[0] * 1.01 + ina[0] * 14.01 + ioa[0] * 16.0 + isa[0] * 32.07 +
              ifl[0] * 19.0 + ibr[0] * 79.9 + icl[0] * 35.45)
    return weight

# ======================================================================
# Find the number of each possible atom in the species CHEM,
# e.g. C's, H's, N's, O's, S's, F's, Cl's, Br's, and .'s  
# ======================================================================
def getatoms(chem, ica, iha, ina, ioa, ira, isa, ifl, ibr, icl):
    tchem = chem
    rjs = [[0 for _ in range(2)] for _ in range(keyparameter.mxring)]
    loring = False
    nc = tchem.find(' ')
    if nc == -1:
        nc = len(tchem)
    
    ifl[0] = 0
    ibr[0] = 0
    icl[0] = 0
    ica[0] = 0
    iha[0] = 0
    ina[0] = 0
    ioa[0] = 0
    isa[0] = 0
    ira[0] = 0
    n = 1

    loring = ('C1' in tchem) or ('Cd1' in tchem) or ('c1' in tchem) or ('-O1' in tchem)
    if loring:
        nring = 2
        rjtool.rjsrm(nring, tchem, rjs)

    i = nc
    while True:
        i -= 1
        if i < 1:
            break
        idx = i - 1

        if idx < len(tchem):
            if tchem[idx] == '2':
                n = 2
            elif tchem[idx] == '3':
                n = 3
            elif tchem[idx] == '4':
                n = 4

        if i > 1:
            prev_idx = idx - 1
            if prev_idx >= 0 and idx + 1 <= len(tchem):
                if tchem[prev_idx:idx+1] == 'Cl':
                    icl[0] += n
                    n = 1
                    i -= 1

        if idx < len(tchem):
            if tchem[idx] == 'F':
                ifl[0] += n
                n = 1
            elif idx + 1 <= len(tchem) and tchem[idx:idx+2] == 'Br':
                ibr[0] += n
                n = 1
            elif tchem[idx] == 'H':
                iha[0] += n
                n = 1
            elif tchem[idx] == 'N':
                ina[0] += n
                n = 1
            elif tchem[idx] == 'O':
                ioa[0] += n
                n = 1
            elif tchem[idx] == 'S':
                isa[0] += 1
                n = 1
            elif tchem[idx] == '.':
                ira[0] += n
                n = 1
            elif (idx + 1 <= len(tchem) and tchem[idx:idx+2] != 'Cl') and \
                 (tchem[idx] == 'C' or tchem[idx] == 'c'):
                ica[0] += 1
                n = 1

# ======================================================================
def oxnum(chem, numoxC):
    noC = [0 for _ in range(keyparameter.mxnode)]
    bond = [[0 for _ in range(keyparameter.mxnode)] for _ in range(keyparameter.mxnode)]
    dbflg = [0]
    nring = [0]
    group = [' ' * keyparameter.mxlgr for _ in range(keyparameter.mxnode)]
    rjg = [[0 for _ in range(2)] for _ in range(keyparameter.mxring)]
    nC = 0
    ngr = 0

    numoxC[0] = 0.0
    for idx in range(len(noC)):
        noC[idx] = 0

    stdgrbond.grbond(chem, group, bond, dbflg, nring)
    ngr = sum(1 for g in group if g.strip() != '')
    rjtool.rjgrm(nring[0], group, rjg)

    for i in range(ngr):
        if group[i][0] != 'C' and group[i][0] != 'c':
            continue
        nC += 1

        if group[i][:3] == 'CH3':
            noC[i] = -3
        elif group[i][:3] == 'CH2':
            noC[i] = -2
        elif group[i][:2] == 'CH':
            noC[i] = -1
        elif group[i][:2] == 'cH':
            noC[i] = -1
        elif group[i][:3] == 'CHO':
            noC[i] = 1
        elif group[i][:2] == 'CO':
            noC[i] = 2

        nfunc = toolbox.countstring(group[i], '(O')
        if nfunc != 0:
            noC[i] += 1
        nfunc = toolbox.countstring(group[i], '(N')
        if nfunc != 0:
            noC[i] += 1

        for j in range(ngr):
            if bond[i][j] == 3:
                noC[i] += 1

    if nC != 0:
        numoxC[0] = sum(noC[:ngr]) / float(nC)
    else:
        numoxC[0] = 0.0