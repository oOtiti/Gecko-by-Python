import os
import keyparameter
import database
import rjtool
import stdgrbond

# =======================================================================
# PURPOSE: read heat of formation of the benson groups
# The data are stored in module database as:     
#   - nbson        : total number of benson group              
#   - bsongrp(:)   : table of benson group                     
#   - bsonval(:)   : heat of formation corresponding to group i
# =======================================================================
def rdbenson():
    ierr = 0
    line = ''

    # initialize
    database.nbson = 0
    for i in range(len(database.bsongrp)):
        database.bsongrp[i] = ' '
    for i in range(len(database.bsonval)):
        database.bsonval[i] = 0.0

    # open the file
    file_path = os.path.join(keyparameter.dirgecko, 'DATA', 'benson.dat')
    try:
        tfu1 = open(file_path, 'r')
    except Exception as e:
        print('--error-- in rdbenson: cannot open file', file_path)
        raise SystemExit("in rdbenson, while opening file") from e

    # read benson group
    while True:
        try:
            line = tfu1.readline()
            if not line:
                ierr = -1
                break
            line = line.rstrip('\n')
        except Exception as e:
            ierr = -1
            break

        if ierr != 0:
            print('--error-- in rdbenson. Missing keyword END ?')
            tfu1.close()
            raise SystemExit("in rdbenson, while reading inputs") from e

        if line[:3] == 'END':
            break
        if line[:1] == '*':
            continue

        if len(line) >= 25 and line[24] != ' ':
            print('--error--, reading bensongrp.dat. Group > 24 characters')
            print(line.strip())
            tfu1.close()
            raise SystemExit("in rdbenson")

        database.nbson += 1
        if database.nbson >= database.mxbg:
            print('--error--, while reading bensongrp.dat')
            print('number of benson grp is greater than mbg')
            tfu1.close()
            raise SystemExit("in rdbenson")

        try:
            grp_str = line[:24].rstrip()
            val_str = line[25:35].strip()
            database.bsongrp[database.nbson - 1] = grp_str.ljust(database.lbg)
            database.bsonval[database.nbson - 1] = float(val_str)
        except Exception as e:
            print('--error-- in rdbenson, while reading bensongrp.dat:')
            print(line.strip())
            tfu1.close()
            raise SystemExit("in rdbenson, while reading inputs") from e

    tfu1.close()

# =======================================================================
# PURPOSE: return the heat (enthalpie) formation for the chemical 
# provided as input based on the Benson group contribution method.
# =======================================================================
def heat(chem):
    base = '   '
    bengrp = ' ' * database.lbg
    tgroup = [' ' * keyparameter.mxlgr for _ in range(keyparameter.mxnode)]
    tempkg = ' ' * keyparameter.mxlgr
    tbond = [[0 for _ in range(keyparameter.mxnode)] for _ in range(keyparameter.mxnode)]
    nca = 0
    point = 0
    numo = 0
    nring = [0]
    dbflg = [0]
    bvalue = 0.0
    i = 0
    j = 0
    ic = 0
    init = 0
    ir = 0
    icd = 0
    io = 0
    ico = 0
    nc = 0
    ib = 0
    ie = 0
    fg = 0
    op = [0, 0, 0, 0]
    cp = [0, 0, 0, 0]
    icdot = 0
    icodot = 0
    ialko = 0
    ipero = 0
    check = 0
    rjg = [[0 for _ in range(2)] for _ in range(keyparameter.mxring)]
    wtflag = 0

    if wtflag != 0:
        print('*heat*')

    # initialize:
    heat_val = 0.0
    bvalue = 0.0
    point = 0
    numo = 0
    base = '   '
    bengrp = ' ' * database.lbg
    tempkg = ' ' * keyparameter.mxlgr
    for x in range(keyparameter.mxnode):
        for y in range(keyparameter.mxnode):
            tbond[x][y] = 0
    for x in range(keyparameter.mxnode):
        tgroup[x] = ' ' * keyparameter.mxlgr

    nc = chem.find(' ')
    if nc == -1:
        nc = 0
    else:
        nc -= 1
    if nc < 1:
        return heat_val

    stdgrbond.grbond(chem, tgroup, tbond, dbflg, nring)
    if nring[0] != 0:
        rjtool.rjgrm(nring[0], tgroup, rjg)

    # count the number of node in chem
    nca = 0
    for i in range(keyparameter.mxnode):
        if tgroup[i][0] != ' ':
            nca += 1

    # treat some one-carbon molecules as special cases:
    if nca < 2:
        if tgroup[0].strip() == 'CH2O':
            heat_val = -26.0
            return heat_val
        elif tgroup[0].strip() == 'CO' and tgroup[1][0] != 'C':
            heat_val = -26.4
            return heat_val
        elif tgroup[0].strip() == 'CH3O':
            heat_val = 4.62
            return heat_val
        elif tgroup[0].strip() == 'CH3.':
            heat_val = 34.9
            return heat_val
        elif tgroup[0].strip() == 'CHO.':
            heat_val = 9.99
            return heat_val

    # treat multinodes:
    for i in range(nca):
        tempkg = tgroup[i]
        if tempkg[:2].strip() == '':
            continue

        #  -O- node
        # -------------
        if tempkg[:3] == '-O-':
            ic = 0
            icdot = 0
            ico = 0
            io = 0
            icodot = 0
            init = 0
            ir = 0
            base = '   '

            point = 2
            base = 'O'

            for j in range(keyparameter.mxnode):
                ir = 0
                if tbond[i][j] != 0:
                    nc_j = tgroup[j].find(' ')
                    if nc_j == -1:
                        nc_j = len(tgroup[j])
                    else:
                        nc_j -= 1
                    if nc_j >= 0 and tgroup[j][nc_j] == '.':
                        ir = 1

                    if tgroup[j][:2] == 'CO' or tgroup[j][:3] == 'CHO':
                        if ir == 1:
                            icodot += 1
                        else:
                            ico += 1
                    else:
                        if ir == 1:
                            icdot += 1
                        else:
                            ic += 1

            ialko = 0
            ipero = 0
            icd = 0
            nameben(base, point, ic, icdot, io, ico, icodot, icd, init, ialko, ipero, bengrp)
            getben(bengrp, bvalue, check)
            if check != 0:
                print('--error-- in heat. Benson group not found:', bengrp.strip())
                print('for chemical: ', chem.strip())
                raise SystemExit("in heat")
            heat_val += bvalue

        # carbon center
        # -------------
        if tempkg[0] == 'C':
            ic = 0
            icdot = 0
            ico = 0
            io = 0
            icodot = 0
            init = 0
            ir = 0
            icd = 0
            ialko = 0
            ipero = 0
            base = '   '

            # treat some special groups.
            # ketene: DHf C2H2O = -47.6 => DHf ketene group = -47.6 - Cd_(Cd) [6.26] = -53.86
            # reference for ketene: Cohen, N. (1996) J.Phys.Chem.Ref.Data, 25(6), 1411-1481
            if tempkg[:3] == 'CdO':
                bvalue = -53.86
                heat_val += bvalue
                continue

            # set the base
            nc_temp = tempkg.find(' ')
            if nc_temp == -1:
                nc_temp = len(tempkg)
            else:
                nc_temp -= 1
            if nc_temp >= 0 and tempkg[nc_temp] == '.':
                ir = 1
            if tempkg[:2] == 'CO':
                if ir == 1:
                    point = 4
                    base = 'CO*'
                else:
                    point = 3
                    base = 'CO'
            elif tempkg[:3] == 'CHO':
                if ir == 1:
                    print('--error--, Benson group HCO* unexpected in: ', chem.strip())
                    raise SystemExit("in heat")
                else:
                    point = 3
                    base = 'CO'
            elif tempkg[:2] == 'Cd':
                if ir == 1:
                    print('--error--, Benson group Cd* unexpected in: ', chem.strip())
                    raise SystemExit("in heat")
                else:
                    point = 3
                    base = 'Cd'
            else:
                if ir == 1:
                    point = 3
                    base = 'C*'
                else:
                    point = 2
                    base = 'C'

            # find node-node bond :
            for j in range(keyparameter.mxnode):
                ir = 0
                if tbond[i][j] != 0:
                    nc_j = tgroup[j].find(' ')
                    if nc_j == -1:
                        nc_j = len(tgroup[j])
                    else:
                        nc_j -= 1
                    if nc_j >= 0 and tgroup[j][nc_j] == '.':
                        ir = 1

                    if tgroup[j][:2] == 'CO' or tgroup[j][:3] == 'CHO':
                        if ir == 1:
                            icodot += 1
                        else:
                            ico += 1
                    elif tgroup[j][:3] == '-O-':
                        io += 1
                    elif tgroup[j][:2] == 'Cd':
                        icd += 1
                    else:
                        if ir == 1:
                            icdot += 1
                        else:
                            ic += 1

            # find carbon - functional group bonds:
            # first the program finds the open and closing parenthesis (since functional
            # groups are included into parenthesis) and count number of groups.
            fg = 0
            op = [0, 0, 0, 0]
            cp = [0, 0, 0, 0]
            for j in range(point - 1, keyparameter.mxlgr):
                if tempkg[j] == '(':
                    fg += 1
                    if fg <= 4:
                        op[fg - 1] = j
                if tempkg[j] == ')':
                    if fg <= 4:
                        cp[fg - 1] = j

            # functional groups (in parenthesis)
            if fg > 0:
                for j in range(fg):
                    ib = op[j]
                    ie = cp[j]
                    if ib < ie and tempkg[ib:ie + 1] == '(O.)':
                        ialko += 1
                    elif ib < ie and tempkg[ib:ie + 1] == '(OO.)':
                        ipero += 1
                    elif ib < ie and tempkg[ib:ie + 1] == '(OH)':
                        io += 1
                        bengrp = 'O_(' + base[:point - 1]
                        nc_bg = bengrp.find(' ')
                        if nc_bg == -1:
                            nc_bg = len(bengrp)
                        if nc_bg < len(bengrp):
                            bengrp = bengrp[:nc_bg] + ')' + bengrp[nc_bg + 1:]
                        else:
                            bengrp += ')'
                        bengrp = bengrp.ljust(database.lbg)
                        getben(bengrp, bvalue, check)
                        if check != 0:
                            print('--error-- in heat. Benson group not found:', bengrp.strip())
                            print('for chemical: ', chem.strip())
                            raise SystemExit("in heat")
                        heat_val += bvalue
                    elif ib < ie and (tempkg[ib:ie + 1] == '(OOH)' or tempkg[ib:ie + 1] == '(OOOH)'):
                        io += 1
                        bengrp = 'O_(' + base[:point - 1]
                        nc_bg = bengrp.find(' ')
                        if nc_bg == -1:
                            nc_bg = len(bengrp)
                        if nc_bg + 4 < len(bengrp):
                            bengrp = bengrp[:nc_bg] + ')(O)' + bengrp[nc_bg + 4:]
                        else:
                            bengrp = bengrp[:nc_bg] + ')(O)'
                        bengrp = bengrp.ljust(database.lbg)
                        getben(bengrp, bvalue, check)
                        if check != 0:
                            print('--error-- in heat. Benson group not found:', bengrp.strip())
                            print('for chemical: ', chem.strip())
                            raise SystemExit("in heat")
                        heat_val += bvalue
                        bengrp = 'O_(O)'.ljust(database.lbg)
                        getben(bengrp, bvalue, check)
                        if check != 0:
                            print('--error-- in heat. Benson group not found:', bengrp.strip())
                            print('for chemical: ', chem.strip())
                            raise SystemExit("in heat")
                        heat_val += bvalue
                    elif ib < ie and tempkg[ib:ie + 1] == '(ONO2)':
                        init += 1
                        bengrp = 'ONO2_(' + base[:point - 1]
                        nc_bg = bengrp.find(' ')
                        if nc_bg == -1:
                            nc_bg = len(bengrp)
                        if nc_bg < len(bengrp):
                            bengrp = bengrp[:nc_bg] + ')' + bengrp[nc_bg + 1:]
                        else:
                            bengrp += ')'
                        bengrp = bengrp.ljust(database.lbg)
                        getben(bengrp, bvalue, check)
                        if check != 0:
                            print('--error-- in heat. Benson group not found:', bengrp.strip())
                            print('for chemical: ', chem.strip())
                            raise SystemExit("in heat")
                        heat_val += bvalue
                    elif ib < ie and tempkg[ib:ie + 1] == '(OONO2)':
                        io += 1
                        bengrp = 'O_(' + base[:point - 1]
                        nc_bg = bengrp.find(' ')
                        if nc_bg == -1:
                            nc_bg = len(bengrp)
                        if nc_bg + 8 < len(bengrp):
                            bengrp = bengrp[:nc_bg] + ')(ONO2)' + bengrp[nc_bg + 8:]
                        else:
                            bengrp = bengrp[:nc_bg] + ')(ONO2)'
                        bengrp = bengrp.ljust(database.lbg)
                        getben(bengrp, bvalue, check)
                        if check != 0:
                            print('--error-- in heat. Benson group not found:', bengrp.strip())
                            print('for chemical: ', chem.strip())
                            raise SystemExit("in heat")
                        heat_val += bvalue
                        bengrp = 'ONO2_(O)'.ljust(database.lbg)
                        getben(bengrp, bvalue, check)
                        if check != 0:
                            print('--error-- in heat. Benson group not found:', bengrp.strip())
                            print('for chemical: ', chem.strip())
                            raise SystemExit("in heat")
                        heat_val += bvalue
                    elif ib < ie and tempkg[ib:ie + 1] == '(NO2)':
                        io += 1
                        bengrp = 'O_(' + base[:point - 1]
                        nc_bg = bengrp.find(' ')
                        if nc_bg == -1:
                            nc_bg = len(bengrp)
                        if nc_bg + 7 < len(bengrp):
                            bengrp = bengrp[:nc_bg] + ')(NO2)' + bengrp[nc_bg + 7:]
                        else:
                            bengrp = bengrp[:nc_bg] + ')(NO2)'
                        bengrp = bengrp.ljust(database.lbg)
                        getben(bengrp, bvalue, check)
                        if check != 0:
                            print('--error-- in heat. Benson group not found:', bengrp.strip())
                            print('for chemical: ', chem.strip())
                            raise SystemExit("in heat")
                        heat_val += bvalue
                        bengrp = 'NO2_(O)'.ljust(database.lbg)
                        getben(bengrp, bvalue, check)
                        if check != 0:
                            print('--error-- in heat. Benson group not found:', bengrp.strip())
                            print('for chemical: ', chem.strip())
                            raise SystemExit("in heat")
                        heat_val += bvalue
                    else:
                        print('--error-- in heat. Benson group contribution can not')
                        print('be estimated for the group: ', tempkg.strip())
                        print('in the species: ', chem.strip())
                        raise SystemExit("in heat")

            nameben(base, point, ic, icdot, io, ico, icodot, icd, init, ialko, ipero, bengrp)
            getben(bengrp, bvalue, check)
            if check != 0:
                print('--error-- in heat. Benson group not found:', bengrp.strip())
                print('for chemical: ', chem.strip())
                raise SystemExit("in heat")
            heat_val += bvalue

    return heat_val

# =======================================================================
# NAMEBEN is used in "heat" to name benson groups
# =======================================================================
def nameben(base, point, ic, icdot, io, ico, icodot, icd, init, ialko, ipero, bengrp):
    i = 0
    bengrp = ' ' * len(bengrp)

    # form benson groups : start with the base:
    bengrp = base.ljust(len(bengrp))
    if point - 1 < len(bengrp):
        bengrp = bengrp[:point - 1] + '_' + bengrp[point:]
    point += 1

    if icodot != 0:
        for i in range(icodot):
            if point - 1 + 5 <= len(bengrp):
                bengrp = bengrp[:point - 1] + '(CO*)' + bengrp[point - 1 + 5:]
            else:
                bengrp = bengrp[:point - 1] + '(CO*)'
            point += 5

    if ic != 0:
        for i in range(ic):
            if point - 1 + 3 <= len(bengrp):
                bengrp = bengrp[:point - 1] + '(C)' + bengrp[point - 1 + 3:]
            else:
                bengrp = bengrp[:point - 1] + '(C)'
            point += 3

    if icdot != 0:
        for i in range(icdot):
            if point - 1 + 4 <= len(bengrp):
                bengrp = bengrp[:point - 1] + '(C*)' + bengrp[point - 1 + 4:]
            else:
                bengrp = bengrp[:point - 1] + '(C*)'
            point += 4

    if icd != 0:
        for i in range(icd):
            if point - 1 + 4 <= len(bengrp):
                bengrp = bengrp[:point - 1] + '(Cd)' + bengrp[point - 1 + 4:]
            else:
                bengrp = bengrp[:point - 1] + '(Cd)'
            point += 4

    if ico != 0:
        for i in range(ico):
            if point - 1 + 4 <= len(bengrp):
                bengrp = bengrp[:point - 1] + '(CO)' + bengrp[point - 1 + 4:]
            else:
                bengrp = bengrp[:point - 1] + '(CO)'
            point += 4

    if io != 0:
        for i in range(io):
            if point - 1 + 3 <= len(bengrp):
                bengrp = bengrp[:point - 1] + '(O)' + bengrp[point - 1 + 3:]
            else:
                bengrp = bengrp[:point - 1] + '(O)'
            point += 3

    if ipero != 0:
        for i in range(ipero):
            if point - 1 + 5 <= len(bengrp):
                bengrp = bengrp[:point - 1] + '(OO*)' + bengrp[point - 1 + 5:]
            else:
                bengrp = bengrp[:point - 1] + '(OO*)'
            point += 5

    if ialko != 0:
        for i in range(ialko):
            if point - 1 + 4 <= len(bengrp):
                bengrp = bengrp[:point - 1] + '(O*)' + bengrp[point - 1 + 4:]
            else:
                bengrp = bengrp[:point - 1] + '(O*)'
            point += 4

    if init != 0:
        for i in range(init):
            if point - 1 + 6 <= len(bengrp):
                bengrp = bengrp[:point - 1] + '(ONO2)' + bengrp[point - 1 + 6:]
            else:
                bengrp = bengrp[:point - 1] + '(ONO2)'
            point += 6

# =======================================================================
# Return the benson value for the group provided as input
# =======================================================================
def getben(bengrp, bvalue, check):
    i = 0
    bvalue[0] = 0.0
    check[0] = 1
    for i in range(database.nbson):
        if bengrp.strip() == database.bsongrp[i].strip():
            check[0] = 0
            bvalue[0] = database.bsonval[i]
            return