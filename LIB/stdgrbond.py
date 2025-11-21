import keyparameter
import toolbox

# =======================================================================
# PURPOSE: Contruction of the group matrix and the carbon-carbon bond   
# matrix of a molecule (chem) given as input. Characters (/, \) used to 
# mark cis/trans C=C bonds in chem are removed in groups. The zebond 
# matrix (optional argument!) handle the cis/trans configuration.  
#                                                                   
# First groups are constructed. A group starts with a carbon "C", an
# aromatic carbon "c", or an ether bond "-O-", and ends with either 
# a next "C" or "c", or a double-bond "=", or an ether bond "-O-", 
# or  blank "NUL". In each group the trailing parentheses are deleted 
# if it is an opening with no closing and vice-versa. Next the carbon 
# skeleton is built with carbon, oxygen (if ether) and remaining 
# parentheses and double-bonds. The bond matrix is then built based on 
# the skeleton. The program checks the valence for each carbon center.            
# =======================================================================
def grbond(chem, group, bond, dbflg, nring, zebond=None):
    skelet = [' '] * len(chem)
    ig = 0
    ng = 0
    ik = 0
    istop = 0
    nb = 0
    tp = 0
    p = 0
    nc = 0
    ring = [[0 for _ in range(len(bond[0]))] for _ in range(len(bond))]
    i = 0
    j = 0
    k = 0
    l = 0
    nok = 0
    ii = 0
    jj = 0
    kk = 0
    n = 0
    nn = 0
    dbbrk = 0
    lengr = 0
    numgr = 0
    loeter = False
    dbtype = [[0 for _ in range(len(bond[0]))] for _ in range(len(bond))]
    loze = False
    
    progname = 'grbond '
    mesg = ''

    # ----------
    # INITIALIZE
    # ----------
    ig = 0
    ik = 0
    dbflg[0] = 0
    nring[0] = 0
    dbbrk = 0
    istop = 0
    l = 0
    for idx in range(len(skelet)):
        skelet[idx] = ' '
    loeter = False
    for idx in range(len(group)):
        group[idx] = ' ' * len(group[idx])
    for row in range(len(bond)):
        for col in range(len(bond[0])):
            bond[row][col] = 0
    for row in range(len(ring)):
        for col in range(len(ring[0])):
            ring[row][col] = 0
    lengr = len(group[0])
    numgr = len(group)
    space_pos = chem.find(' ')
    if space_pos == -1:
        nc = len(chem)
    else:
        nc = space_pos

    # initialize Z/E table and check for cis/trans C=C bonds
    for row in range(len(dbtype)):
        for col in range(len(dbtype[0])):
            dbtype[row][col] = 0
    loze = False
    if zebond is not None:
        for row in range(len(zebond)):
            for col in range(len(zebond[0])):
                zebond[row][col] = 0
        if len(dbtype) * len(dbtype[0]) != len(zebond) * len(zebond[0]):
            mesg = "size of zebond table does not match size of bond"
            toolbox.stoperr(progname, mesg, chem)
    if '/' in chem:
        loze = True
    elif '\\' in chem:
        loze = True
 
    if chem[0] != 'C' and chem[0] != 'c':
        if len(chem) >= 2 and chem[0:2] != '-O':
            if len(chem) >= 4 and chem[0:4] == '=Cd1':
                dbbrk = 1
            else:
                mesg = "input to grbond must begin with C, c or -O "
                toolbox.stoperr(progname, mesg, chem)

    # -------------------------------------------------------------
    # CONSTRUCT CARBON-CENTERED GROUPS, SKELETON AND DIAGONAL TERM
    # -------------------------------------------------------------
    for i in range(nc):
        # identify group
        # --------------
        # find where the group starts
        nok = 0
        if chem[i] == 'C' or chem[i] == 'c':
            nok = 1
        if len(chem) >= i+2 and chem[i:i+2] == '-O':
            nok = 1
        if len(chem) >= i+2 and chem[i:i+2] == 'Cl':
            nok = 0
        if nok == 0:
            continue
        ig += 1

        # find where the group ends
        loopgr_break = False
        for j in range(i+1, min(i+lengr+1, nc)):
            nok = 0
            if chem[j] == 'C' or chem[j] == 'c':
                nok = 1
            if len(chem) >= j+2 and chem[j:j+2] == '-O':
                nok = 1
            if len(chem) >= j+2 and chem[j:j+2] == 'Cl':
                nok = 0
            if chem[j] == '=' or chem[j] == ' ':
                nok = 1
            if nok == 1:
                l = j - i
                istop = j - 1
                loopgr_break = True
                break
        if not loopgr_break:
            istop = min(i+lengr, nc-1)
            l = istop - i + 1
        end_idx = min(istop + 1, len(chem))
        group[ig-1] = chem[i:end_idx].ljust(len(group[ig-1]))

        # delete trailing parentheses from group
        while True:
            if l >= 1 and group[ig-1][l-1:l] == '(':
                group[ig-1] = group[ig-1][:l-1] + ' ' + group[ig-1][l:]
                l -= 1
                continue
            elif l >= 1 and group[ig-1][l-1:l] == ')':
                p = 0
                for j in range(1, l+1):
                    if group[ig-1][j-1:j] == '(':
                        p += 1
                    if group[ig-1][j-1:j] == ')':
                        p -= 1
                if p < 0:
                    group[ig-1] = group[ig-1][:l-1] + ' ' + group[ig-1][l:]
                    l -= 1
                    continue
            break

        # build skeleton
        # --------------
        ik += 1
        ik_idx = ik - 1

        # oxygen skeleton
        if group[ig-1][:2] == '-O':
            skelet[ik_idx] = 'O'
            j = 3
            # ring-joining oxygen (Oxygen is never a multiple center)
            for n in range(1, keyparameter.mxring + 1):
                if len(group[ig-1]) >= 3 and group[ig-1][2:3] == keyparameter.digit[n-1]:
                    nring[0] += 1
                    ik += 1
                    ik_idx_new = ik - 1
                    skelet[ik_idx_new] = keyparameter.digit[n-1]
                    j = 4
            if len(group[ig-1]) >= j:
                group[ig-1] = group[ig-1][:j-1] + '-' + group[ig-1][j:]
            loeter = True
        else:
            # carbon (either flavor)
            if group[ig-1][0] == 'c' or group[ig-1][0] == 'C':
                skelet[ik_idx] = group[ig-1][0]
                jj = 2
                if len(group[ig-1]) >= 2 and group[ig-1][1:2] == 'd':
                    jj = 3

                # ring-joining carbon
                rjloop_break = False
                for j in range(jj, l+1):
                    for n in range(1, keyparameter.mxring + 1):
                        if len(group[ig-1]) >= j and group[ig-1][j-1:j] == keyparameter.digit[n-1]:
                            ik += 1
                            ik_idx_new = ik - 1
                            skelet[ik_idx_new] = group[ig-1][j-1:j]
                            nring[0] += 1
                            rjloop_break = True
                            break
                    if rjloop_break:
                        rjloop_break = False
                        continue
                    break

        # put parenthesis if any
        for j in range(i + l, istop + 1):
            ik += 1
            ik_idx_new = ik - 1
            if ik_idx_new < len(skelet):
                skelet[ik_idx_new] = chem[j]
        if istop + 1 < len(chem) and chem[istop+1] == '=':
            ik += 1
            ik_idx_new = ik - 1
            if ik_idx_new < len(skelet):
                skelet[ik_idx_new] = '='

        # get the number of bonds at each node. Construct diagonal of bond matrix
        # (aromaticity is regarded as valence 3 for bookkeeping)
        # --------------------------------------------
        nb = 0
        p = 0
        for j in range(2, l+1):
            # number of functional group inside parenthesis
            if group[ig-1][j-1:j] == '(' and p == 0:
                nb += 1
            if group[ig-1][j-1:j] == '(':
                p += 1
            if group[ig-1][j-1:j] == ')':
                p -= 1

            # H-,O-atom and radical attached to a carbon:
            if p != 0:
                continue
            if group[ig-1][j-1:j] == 'H':
                nb += 1
            if j >= 2 and group[ig-1][j-2:j-1] == 'H' and group[ig-1][j-1:j] == '3':
                nb += 2
            if j >= 2 and group[ig-1][j-2:j-1] == 'H' and group[ig-1][j-1:j] == '2':
                nb += 1
            if group[ig-1][j-1:j] == 'O':
                nb += 2
            if group[ig-1][j-1:j] == '.':
                nb += 1

        # assign to the diagonal term of the bond matrix
        bond[ig-1][ig-1] = nb

    # store number of groups and rings
    ng = ig
    if nring[0] % 2 != 0:
        mesg = "ERROR in grbond: RING NOT CLOSED "
        toolbox.stoperr(progname, mesg, chem)
    nring[0] = nring[0] // 2

    # ----------------------------------------------------------
    # FIND TERMS OF BOND MATRIX (C-C, c-c, c-C, AND C-O-C BONDS)
    # ----------------------------------------------------------
    ig = 0
    for i in range(ik - 1):
        nok = 0
        if skelet[i] == 'C' or skelet[i] == 'c' or skelet[i] == 'O':
            nok = 1
        if nok == 0:
            continue
 
        ig += 1
        k = ig
        p = 0
        tp = 0
        nb = 1
        brchloop_break = False
        for j in range(i + 1, ik):
            if skelet[j] == '(':
                tp = 1
                p += 1
            elif skelet[j] == ')':
                tp = 0
                p -= 1
                if p < 0:
                    continue
            elif skelet[j] == 'C' or skelet[j] == 'c':
                k += 1
                if p == 0 and tp == 0:
                    bond[ig-1][k-1] = nb
                    bond[k-1][ig-1] = nb
                    brchloop_break = True
                    break
                elif p == 1 and tp == 1:
                    bond[ig-1][k-1] = nb
                    bond[k-1][ig-1] = nb
                    nb = 1
                    tp = 0
                    continue
                elif p != tp:
                    nb = 1
            elif skelet[j] == 'O':
                k += 1
                if p == 0 and tp == 0:
                    bond[ig-1][k-1] = nb
                    bond[k-1][ig-1] = nb
                    brchloop_break = True
                    break
                elif p == 1 and tp == 1:
                    bond[ig-1][k-1] = nb
                    bond[k-1][ig-1] = nb
                    nb = 1
                    tp = 0
                    continue
            elif j >= 1 and skelet[j-1:j+1] == '(=' and p == 1:
                nb = 2
                dbflg[0] = 1
            elif skelet[j] == '=' and p == 0:
                nb = 2
                dbflg[0] = 1
        if brchloop_break:
            continue

    # ----------------------------------------------------------
    # FIND BOND MATRIX TERMS FOR RING-JOINING CARBONS
    # ----------------------------------------------------------
    # i  = group index of current node
    # ii = group index of subsequent node 
    # j  = skeleton index of current node
    # jj = skeleton index of subsequent node
    # k  = skeleton index of current ring-join character
    # kk = skeleton index of subsequent ring-join character
    # ----------------------------------------------------------
    if nring[0] > 0:
        i = 0
        ii = 0
        j = 0
        jj = 0
        k = 0
        kk = 0

        for j in range(ik - 1):
            if skelet[j] == 'C' or skelet[j] == 'c' or skelet[j] == 'O':
                k = j
                i += 1
                nb = 0

                while True:
                    k += 1
                    for n in range(1, nring[0] + 1):
                        if k < len(skelet) and skelet[k] == keyparameter.digit[n-1]:
                            if ring[n-1][0] == 1:
                                break
                            ii = i
                            nb = 1
                            for jj in range(k + 1, ik):
                                if skelet[jj] == 'C' or skelet[jj] == 'c' or skelet[jj] == 'O':
                                    ii += 1
                                    kk = jj + 1
                                    if kk < len(skelet) and skelet[kk] == keyparameter.digit[n-1]:
                                        ring[n-1][0] = 1
                                        bond[i-1][ii-1] = nb
                                        bond[ii-1][i-1] = nb
                                        if keyparameter.digit[n-1] == '1' and dbbrk == 1:
                                            bond[i-1][ii-1] += 1
                                            bond[ii-1][i-1] += 1
                                        break
                                    else:
                                        for nn in range(1, nring[0] + 1):
                                            if kk < len(skelet) and skelet[kk] == keyparameter.digit[nn-1]:
                                                kk += 1
                                                if kk < len(skelet) and skelet[kk] == keyparameter.digit[n-1]:
                                                    ring[n-1][0] = 1
                                                    bond[i-1][ii-1] = nb
                                                    bond[ii-1][i-1] = nb
                                                    break
                                        break
                            break
                    if k >= len(skelet):
                        break

    # -----------
    # FINAL CHECK
    # -----------
    # check valence = 4 on each carbon, or '3' if aromatic carbon:
    for i in range(ng):
        nb = sum(bond[i][:ng])

        if (group[i][0] == 'C' and nb != 4) or (group[i][0] == 'c' and nb != 3):
            print('--error-- in grbond, check valence for species:')
            print(chem.strip())
            print(f'               skeleton = {("".join(skelet)).strip()}')
            print('          bond - matrix = ')
            for k in range(ng):
                print(' '.join(f'{bond[k][l]:3d}' for l in range(ng)))
            mesg = "check valence  "
            toolbox.stoperr(progname, mesg, chem)
        bond[i][i] = 0

    # check that no 'Cd' group exists without a Cd=Cd structure
    dbflg[0] = 0
    for i in range(ng):
        if 'Cd' in group[i]:
            for j in range(numgr):
                if bond[i][j] == 2:
                    dbflg[0] = 1
            if dbflg[0] != 1:
                print('--error-- in routine grbond. >C=C< expected in :')
                print(chem.strip())
                print(f'               skeleton = {("".join(skelet)).strip()}')
                print('          bond - matrix = ')
                for k in range(ng):
                    print(' '.join(f'{bond[k][l]:3d}' for l in range(ng)))
                mesg = ">C=C< expected "
                toolbox.stoperr(progname, mesg, chem)

    # check that the -O- function is not a terminal position
    if loeter:
        for i in range(ng):
            if group[i][:2] == '-O':
                nok = 0
                for j in range(ng):
                    nok += bond[i][j]
                if nok != 2:
                    print('--error-- in grbond. Wrong bond # at -O- site')
                    print('Expected valence is 2 at that site. ')
                    print(f'Species= {chem.strip()}')
                    print(f'               skeleton = {("".join(skelet)).strip()}')
                    print('          bond - matrix = ')
                    for k in range(ng):
                        print(' '.join(f'{bond[k][l]:3d}' for l in range(ng)))
                    mesg = "Wrong bond # "
                    toolbox.stoperr(progname, mesg, chem)

        # C-O-C bonds turned to 3 (instead of 1) to be recognized later. 
        for i in range(ng):
            if group[i][:2] == '-O':
                for j in range(ng):
                    if bond[i][j] == 1:
                        bond[i][j] = 3
                        bond[j][i] = 3

    # remove the cis/trans character (/ and \) in groups and get zebond table
    if loze:
        if dbflg[0] == 0:
            mesg = "cis-trans (/ or \) found without C=C structure"
            toolbox.stoperr(progname, mesg, chem)
        get_zebond(chem, bond, ng, group, dbtype)
        if zebond is not None:
            for i in range(len(dbtype)):
                for j in range(len(dbtype[0])):
                    zebond[i][j] = dbtype[i][j]

# =======================================================================
# Purpose: 
# 1. Remove the characters (/ or \) in groups used to indicate the cis
#    or trans configuration on a C=C bond.
# 2. Create a table of cis (Z) or trans (E) C=C bond type. The element
#    zebond(i,j) is set to 1 if Cdi and Cdj are connected via a cis 
#    configuration and set 2 if Cdi and Cdj are connected via a trans 
#    configuration.
# The cis/trans configuration is kept only for RCH=CHR structures here.    
# =======================================================================
def get_zebond(chem, bond, ng, group, zebond):
    i = 0
    j = 0
    lengr = 0
    iaslash = 0
    islash = 0
    jaslash = 0
    jslash = 0
    progname = 'get_zebond '
    mesg = ''

    for row in range(len(zebond)):
        for col in range(len(zebond[0])):
            zebond[row][col] = 0
    lengr = len(group[0])
  
    # examine the various cases (//, \\, /\, \/) of Cd=Cd bonds
    for i in range(ng - 1):
        if group[i][:2] == 'Cd':
            islash = group[i].find('/')
            iaslash = group[i].find('\\')       
 
            # start with >Cdi/= 
            if islash != -1:
                group[i] = group[i][:islash] + group[i][islash+1:]
                group[i] = group[i].ljust(lengr)
                for j in range(i + 1, ng):
                    if bond[i][j] == 2:
                        jslash = group[j].find('/')
                        jaslash = group[j].find('\\')
                        if (jslash == -1) and (jaslash == -1):
                            mesg = "missing cis-trans character on C=C structure"
                            toolbox.stoperr(progname, mesg, chem)

                        #>Cdi/=Cdj/< (trans bond)
                        if jslash != -1:
                            group[j] = group[j][:jslash] + group[j][jslash+1:]
                            group[j] = group[j].ljust(lengr)
                            zebond[i][j] = 2
                            zebond[j][i] = 2
                        #>Cdi/=Cdj\< (cis bond)
                        elif jaslash != -1:
                            group[j] = group[j][:jaslash] + group[j][jaslash+1:]
                            group[j] = group[j].ljust(lengr)
                            zebond[i][j] = 1
                            zebond[j][i] = 1
                        else:
                            mesg = "cis/trans structure cannot be identified"
                            toolbox.stoperr(progname, mesg, chem)
                        break
 
            # start with >Cdi\=
            elif iaslash != -1:
                group[i] = group[i][:iaslash] + group[i][iaslash+1:]
                group[i] = group[i].ljust(lengr)
                for j in range(i + 1, ng):
                    if bond[i][j] == 2:
                        jslash = group[j].find('/')
                        jaslash = group[j].find('\\')
                        if (jslash == -1) and (jaslash == -1):
                            mesg = "missing cis-trans character on C=C structure"
                            toolbox.stoperr(progname, mesg, chem)
                        
                        #>Cdi\=Cdj/< (cis case)
                        if jslash != -1:
                            group[j] = group[j][:jslash] + group[j][jslash+1:]
                            group[j] = group[j].ljust(lengr)
                            zebond[i][j] = 1
                            zebond[j][i] = 1
                        #>Cdi\=Cdj\< (trans case)
                        elif jaslash != -1:
                            group[j] = group[j][:jaslash] + group[j][jaslash+1:]
                            group[j] = group[j].ljust(lengr)
                            zebond[i][j] = 2
                            zebond[j][i] = 2
                        else:
                            mesg = "cis/trans structure cannot be identified"
                            toolbox.stoperr(progname, mesg, chem)
                        break

    # cis/trans only considered for RCH=CH-R. Kill other cases
    for i in range(ng - 1):
        for j in range(ng):
            if zebond[i][j] != 0:
                if group[i].strip() == 'CdH2':
                    zebond[i][j] = 0
                    zebond[j][i] = 0
                if group[j].strip() == 'CdH2':
                    zebond[i][j] = 0
                    zebond[j][i] = 0
                if 'H' not in group[i]:
                    zebond[i][j] = 0
                    zebond[j][i] = 0
                if 'H' not in group[j]:
                    zebond[i][j] = 0
                    zebond[j][i] = 0