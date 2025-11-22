from keyparameter import mxnode, mxlgr, mxring, mxcp
from keyflag import sar_only_fg
from atomtool import cnum, onum
from stdtool import ckgrppt, lntree, mkcopy, revers, prioty, dwrite
from stdgrbond import grbond
from stdratings import ratings
from ringtool import uniqring
from toolbox import stoperr

def stdchm(chem):
    if sar_only_fg:
        return
    
    copy = [''] * mxcp
    clngth = [[0 for _ in range(mxnode)] for _ in range(mxnode)]
    lobond = [[False for _ in range(mxnode)] for _ in range(mxnode)]
    path = [[[0 for _ in range(mxnode)] for _ in range(mxnode)] for _ in range(mxnode)]
    loze = False
    if '/' in chem or '\\' in chem:
        loze = True
        
    nca = cnum(chem) + onum(chem)
    if nca < 2:
        return
        
    nc = len(chem.strip())
    p = 0
    for i in range(nc):
        if chem[i] == '(':
            p += 1
        if chem[i] == ')':
            p -= 1
    if p != 0:
        mesg = "parentheses mismatch (error 1)"
        stoperr('stdchm', mesg, chem)
        
    bond = [[0 for _ in range(mxnode)] for _ in range(mxnode)]
    group = [''] * mxnode
    dbflg = 0
    nring = 0
    zebond = [[0 for _ in range(mxnode)] for _ in range(mxnode)]
    grbond(chem, group, bond, dbflg, nring, zebond)
    
    ncis = 0
    ntrans = 0
    if loze:
        add_zetag(chem, nca, zebond, group, ncis, ntrans)
        
    rank = [0] * mxnode
    ratings(nca, group, bond, nring, rank)
    
    if nring > 0:
        rjg = [[0, 0] for _ in range(mxring)]
        uniqring(nring, nca, group, bond, rank, rjg)
        for k in range(nring):
            i = rjg[k][0]
            j = rjg[k][1]
            bond[i][j] = 0
            bond[j][i] = 0
            
    for i in range(1, nca+1):
        for j in range(1, nca+1):
            if bond[i][j] != 0:
                lobond[i][j] = True
                
    for i in range(1, nca+1):
        locat = group[i].find(')(')
        if locat != -1:
            ckgrppt(locat, group[i])
            
    lntree(bond, 1, 2, nca, clngth, path)
    
    ncp = 0
    for i in range(1, nca+1):
        if clngth[1][i] != 0:
            leaf = path[1][i][clngth[1][i]]
            last = path[1][i][clngth[1][i]-1]
            lntree(bond, leaf, last, nca, clngth, path)
            
            for j in range(1, nca+1):
                if clngth[leaf][j] != 0:
                    ncp += 1
                    ptr = 1
                    for k in range(1, nca+1):
                        ig = path[leaf][j][k]
                        if ig != 0:
                            pg = 0
                            ng = 0
                            if k > 1:
                                pg = path[leaf][j][k-1]
                            if k < nca:
                                ng = path[leaf][j][k+1]
                            mkcopy(lobond, group, nca, rank, nring, ig, pg, ng, ptr, copy[ncp-1])
                    revers(copy[ncp-1], copy[ncp])
                    ncp += 1
                    
    if ncp > 1:
        nl = ncp
        i = 1
        while i < nl:
            j = i
            while True:
                j += 1
                if j > nl:
                    break
                if copy[i-1] == copy[j-1]:
                    copy[j-1] = ''
                    for k in range(j, nl):
                        copy[k-1] = copy[k]
                    if copy[i-1] != '':
                        ncp -= 1
                        j -= 1
            i += 1
            
    if ncp == 0:
        mesg = "no copies left "
        stoperr('stdchm', mesg, chem)
        
    if ncp == 1:
        tchem = copy[0]
    else:
        tchem = ''
        prioty(group, rank, copy, ncp, nring, tchem)
        
    if dbflg != 0:
        dwrite(tchem)
        
    nc = len(tchem.strip())
    p = 0
    for i in range(nc):
        if tchem[i] == '(':
            p += 1
        if tchem[i] == ')':
            p -= 1
    if p != 0:
        mesg = "parentheses mismatch (error 3) "
        stoperr('stdchm', mesg, chem)
        
    if loze:
        rm_zetag(tchem, ncis, ntrans)
        
    chem = tchem

def add_zetag(chem, nca, zebond, group, ncis, ntrans):
    lgr = len(group[1])
    ncis = 0
    ntrans = 0
    k=0
    zechar = ''
    for i in range(1, nca):
        for j in range(i, nca+1):
            if zebond[i][j] == 0:
                continue
                
            if zebond[i][j] == 1:
                ncis += 1
                if group[i][:2] != 'Cd' and group[j][:2] != 'Cd':
                    mesg = "cis/trans info provided for group without 'Cd'"
                    stoperr('add_zetag', mesg, chem)
                    
                if ncis == 1:
                    zechar = 'a'
                elif ncis == 2:
                    zechar = 'b'
                elif ncis == 3:
                    zechar = 'c'
                else:
                    mesg = "too many of cis C=C bond"
                    stoperr('add_zetag', mesg, chem)
                    
            elif zebond[i][j] == 2:
                ntrans += 1
                if group[i][:2] != 'Cd' and group[j][:2] != 'Cd':
                    mesg = "cis/trans info provided for group without 'Cd'"
                    stoperr('add_zetag', mesg, chem)
                    
                if ntrans == 1:
                    zechar = 'x'
                elif ntrans == 2:
                    zechar = 'y'
                elif ntrans == 3:
                    zechar = 'z'
                else:
                    mesg = "too many of trans C=C bond"
                    stoperr('add_zetag', mesg, chem)
            else:
                mesg = "cis/trans bond number is out of range (1,2)"
                stoperr('add_zetag', mesg, chem)
                
            for k in range(3, lgr+1):
                if group[i][k-1] != '1' and group[i][k-1] != '2' and group[i][k-1] != '3':
                    break
            group[i] = group[i][:k-1] + zechar + group[i][k-1:lgr]
            
            for k in range(3, lgr+1):
                if group[j][k-1] != '1' and group[j][k-1] != '2' and group[j][k-1] != '3':
                    break
            group[j] = group[j][:k-1] + zechar + group[j][k-1:lgr]

def rm_zetag(chem, ncis, ntrans):
    zechar = ''
    if ncis != 0:
        for i in range(1, ncis+1):
            if i == 1:
                zechar = 'a'
            elif i == 2:
                zechar = 'b'
            elif i == 3:
                zechar = 'c'
                
            k = chem.find(zechar)
            if k == -1:
                mesg = "expected 1st cis/trans character not identified"
                stoperr('rm_zetag', mesg, chem)
            else:
                chem = chem[:k] + '/' + chem[k+1:]
                
            k = chem.find(zechar)
            if k == -1:
                mesg = "expected 2nd cis/trans character not identified"
                stoperr('rm_zetag', mesg, chem)
            else:
                chem = chem[:k] + '\\' + chem[k+1:]
                
    if ntrans != 0:
        for i in range(1, ntrans+1):
            if i == 1:
                zechar = 'x'
            elif i == 2:
                zechar = 'y'
            elif i == 3:
                zechar = 'z'
                
            k = chem.find(zechar)
            if k == -1:
                mesg = "expected 1st cis/trans character not identified"
                stoperr('rm_zetag', mesg, chem)
            else:
                chem = chem[:k] + '/' + chem[k+1:]
                
            k = chem.find(zechar)
            if k == -1:
                mesg = "expected 2nd cis/trans character not identified"
                stoperr('rm_zetag', mesg, chem)
            else:
                chem = chem[:k] + '/' + chem[k+1:]