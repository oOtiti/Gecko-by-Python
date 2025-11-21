def srh5(aseek, alist, nlist):
    jhi = nlist + 1
    jlo = 1
    jold = 0
    j = 0
    srh5 = 0

    searchloop = True
    while searchloop:
        j = (jhi + jlo) // 2
        if j == jold:
            searchloop = False
            break
        
        jold = j
        if aseek > alist[j-1]:
            jlo = j
            continue
        if aseek == alist[j-1]:
            srh5 = j
            return srh5
        jhi = j
    
    srh5 = -j
    return srh5


def srch(nrec, chem, dict):
    jhi = nrec + 1
    jlo = 1
    jold = 0
    j = 0
    srch = 0

    searchloop = True
    while searchloop:
        j = (jhi + jlo) // 2
        if j == jold:
            searchloop = False
            break
        
        jold = j
        if chem > dict[j-1][9:129]:
            jlo = j
            continue
        if chem == dict[j-1][9:129]:
            srch = j
            return srch
        jhi = j
    
    srch = -j
    return srch


def search_ipos(iseek, ilist):
    search_ipos = 0
    if iseek < 1:
        print("--error-- in search_ipos, requested search using <1 int.")
        raise Exception("in search_ipos")
    
    for i in range(1, len(ilist) + 1):
        if iseek > ilist[i-1]:
            search_ipos = i
            return search_ipos
    
    if search_ipos < 1:
        print("--error-- in search_ipos, no slot found")
        raise Exception("in search_ipos")
    
    return search_ipos