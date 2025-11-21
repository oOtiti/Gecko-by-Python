def primes(nca, bond, rank):
    from keyparameter import prim

    iterloop = True
    while iterloop:
        cprim = [0] * len(rank)
        pp = [1] * len(rank)
        
        for i in range(1, nca + 1):
            cprim[i-1] = prim[rank[i-1]]
            pp[i-1] = 1
        
        for i in range(1, nca + 1):
            for j in range(1, nca + 1):
                if i != j and bond[i-1, j-1] > 0:
                    pp[i-1] = pp[i-1] * cprim[j-1]
        
        maxr = 0
        prank = [[0] * 3 for _ in range(len(rank))]
        for i in range(1, nca + 1):
            prank[i-1][0] = i
            prank[i-1][1] = rank[i-1]
            prank[i-1][2] = pp[i-1]
            maxr = max(maxr, rank[i-1])
            rank[i-1] = 0
        
        sortrank(nca, prank, rank)
        
        iterloop = False
        for i in range(1, nca + 1):
            if rank[i-1] != prank[i-1][1]:
                iterloop = True
                break


def sortrank(nca, prank, rank):
    if len(prank[0]) != 3:
        print("in sortrank, size array issue")
        raise Exception("in sortrank")
    
    pr = [row[:] for row in prank]
    
    i = 1
    while i < nca:
        ii = i + 1
        if pr[i-1][1] <= pr[ii-1][1]:
            i += 1
            continue
        for j in range(3):
            store = pr[ii-1][j]
            pr[ii-1][j] = pr[i-1][j]
            pr[i-1][j] = store
        i -= 1
        if i == 0:
            i = 1
        continue
    
    i = 1
    while i < nca:
        ii = i + 1
        if pr[i-1][1] == pr[ii-1][1]:
            if pr[i-1][2] <= pr[ii-1][2]:
                i += 1
                continue
            for j in range(3):
                store = pr[ii-1][j]
                pr[ii-1][j] = pr[i-1][j]
                pr[i-1][j] = store
            i -= 1
        else:
            i += 1
            continue
        if i == 0:
            i = 1
        continue
    
    ctr = 1
    rank[pr[0][0]-1] = ctr
    for ii in range(2, nca + 1):
        i = ii - 1
        if pr[ii-1][1] == pr[i-1][1] and pr[ii-1][2] == pr[i-1][2]:
            pass
        else:
            ctr += 1
        rank[pr[ii-1][0]-1] = ctr