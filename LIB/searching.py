import sys
#=======================================================================
# PURPOSE: Binary tree search for string      
#
# The value returned (srh5) is: 
#    <0 if not found (|srh5|==pointer where aseek should be inserted) 
#
# >0 if found (|srh5|==pointer where aseek is inserted in alist)   
#=======================================================================
def srh5(aseek, alist, nlist):
   
    srh5_val = 0
    jold = 0
    jlo = 1
    # jhi  = nlist + 1
    jhi = nlist + 1

    # search loop
    while True:
        j = (jhi + jlo) // 2
        if j == jold:
            break 

        jold = j
        
        # Note: Fortran is 1-based, Python is 0-based. 
        # Adjusting access to alist by -1 to maintain logic of j.
        if aseek > alist[j-1]:
            jlo = j
            continue # CYCLE searchloop
        
        if aseek == alist[j-1]: 
            srh5_val = j
            return srh5_val
            
        jhi = j

    # string not found 
    srh5_val = -j
    return srh5_val
#=======================================================================
#
# PURPOSE: Binary tree search for formula in dictionary      
# The value returned (srch) is: 
#
# <0 if not found (|srch|==pointer where chem should be inserted) 
#
# >0 if found (|srch|==pointer where chem is inserted in dict)
#
# WARNING (BA) : "same" routine as srh5 => merge later    
#=======================================================================
def srch(nrec, chem, dict_arr):
    
    srch_val = 0
    jold = 0
    jlo = 1
    # jhi  = nrec + 1
    jhi = nrec + 1

    # search loop
    while True:
        j = (jhi + jlo) // 2
        if j == jold:
            break 

        jold = j
        
        # Note: Fortran substring (10:129) is 1-based inclusive.
        # Python slice [9:129] starts at index 9 (10th char) up to 129 (exclusive).
        # dict access is adjusted by [j-1].
        if chem > dict_arr[j-1][9:129]:
            jlo = j
            continue # CYCLE searchloop
        
        if chem == dict_arr[j-1][9:129]: 
            srch_val = j
            return srch_val
            
        jhi = j

    # string not found 
    srch_val = -j
    return srch_val
#=======================================================================
# PURPOSE: return the rank (position) of an integer provided as input 
# (iseek) in a sorted list of integer (from large to small number). The
# list provided as input (ilist) is expected to be short: a simple 
# search is performed.
# Turn the function using binary tree search for 
# long list.
# !=======================================================================
def search_ipos(iseek, ilist):

    search_ipos_val = 0
    if iseek < 1:
        print("--error-- in search_ipos, requested search using <1 int.")
        # STOP "in search_ipos"
        sys.exit("in search_ipos")
    
    for i in range(1, len(ilist) + 1):
        if iseek > ilist[i-1]:
            search_ipos_val = i
            return search_ipos_val
 
    if search_ipos_val < 1:
        print("--error-- in search_ipos, no slot found")
        # STOP "in search_ipos"
        sys.exit("in search_ipos")
    
    return search_ipos_val
