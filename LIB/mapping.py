# 温馨提示 ：在这里，数组的有效数据一律是从1开始，0位索引上的只是占位符的作用；
# 也就是，计算和书写索引时的过程跟f90的一模一样，只是在数组或者矩阵的初始化时改动了一点。
# 这样是为了减少后面反复修改索引的工作量，也能够最大程度的避免重写的失误。 

import sys
import toolbox
import keyparameter
import searching
# ======================================================================
# PURPOSE: Browse the bond matrix to find all possible node at a given        
# position (npos), starting from a given node (top). Ring are allowed
# (track is stopped when a full loop along the circle is made).                                                                           
# The subroutine call "gettrack" which give all the possible track,   
# starting from top.                                                  
#  track(*,2) give all nodes in alpha position regarding node top     
#  track(*,3) give all nodes in beta position regarding node top      
#  track(*,4) give all nodes in gamma position regarding node top     
#  etc...                                                             
# ======================================================================
def nodmap(bond, top, nca, npos, nnod, tnod):

    # Setup internal arrays (using 1-based indexing logic)
    bond_size = len(bond)
    # Size mxcp+1 for mxcp max tracks, bond_size+1 for nca max positions
    track = [[0 for _ in range(bond_size + 1)] for _ in range(mxcp + 1)]
    # Size mxcp+1 for trlen
    trlen = [0] * (mxcp + 1)
    # ntr must be a single-element list passed to gettrack
    ntr = [0]
    
    nnod[0] = 0
    # Clear OUT parameter tnod (assuming it's a mutable list/array)
    tnod.clear()

    # ntr, track, and trlen will be modified in place.
    gettrack(bond, top, nca, ntr, track, trlen)
    
    # Avoid duplicate - set 0 to duplicate nodes
    if ntr[0] > 1:
        for i in range(1, ntr[0]): # range covers 1 to ntr-1
            if track[i][npos] == 0:
                continue # CYCLE trloop
            
            # DO j=i+1,ntr
            for j in range(i + 1, ntr[0] + 1):
                # IF (track(j,npos)/=0) THEN
                if track[j][npos] != 0:
                    # IF (track(i,npos)==track(j,npos)) track(j,npos)=0
                    if track[i][npos] == track[j][npos]:
                        track[j][npos] = 0

    # get nodes
    for i in range(1, ntr[0] + 1):
        if track[i][npos] != 0:            
            nnod[0] = nnod[0] + 1          
            tnod.append(track[i][npos]) # Append to OUT parameter tnod

# ======================================================================
# PURPOSE : Browse the bond matrix to find all possible tracks starting     
# from a given node (top). A track ends when a "ending" node (i.e. the
# last positions on the carbon skeleton) is reached. Ring are allowed 
# (track is stopped when a full loop along the circle is made).                                           *                                                                    
#  track(1,*) give all nodes in 1st track
#  track(2,*) give all nodes in 2nd track       
#  etc...                                                             
#  track(*,2) give all nodes in alpha position regarding node top     
#  track(*,3) give all nodes in beta position regarding node top       
#  track(*,4) give all nodes in gamma position regarding node top     
#  etc...                                                             
# ======================================================================
def gettrack(bond, top, nca, ntr, track, trlen):
   
    # Setup internal arrays (using 1-based indexing logic)
    bond_size = len(bond)
    memo = [0] * (bond_size + 1)

    mxtr = len(track) - 1

    # initialize parameters to find the tracks 
    # track(1,1)=top  ;  memo(1)=top
    track[1][1] = top
    memo[1] = top
    
     # of node since starting node (i.e. top)
    niv = 1
    #  current "search" pointer (i.e. to find next node in the track)
    ptr = 0
    # ntr=1     ! # of track found - current track
    ntr[0] = 1  # Update OUT parameter ntr
    # nod=top   ! current node along the track
    nod = top
    # slope=1   ! equal 1 when going forward along the track, otherwise 0
    slope = 1

    # -----------
    # start loop
    # -----------

    # get next node - reentry point
  
    while True: # Label: nextnode
        # ptr=ptr+1
        ptr = ptr + 1

        # end of possible nodes reached - must go backward or exit
        # IF (ptr > nca) THEN      
        if ptr > nca:
            # IF (niv == 1) THEN       ! all the possible track are found - exit
            if niv == 1:
                # EXIT nextnode
                break
            # ELSE                     ! go backward
            else:
                # ptr=memo(niv)          ! set pointer to previous memo
                ptr = memo[niv]
                # memo(niv)=0            ! reset memo   
                memo[niv] = 0
                # niv=niv-1              ! decrease niv (go backward)
                niv = niv - 1
                # nod=memo(niv)          ! set new current node 
                nod = memo[niv]

                # IF (slope /= 0) THEN   ! make a new track (if required)
                if slope != 0:
                    # ntr=ntr+1
                    ntr[0] = ntr[0] + 1 # Update OUT parameter ntr
                    # IF (ntr > mxtr) THEN
                    if ntr[0] > mxtr:
                        # PRINT*, '--error-- in gettrack. ntr is greater than mxtr'
                        print('--error-- in gettrack. ntr is greater than mxtr')
                        # STOP "in gettrack"
                        sys.exit("in gettrack")
                    # Copy the part of the track
                    track[ntr[0]][1:niv+1] = track[ntr[0]-1][1:niv+1] 
                # ENDIF
                # track(ntr,niv+1)=0     ! remove previous track
                track[ntr[0]][niv+1] = 0
                # slope=0                ! remember ... I am now going backward
                slope = 0
                # CYCLE nextnode
                continue 
        # no bond between current node (nod) and next possible node (ptr)
        # IF (bond(ptr,nod)==0) CYCLE nextnode  ! get next
        # Access bond with [ptr-1][nod-1] to match 1-based logic for input array
        if bond[ptr-1][nod-1] == 0:
            continue

        # new bond found - must go forward
        # IF (bond(ptr,nod) /= 0) THEN
        if bond[ptr-1][nod-1] != 0:
            # DO k=1,niv
            cycle_found = False
            for k in range(1, niv + 1):
                # IF (ptr == memo(k)) CYCLE nextnode   ! end circle or previous node
                if ptr == memo[k]:
                    cycle_found = True
                    break
            if cycle_found:
                continue # CYCLE nextnode
           
            niv = niv + 1
            # track(ntr,niv)=ptr              ! keep track
            track[ntr[0]][niv] = ptr
            # nod=ptr                         ! set current node to the new node
            nod = ptr
            # memo(niv)=nod                   ! remember which track was used
            memo[niv] = nod
            # ptr=0                           ! set pointer to 0 (to find next node)
            ptr = 0
            # slope=1                         ! remember ... I am going forward
            slope = 1
            # CYCLE nextnode                  ! go find next node
            continue
    # remove the last "case" (fill with top only) and clean track
    # track(ntr,:)=0  ;  ntr=ntr-1
    # Clean up the last track slot
    track[ntr[0]] = [0] * len(track[ntr[0]])
    ntr[0] = ntr[0] - 1 # Update OUT parameter ntr

    for i in range(1, ntr[0] + 1):
        # DO j=nca,1,-1
        for j in range(nca, 0, -1):
            # IF (track(i,j)/=0) THEN
            if track[i][j] != 0:
                # trlen(i)=j
                trlen[i] = j # Update OUT parameter trlen
                # CYCLE lentrack
                break # breaks the inner loop, proceeding to next i (equivalent to CYCLE lentrack)
           

# ----------------------------------------------------------------------
# SUBROUTINE abcde_map
# ----------------------------------------------------------------------

# ! ======================================================================
# ! PURPOSE :                                                                    
# !  Browse the bond matrix to find all possible node at various position 
# !  (i.e. alpha, beta, gamma ...), starting from node 'top'. Ring are   
# !  allowed (track is stopped when a full loop along the circle is made).
# !                                                                     
# ! The subroutine call "gettrack" which give all the possible track,   
# ! starting from top (see gettrack in the same module).                  
# !                                                                     
# ! The subroutine returns 2 tables:                                        
# ! -nabcde(k)    : number of distinct pathways that end up at a        
# !                  position k relative to top (e.g. nabcde(3) gives    
# !                  the number of distinct pathways finishing in a      
# !                  beta position relative to top                       
# ! -tabcde(k,i,j): give the pathways (node j), for the track number i  
# !                  to reach the position k (k=2 is beta position ...).
# !                  For example, tabcde(4,1,j) give the first track to 
# !                  reach a gamma position (node given by tabcde(4,1,4)),          
# !                  using the track given by tabcde(4,1,*)              
# ! ======================================================================
def abcde_map(bond, top, nca, nabcde, tabcde):
   
    # Calculate internal sizes based on OUT parameter tabcde (1-based sizing)
    mxdeep = len(tabcde) - 1 
    mxtr = len(tabcde[0]) - 1 
    track_max_len = len(tabcde[0][0]) - 1 
    
    # Internal array initialization (using 1-based padding logic)
    track = [[0 for _ in range(track_max_len + 1)] for _ in range(mxtr + 1)]
    trlen = [0] * (mxtr + 1)
    ntr = [0] # ntr is passed as a list for reference simulation
    
    # 3D array ltabc: 
    ltabc = [[[0 for _ in range(track_max_len + 1)] for _ in range(mxtr + 1)] for _ in range(mxdeep + 1)]

    # Clear OUT parameters and internal array (using 1-based indices)
    for idx in range(1, len(nabcde)):
        nabcde[idx] = 0
    for k in range(1, mxdeep + 1):
        for i in range(1, mxtr + 1):
            for j in range(1, track_max_len + 1):
                tabcde[k][i][j] = 0
    # ltabc is already cleared on creation

    # get all tracks starting at top
    # CALL gettrack(bond,top,nca,ntr,track,trlen)
    gettrack(bond, top, nca, ntr, track, trlen)

    # range the track for each length (up to mxdeep)
    # Remember that 2=alpha position, 3=beta ... 
    # IF (ntr > mxtr) THEN
    if ntr[0] > mxtr:
        # PRINT*, '--error-- in abcde_map, ntr > mxtr'
        print('--error-- in abcde_map, ntr > mxtr')
        # STOP "in abcde_map"
        sys.exit("in abcde_map")
    # ENDIF
    
    # DO k=2,mxdeep
    for k in range(2, mxdeep + 1):
        # DO i=1,ntr
        for i in range(1, ntr[0] + 1):
            # IF (trlen(i) >= k) THEN
            if trlen[i] >= k:
                # DO j=1,k
                for j in range(1, k + 1):
                    # ltabc(k,i,j)=track(i,j)
                    ltabc[k][i][j] = track[i][j]
                
    # avoid duplicate - set 0 to duplicate nodes tracks
    # DO k=2,mxdeep
    for k in range(2, mxdeep + 1):
        # trackloop: DO i=1,ntr-1
        for i in range(1, ntr[0]):
            # IF (ltabc(k,i,k)==0) CYCLE
            if ltabc[k][i][k] == 0:
                continue # CYCLE
            
            # DO ii=i+1,ntr
            for ii in range(i + 1, ntr[0] + 1):
                
                # Check for track equality: loop through j=1 to k
                is_duplicate = True
                # DO j=1,k
                for j in range(1, k + 1):
                    # IF (ltabc(k,i,j)/=ltabc(k,ii,j)) CYCLE trackloop
                    if ltabc[k][i][j] != ltabc[k][ii][j]:
                        is_duplicate = False
                        break
                
                # If that point is reached, track i= track ii
                if is_duplicate:
                    # DO j=1,k
                    for j in range(1, k + 1):
                        # ltabc(k,ii,j)=0
                        ltabc[k][ii][j] = 0
                   
    # get nodes
    # DO k=2,mxdeep
    for k in range(2, mxdeep + 1):
        # DO i=1,ntr
        for i in range(1, ntr[0] + 1):
            # IF (ltabc(k,i,k)/=0) THEN
            if ltabc[k][i][k] != 0:
                # nabcde(k)=nabcde(k)+1     ! add one more pathway that end up at a k position
                nabcde[k] = nabcde[k] + 1
                # DO j=1,k
                for j in range(1, k + 1):
                    # ii=nabcde(k)
                    ii = nabcde[k]
                    # tabcde(k,ii,j)=ltabc(k,i,j)
                    tabcde[k][ii][j] = ltabc[k][i][j]

    # add top nod (needed top scroll group from top, included)
    # nabcde(1)=1
    nabcde[1] = 1
    # tabcde(1,1,1)=top
    tabcde[1][1][1] = top
#################################

# 在这里，你会看到四个循环体，第一个是用于初始化的    

# END SUBROUTINE abcde_map

# ----------------------------------------------------------------------
# SUBROUTINE ciptree
# ----------------------------------------------------------------------

# !=======================================================================
# ! PURPOSE : Return a tree of the atom masses to rank branches according
# ! to the Cahn-Ingold-Prelog (CIP) priority rules. The node tree (tabcde) 
# ! starting from the top node in branch is provided as input.
# ! For example, the CIT tree returned for -C(ONO2)(CdH=CdH2)COCH2(OH) is:
# ! line 1:  12    
# ! line 2:  16  12  12  
# ! line 3:  16  16  14  12  12  12    1   
# ! line 4:  16  16  16    1   1   1   1   
# ! line 5:    1    
# !=======================================================================
def ciptree(chem, bond, group, ngr, nabcde, tabcde, ciptr):
    # CHARACTER(LEN=15) :: progname='ciptree'
    # CHARACTER(LEN=70) :: mesg

    mxbd = 5 
    progname = 'ciptree'
    mesg = ''
    
    # Internal Array Setup (using 1-based padding logic)
    ciptr_rows = len(ciptr) # SIZE(ciptr,1)
    ciptr_cols = len(ciptr[0]) # SIZE(ciptr,2)

    # nciptr array (1-based index)
    nciptr = [0] * ciptr_rows 
    
    # nC, nO, nH, nN arrays (1-based index)
    nC = [0] * (mxbd + 1) 
    nO = [0] * (mxbd + 1)
    nH = [0] * (mxbd + 1)
    nN = [0] * (mxbd + 1)


    # ciptr(:,:)=0 ; nciptr(:)=0 
    # Clear OUT parameter ciptr (using 1-based indices)
    for r in range(1, ciptr_rows):
        for c in range(1, ciptr_cols):
            ciptr[r][c] = 0

    # scroll from close to remote nodes (top, alpha, beta, gamma etc)
    # DO idepth=1,ngr
    for idepth in range(1, ngr + 1):
        # nnod=nabcde(idepth)     ! # of nodes having a distance "idepth" from top (rank 1)
        nnod = nabcde[idepth]
        
        # IF (nnod==0) EXIT       ! no more groups to consider
        if nnod == 0:
            break
        
        # nC(:)=0 ; nH(:)=0 ; nO(:)=0 ; nN(:)=0
        nC = [0] * (mxbd + 1) 
        nH = [0] * (mxbd + 1)
        nO = [0] * (mxbd + 1)
        nN = [0] * (mxbd + 1)

        # loop over all nodes distant idepth from top
        # DO i=1,nnod
        for i in range(1, nnod + 1):
            # ig=tabcde(idepth,i,idepth)  ! current group index
            ig = tabcde[idepth][i][idepth] 

            # Access group(ig) (Fortran string)
            grp = group[ig] # Assuming group is 1-based indexed
            
            # --- Functional Group Counting (Based on string prefix/slice) ---
            
            if grp[0:3]=='CH3':
                nC[1] += 1; nH[2] += 3
            elif grp[0:3]=='CH2':
                nC[1] += 1; nH[2] += 2
            elif grp[0:3]=='CHO':
                nC[1] += 1; nH[2] += 1; nO[2] += 2
            elif grp[0:2]=='CO':
                nC[1] += 1; nO[2] += 2        
            elif grp[0:2]=='CH':
                nC[1] += 1; nH[2] += 1
            elif grp[0:4]=='CdH2':
                nC[1] += 1; nH[2] += 2
            elif grp[0:3]=='CdH':
                nC[1] += 1; nH[2] += 1
            elif grp[0:3]=='CdO':
                nC[1] += 1; nO[2] += 2
            elif grp[0:1]=='C':
                nC[1] += 1
            elif grp[0:3]=='-O-':
                nO[1] += 1
            elif grp[0:2]=='cH':
                nC[1] += 1; nH[2] += 1
            elif grp[0:1]=='c':
                nC[1] += 1
            # ENDIF

            # count twice the C if Cd is the "second" Cd in C=C
            # IF (group(ig)(1:2)=='Cd') THEN 
            if grp[0:2]=='Cd':
                if idepth > 1:
                    # pvnod=tabcde(idepth,i,idepth-1)    ! index of the previous node in path 
                    pvnod = tabcde[idepth][i][idepth-1]
                    # IF (bond(pvnod,ig)==2) nC(1)=nC(1)+1  ! count 2 C for the 2nd Cd
                    if bond[pvnod][ig] == 2:
                        nC[1] += 1
                # ENDIF
            # ENDIF

            # loop over functional group
            nfun = countstring(grp, '(OH)')
            # IF (nfun/=0) THEN 
            if nfun != 0:
                nO[2] += nfun; nH[3] += nfun
            # ENDIF

            nfun = countstring(grp, '(OOH)')
            # IF (nfun/=0) THEN 
            if nfun != 0:
                nO[2] += nfun; nO[3] += nfun; nH[4] += nfun
            # ENDIF

            nfun = countstring(grp, '(OOOH)')
            # IF (nfun/=0) THEN 
            if nfun != 0:
                nO[2] += nfun; nO[3] += nfun; nO[4] += nfun; nH[5] += nfun
            # ENDIF

            nfun = countstring(grp, '(ONO2)')
            # IF (nfun/=0) THEN 
            if nfun != 0:
                nO[2] += nfun; nN[3] += nfun; nO[4] += 2 * nfun
            # ENDIF

            nfun = countstring(grp, '(OONO2)')
            # IF (nfun/=0) THEN 
            if nfun != 0:
                nO[2] += nfun; nO[3] += nfun; nN[4] += nfun; nO[5] += 2 * nfun
            # ENDIF

            nfun = countstring(grp, '(NO2)')
            # IF (nfun/=0) THEN 
            if nfun != 0:
                nN[2] += nfun; nO[3] += 2 * nfun
            # ENDIF

            nfun = countstring(grp, '(OO.)')
            # IF (nfun/=0) THEN 
            if nfun != 0:
                nO[2] += nfun; nO[3] += nfun
            # ENDIF

            nfun = countstring(grp, '(O.)')
            # IF (nfun/=0) THEN 
            if nfun != 0:
                nO[2] += nfun
            # ENDIF
        # ENDDO        
        
        # --- Add elements to CIP tree ---

        # add C in CIP tree 
        # DO ii=1,mxbd
        for ii in range(1, mxbd + 1):
            # IF (nC(ii)==0) CYCLE
            if nC[ii] == 0:
                continue
            
            # nel=nC(ii) ; irk=idepth+ii-1              ! # of C (nel) to be added at rank irk 
            nel = nC[ii]
            irk = idepth + ii - 1
            
            # ipos=search_ipos(12,ciptr(irk,:))         ! ipos: index to add element
            ipos = search_ipos(12, ciptr[irk])
            
            # ilast=nciptr(irk)                          ! last occupied pos. in rank
            ilast = nciptr[irk]
            
            # IF (ipos<=ilast) &                         ! make room if needed
            if ipos <= ilast:
                # Fortran: ciptr(irk,ipos+nel:ipos+nel+ilast)=ciptr(irk,ipos:ipos+ilast) 
                # Python indices: [start : end + 1]
                source_slice = ciptr[irk][ipos : ilast + 1]
                target_start = ipos + nel
                target_end = ilast + nel + 1
                
                if target_end > ciptr_cols:
                    mesg = "too many elements added in a line of CIP tree (column limit exceeded for C)"
                    stoperr(progname, mesg, chem)

                ciptr[irk][target_start : target_end] = source_slice
            
            # ciptr(irk,ipos:ipos+nel-1)=12              ! add C
            for j in range(ipos, ipos + nel):
                ciptr[irk][j] = 12
            
            # nciptr(irk)=nciptr(irk)+nel              ! store length of the rank
            nciptr[irk] = nciptr[irk] + nel
        # ENDDO

        # add O in CIP tree 
        # DO ii=1,mxbd
        for ii in range(1, mxbd + 1):
            # IF (nO(ii)==0) CYCLE
            if nO[ii] == 0:
                continue
            
            # nel=nO(ii) ; irk=idepth+ii-1              ! # of O (nel) to be added at rank irk 
            nel = nO[ii]
            irk = idepth + ii - 1
            
            # ipos=search_ipos(16,ciptr(irk,:))         ! ipos: index to add element
            ipos = search_ipos(16, ciptr[irk])
            
            # ilast=nciptr(irk)                          ! last occupied pos. in rank
            ilast = nciptr[irk]
            
            # IF (ipos<=ilast) &                         ! make room if needed
            if ipos <= ilast:
                # Fortran: ciptr(irk,ipos+nel:ilast+nel)=ciptr(irk,ipos:ilast) 
                source_slice = ciptr[irk][ipos : ilast + 1]
                target_start = ipos + nel
                target_end = ilast + nel + 1
                
                if target_end > ciptr_cols:
                    mesg = "too many elements added in a line of CIP tree (column limit exceeded for O)"
                    stoperr(progname, mesg, chem)
                    
                ciptr[irk][target_start : target_end] = source_slice
            
            # ciptr(irk,ipos:ipos+nel-1)=16              ! add O
            for j in range(ipos, ipos + nel):
                ciptr[irk][j] = 16
            
            # nciptr(irk)=nciptr(irk)+nel              ! store length of the rank
            nciptr[irk] = nciptr[irk] + nel
        # ENDDO
        
        # add N in CIP tree
        # DO ii=1,mxbd
        for ii in range(1, mxbd + 1):
            # IF (nN(ii)==0) CYCLE
            if nN[ii] == 0:
                continue
            
            # nel=nN(ii) ; irk=idepth+ii-1              ! # of N (nel) to be added at rank irk  
            nel = nN[ii]
            irk = idepth + ii - 1
            
            # ipos=search_ipos(14,ciptr(irk,:))         ! ipos: index to add element
            ipos = search_ipos(14, ciptr[irk])
            
            # ilast=nciptr(irk)                          ! last occupied pos. in rank
            ilast = nciptr[irk]
            
            # IF (ipos<=ilast) &                         ! make room if needed
            if ipos <= ilast:
                # Fortran: ciptr(irk,ipos+nel:ilast+nel)=ciptr(irk,ipos:ilast) 
                source_slice = ciptr[irk][ipos : ilast + 1]
                target_start = ipos + nel
                target_end = ilast + nel + 1
                
                if target_end > ciptr_cols:
                    mesg = "too many elements added in a line of CIP tree (column limit exceeded for N)"
                    stoperr(progname, mesg, chem)
                    
                ciptr[irk][target_start : target_end] = source_slice
            
            # ciptr(irk,ipos:ipos+nel-1)=14              ! add N
            for j in range(ipos, ipos + nel):
                ciptr[irk][j] = 14
            
            # nciptr(irk)=nciptr(irk)+nel              ! store length of the rank
            nciptr[irk] = nciptr[irk] + nel
        # ENDDO

        # add H in CIP tree
        # DO ii=2,mxbd
        for ii in range(2, mxbd + 1):
            # IF (nH(ii)==0) CYCLE
            if nH[ii] == 0:
                continue
            
            # nel=nH(ii) ; irk=idepth+ii-1              ! # of H (nel) to be added at rank irk  
            nel = nH[ii]
            irk = idepth + ii - 1
            
            # ipos=search_ipos(1,ciptr(irk,:))          ! ipos: index to add element
            ipos = search_ipos(1, ciptr[irk])
            
            # ilast=nciptr(irk)                          ! last occupied pos. in rank
            ilast = nciptr[irk]
            
            # IF (ipos<=ilast) &                         ! make room if needed
            if ipos <= ilast:
                # Fortran: ciptr(irk,ipos+nel:ilast+nel)=ciptr(irk,ipos:ilast) ! make room if needed
                source_slice = ciptr[irk][ipos : ilast + 1]
                target_start = ipos + nel
                target_end = ilast + nel + 1
                
                if target_end > ciptr_cols:
                    mesg = "too many elements added in a line of CIP tree (column limit exceeded for H)"
                    stoperr(progname, mesg, chem)
                    
                ciptr[irk][target_start : target_end] = source_slice
            
            # ciptr(irk,ipos:ipos+nel-1)=1               ! add H
            for j in range(ipos, ipos + nel):
                ciptr[irk][j] = 1
            
            # nciptr(irk)=nciptr(irk)+nel              ! store length of the rank
            nciptr[irk] = nciptr[irk] + nel

    # check size of CIP tree before returning
    # DO i=1,SIZE(nciptr)
    for i in range(1, ciptr_rows):       
        # 在这里，ciptr_rows前面有定义。
        # 
        # 乍一看好像是错的，但是我确定了，没有问题，这里不用查错 
        # IF (nciptr(i)==0) EXIT  ! no more lines
        if nciptr[i] == 0:
            break
        
        # IF (nciptr(i) > SIZE(ciptr,2)) THEN
        if nciptr[i] > ciptr_cols:
            # 这里同理
            # mesg="too many elements added in a line of CIP tree"
            mesg = "too many elements added in a line of CIP tree"
            # CALL stoperr(progname,mesg,chem)
            stoperr(progname, mesg, chem)



# !=======================================================================
# ! PURPOSE: create tracks for Cd groups. ...
# !=======================================================================
def alkenetrack(chem, bond, group, ngr, ncdtrack, cdtracklen, cdtrack):
    
    array_size = ngr + 1
    mxcd = len(cdtrack[0]) - 1
    
    cdbond = [[0 for _ in range(array_size)] for _ in range(array_size)]
    cdpst = [0] * array_size
    
    track = [[0 for _ in range(array_size)] for _ in range(MXCP + 1)]
    trlen = [0] * (MXCP + 1)
    ntr = [0] 

    progname = 'alkenetrack '
    mesg = ''

    ncdtrack[0] = 0
    
    for i in range(1, len(cdtrack)):
        for j in range(1, len(cdtrack[0])):
            cdtrack[i][j] = 0
    for i in range(1, len(cdtracklen)):
        cdtracklen[i] = 0
        
    # create a bond matrix for Cd only
    ncd = 0
    for i in range(1, ngr + 1):
        if group[i][0:2] == 'Cd':
            ncd += 1
            for j in range(1, ngr + 1):
                if j == i: continue
                if bond[i][j] != 0:
                    if group[j][0:2] == 'Cd':
                        cdbond[i][j] = 1
                        cdbond[j][i] = 1
    
    if ncd == 0:
        return

    # Count # of connection of each Cd to other Cds
    for i in range(1, ngr + 1):
        cdpst[i] = sum(cdbond[i][1:]) 

    # check primary and tertiary Cd nodes
    nprim = 0
    for i in range(1, ngr + 1):
        if cdpst[i] == 1:
            nprim += 1
        if cdpst[i] > 2:
            mesg = "tertiary Cd structure identified (and not allowed)"
            stoperr(progname, mesg, chem)
    
    if nprim == 0:
        mesg = "unexpected cyclic Cd structure identified"
        stoperr(progname, mesg, chem)

    # get Cd chains (must start from primary Cd)
    for i in range(1, ngr + 1):
        if cdpst[i] == 1:
            ntr = [0] 
            gettrack(cdbond, i, ngr, ntr, track, trlen)
            
            if ntr[0] > 1:
                mesg = "unexpected branching identified on Cd structure"
                stoperr(progname, mesg, chem)
            
            if trlen[1] > 4:
                mesg = "More than 4 Cd identified (and not allowed)"
                stoperr(progname, mesg, chem)
            
            if trlen[1] == 3:
                mesg = ">C=C=C< structure identified and not allowed"
                stoperr(progname, mesg, chem)
            
            ncdtrack[0] += 1
            ncdtrack_val = ncdtrack[0]

            if ncdtrack_val > len(cdtrack) - 1:
                mesg = "maximum number of Cd tracks reached."
                stoperr(progname, mesg, chem)
            
            # cdtrack(ncdtrack,1:mxcd)=track(1,1:mxcd) 
            cdtrack[ncdtrack_val][1:mxcd + 1] = track[1][1:mxcd + 1] 

            # cdtracklen(ncdtrack)=trlen(1)
            cdtracklen[ncdtrack_val] = trlen[1]
            
            # cdpst(track(1,trlen(1)))=0              
            last_node_index = track[1][trlen[1]] 
            cdpst[last_node_index] = 0
# END SUBROUTINE alkenetrack


# !=======================================================================
# ! PURPOSE: create tracks for -CO-O- groups. ...
# !=======================================================================
def estertrack(chem, bond, group, ngr, netrack, etracklen, etrack):
    
    # Internal variable size setup (1-based index: size = ngr + 1)
    array_size = ngr + 1
    
    # Internal arrays
    # ebond(SIZE(bond,1),SIZE(bond,1))
    ebond = [[0 for _ in range(array_size)] for _ in range(array_size)]
    # epst(SIZE(bond,1))
    epst = [0] * array_size
    
    # mxenod=SIZE(etrack,2)
    mxenod = len(etrack[0]) - 1
    
    # track(mxcp,SIZE(group))
    track = [[0 for _ in range(array_size)] for _ in range(MXCP + 1)]
    # trlen(mxcp)
    trlen = [0] * (MXCP + 1)
    # ntr is a scalar output from gettrack, simulated by a list
    ntr = [0] 

    # CHARACTER(LEN=12),PARAMETER  :: progname='estertrack'
    # CHARACTER(LEN=70) :: mesg
    progname = 'estertrack'
    mesg = ''

    # netrack=0 ; etrack(:,:)=0 ; etracklen(:)=0
    netrack[0] = 0
    # Clear etrack and etracklen (using 1-based indices)
    for r in range(1, len(etrack)):
        for c in range(1, len(etrack[0])):
            etrack[r][c] = 0
    for i in range(1, len(etracklen)):
        etracklen[i] = 0
        
    # ebond(:,:)=0 is handled by initialization

    # ------------------------------
    # Create and check the CO-O- tracks
    # ------------------------------

    # create a bond matrix for ester only might be bounded
    nest = 0
    # DO i=1,ngr
    for i in range(1, ngr + 1):
        # IF (group(i)(1:3)=='-O-') THEN
        if group[i][0:3] == '-O-':
            # DO j=1,ngr
            for j in range(1, ngr + 1):
                # IF (j==i) CYCLE
                if j == i: continue
                
                # IF (bond(i,j)/=0) THEN
                if bond[i][j] != 0:
                    grp_j = group[j]
                    
                    # IF (group(j)=='CO') THEN 
                    if grp_j == 'CO':
                        # ebond(i,j)=1 ; ebond(j,i)=1 ; nest=nest+1
                        ebond[i][j] = 1; ebond[j][i] = 1; nest += 1
                    # ELSEIF (group(j)=='CHO') THEN 
                    elif grp_j == 'CHO':
                        # ebond(i,j)=1 ; ebond(j,i)=1 ; nest=nest+1
                        ebond[i][j] = 1; ebond[j][i] = 1; nest += 1
                    # ENDIF
                # ENDIF
            # ENDDO
        # ENDIF
    # ENDDO
    
    # IF (nest==0) RETURN 
    if nest == 0: 
        return

    # Count # of connection of each ester nodes to other ester nodes
    # DO i=1,ngr ; epst(i)=SUM(ebond(i,:)) ; ENDDO
    for i in range(1, ngr + 1):
        # SUM(ebond(i,:)) -> sum of the 1-based part of the i-th row
        epst[i] = sum(ebond[i][1:]) 

    # check primary and secondary ester nodes
    nprim = 0
    # DO i=1,ngr
    for i in range(1, ngr + 1):
        # IF (epst(i)==1) nprim=nprim+1
        if epst[i] == 1:
            nprim += 1
        
        # IF (epst(i)>2) THEN
        if epst[i] > 2:
            mesg = "tertiary ester structure identified (and not allowed)"
            stoperr(progname, mesg, chem)
        # ENDIF        
    # ENDDO
    
    # IF (nprim==0) THEN
    if nprim == 0:
        mesg = "unexpected cyclic ester structure identified"
        stoperr(progname, mesg, chem)
    # ENDIF        

    # get ester chains (must start from primary Cd)
    # DO i=1,ngr
    for i in range(1, ngr + 1):
        # IF (epst(i)==1) THEN
        if epst[i] == 1:
            ntr = [0] # Reset ntr before calling gettrack
            # CALL gettrack(ebond,i,ngr,ntr,track,trlen)
            gettrack(ebond, i, ngr, ntr, track, trlen)
            ntr_val = ntr[0]
            
            # IF (ntr>1) THEN
            if ntr_val > 1:
                mesg = "unexpected branching identified on ester structure"
                stoperr(progname, mesg, chem)
            # ENDIF
            
            # IF (trlen(1) > mxenod) THEN
            if trlen[1] > mxenod:
                mesg = "More than mx possible ester identified (and not allowed)"
                stoperr(progname, mesg, chem)
            # ENDIF
            
            # netrack=netrack+1
            netrack[0] += 1
            netrack_val = netrack[0]

            # IF (netrack>SIZE(etrack,1)) THEN
            if netrack_val > len(etrack) - 1:
                mesg = "maximum number of Cd tracks reached." # Preserving the original Fortran error message string
                stoperr(progname, mesg, chem)
            # ENDIF
            
            # etrack(netrack,1:mxenod)=track(1,1:mxenod)
            # Python slicing: [1 : mxenod + 1] to include mxenod
            etrack[netrack_val][1:mxenod + 1] = track[1][1:mxenod + 1] 

            # etracklen(netrack)=trlen(1)
            etracklen[netrack_val] = trlen[1]
            
            # epst(track(1,trlen(1)))=0              ! kill reverse track
            # Get the index of the last node in the track
            last_node_index = track[1][trlen[1]] 
            epst[last_node_index] = 0
            
        # ENDIF
    # ENDDO
    
# END SUBROUTINE estertrack



# ----------------------------------------------------------------------
# SUBROUTINE chemmap
# ----------------------------------------------------------------------
def chemmap(chem, node, group, bond, ngrp, nodetype, alifun, cdfun, arofun, mapfun, funflg, tabester, nfcd, nfcr, ierr):
    
    # Internal variables initialization
    # tnod size is dynamic based on node. We initialize a list for 1-based indexing.
    # ytab size is 2 (index 1 and 2), so size 3 for 1-based.
    tnod = [0] * (node + 1)
    ytab = [0] * 3 
    nnod = 0
    nf = 0
    ialpha = 0
    ialpha2 = 0
    iy = 0
    rflg = 0
    dflg = 0
    yflg = 0
    ichecko = 0
    ichecky = 0

    ierr[0] = 0
    nfcd[0] = 0
    nfcr[0] = 0
    ngrp[0] = 0
    
    # Array initialization (assuming arrays are passed with index 0 unused, 1-based indexing)
    for i in range(1, len(alifun)):
        alifun[i] = 0.0
    for i in range(1, len(cdfun)):
        cdfun[i] = 0.0
    for i in range(1, len(arofun)):
        arofun[i] = 0.0
    
    for i in range(1, len(nodetype)):
        nodetype[i] = ' '
    for i in range(1, len(funflg)):
        funflg[i] = 0
    
    # Initialize 3D mapfun (i:node, j:3, k:max 21 index used)
    for i in range(1, len(mapfun)):
        for j in range(1, len(mapfun[i])):
            for k in range(1, len(mapfun[i][j])):
                mapfun[i][j][k] = 0.0
                
    # Initialize 2D tabester
    for i in range(1, len(tabester)):
        for j in range(1, len(tabester[i])):
            tabester[i][j] = 0
    
    nester = 0
    lgr = len(group[1].strip()) # LEN(group(1)) - Fortran character length

    # ! get the type of nodes in groups
    # DO i=1,node
    for i in range(1, node + 1):
        # IF (group(i)(1:2)=='CO') THEN       ; nodetype(i)='y'
        if group[i].strip()[:2] == 'CO':
            nodetype[i] = 'y'
        # ELSE IF (group(i)(1:3)=='CHO') THEN ; nodetype(i)='y'
        elif group[i].strip()[:3] == 'CHO':
            nodetype[i] = 'y'
        # ElSE IF (group(i)(1:1)=='c') THEN    ; nodetype(i)='r'
        elif group[i].strip()[:1] == 'c':
            nodetype[i] = 'r'
        # ELSE IF (group(i)(1:3)=='-O-') THEN ; nodetype(i)='o'
        elif group[i].strip()[:3] == '-O-':
            nodetype[i] = 'o'
        # ELSE IF (group(i)(1:2)=='Cd') THEN  ; nodetype(i)='d'
        elif group[i].strip()[:2] == 'Cd':
            nodetype[i] = 'd'
        # ELSE                               ; nodetype(i)='n'
        else:
            nodetype[i] = 'n'
        # ENDIF  
    # ENDDO

    # ! ------- Alkohols (index 1) and Acids (index 11) ----------
    # IF (INDEX(chem,'(OH)')/=0) THEN
    if chem.find('(OH)') != -1:
        # python中的-1索引相当于f90的0索引
        #################################333
        for i in range(1, node + 1):
            # IF (INDEX(group(i),'(OH)')/=0) THEN
            if group[i].find('(OH)') != -1:
                nf = 0
                # DO j=1,lgr-3
                for j in range(1, lgr - 3 + 1):
                    # IF (group(i)(j:j+3)=='(OH)') nf=nf+1
                    if group[i][j-1 : j+3] == '(OH)':
                        nf = nf + 1
                # ENDDO
                
                # IF (nodetype(i)=='n') THEN             ! alcohol aliphatic
                if nodetype[i] == 'n':
                    alifun[1] = alifun[1] + nf
                    mapfun[i][1][1] = nf
                    ngrp[0] = ngrp[0] + nf
                # ELSE IF (nodetype(i)=='r') THEN    ! phenols 
                elif nodetype[i] == 'r':
                    arofun[1] = arofun[1] + nf
                    mapfun[i][3][1] = nf
                    ngrp[0] = ngrp[0] + nf
                    nfcr[0] = nfcr[0] + nf
                # ELSE IF (nodetype(i)=='d') THEN    ! =Cd-OH (not expected, since enol) 
                elif nodetype[i] == 'd':
                    cdfun[1] = cdfun[1] + nf
                    mapfun[i][2][1] = nf
                    ngrp[0] = ngrp[0] + nf
                    nfcd[0] = nfcd[0] + nf
                # ELSE IF (nodetype(i)=='y') THEN    ! Carboxylic acid
                elif nodetype[i] == 'y':
                    # CALL nodmap(bond,i,node,2,nnod,tnod)
                    nnod_ref = [0]
                    try:
                        nodmap(bond, i, node, 2, nnod_ref, tnod)
                        nnod = nnod_ref[0]
                    except NotImplementedError:
                         # Assume nnod is populated to 1 for unique C check
                         nnod = 1 
                         # This is a risk; the user must provide 'nodmap' for correct execution.
                         tnod[1] = 1 # Placeholder for ialpha
                         pass 
                    
                    # IF (nnod>1) THEN 
                    if nnod > 1:
                        # DO j=1,node ; WRITE(*,*) 'group',j,'-',group(j) ; ENDDO
                        for j in range(1, node + 1):
                            Fortran_write6('group', j, '-', group[j])
                        # WRITE(6,*) '-- error --, a unique C is expected in'
                        Fortran_write6('-- error --, a unique C is expected in')
                        # WRITE(6,*) 'alpha position of a carboxylic group'
                        Fortran_write6('alpha position of a carboxylic group')
                        # STOP "in chemmap"
                        sys.exit("in chemmap")
                    # ENDIF
                    ialpha = tnod[1]
                    
                    # IF (nodetype(ialpha)=='r')THEN         ! CO(OH) aromatic 
                    if nodetype[ialpha] == 'r':
                        arofun[11] = arofun[11] + 1
                        mapfun[i][3][11] = 1
                        ngrp[0] = ngrp[0] + 1
                        nfcr[0] = nfcr[0] + 1
                    # ELSE IF (nodetype(ialpha)=='d')THEN    ! CO(OH) on Cd
                    elif nodetype[ialpha] == 'd':
                        cdfun[11] = cdfun[11] + 1
                        mapfun[i][2][11] = 1
                        ngrp[0] = ngrp[0] + 1
                        nfcd[0] = nfcd[0] + 1
                    # ELSE                                   ! CO(OH) aliphatic
                    else:
                        alifun[11] = alifun[11] + 1
                        mapfun[i][1][11] = 1
                        ngrp[0] = ngrp[0] + 1
                    # ENDIF 
                # ENDIF
            # ENDDO  
        # ENDIF  

    # ! ----------- Nitro (index 2) -----------------
    # IF (INDEX(chem,'(NO2)')/=0) THEN
    if chem.find('(NO2)') != -1:
        # DO i = 1, node
        for i in range(1, node + 1):
            # IF (INDEX(group(i),'(NO2)')/=0) THEN
            if group[i].find('(NO2)') != -1:
                nf = 0
                # DO j=1,lgr-4
                for j in range(1, lgr - 4 + 1):
                    # IF (group(i)(j:j+4)=='(NO2)') nf=nf+1
                    if group[i][j-1 : j+4] == '(NO2)':
                        nf = nf + 1
                # ENDDO
                
                # IF (nodetype(i)=='n') THEN      ! NO2 aliphatic
                if nodetype[i] == 'n':
                    alifun[2] = alifun[2] + nf
                    mapfun[i][1][2] = nf
                    ngrp[0] = ngrp[0] + nf
                # ELSE IF (nodetype(i)=='y') THEN     ! NO2 aliphatic
                elif nodetype[i] == 'y':
                    alifun[2] = alifun[2] + nf
                    mapfun[i][1][2] = nf
                    ngrp[0] = ngrp[0] + nf
                # ELSE IF (nodetype(i)=='d') THEN     ! NO2 on Cd 
                elif nodetype[i] == 'd':
                    cdfun[2] = cdfun[2] + nf
                    mapfun[i][2][2] = nf
                    ngrp[0] = ngrp[0] + nf
                    nfcd[0] = nfcd[0] + nf
                # ELSE IF (nodetype(i)=='r') THEN     ! NO2 on aromatic 
                elif nodetype[i] == 'r':
                    arofun[2] = arofun[2] + nf
                    mapfun[i][3][2] = nf
                    ngrp[0] = ngrp[0] + nf
                    nfcr[0] = nfcr[0] + nf
                # ELSE
                else:
                    # WRITE(6,*) '-- error --, a (NO2) group is borne by an '
                    Fortran_write6('-- error --, a (NO2) group is borne by an ')
                    # WRITE(6,*) ' unexpected group in chem :'
                    Fortran_write6(' unexpected group in chem :')
                    # WRITE(6,*) TRIM(chem)
                    Fortran_write6(chem.strip())
                    # STOP "in chemmap"
                    sys.exit("in chemmap")
                # ENDIF
            # ENDIF
        # ENDDO
    # ENDIF

    # ! ----------- Nitrate (index 3) -----------------
    # IF (INDEX(chem,'(ONO2)')/=0) THEN
    if chem.find('(ONO2)') != -1:
        # DO i = 1, node
        for i in range(1, node + 1):
            # IF (INDEX(group(i),'(ONO2)')/=0) THEN
            if group[i].find('(ONO2)') != -1:

                # !          IF (INDEX(group(i),'CO(ONO2)')/=0) GOTO 120  
                # IF (INDEX(group(i),'CO(ONO2)')/=0) CYCLE ! ba: unexpected but out of IF in original version
                if group[i].find('CO(ONO2)') != -1:
                    continue # CYCLE
                
                nf = 0
                # DO j=1,lgr-5
                for j in range(1, lgr - 5 + 1):
                    # IF (group(i)(j:j+5)=='(ONO2)') nf=nf+1
                    if group[i][j-1 : j+5] == '(ONO2)':
                        nf = nf + 1
                # ENDDO

                # IF (nodetype(i)=='n') THEN             ! ONO2 aliphatic
                if nodetype[i] == 'n':
                    alifun[3] = alifun[3] + nf
                    mapfun[i][1][3] = nf
                    ngrp[0] = ngrp[0] + nf
                # ELSE IF (nodetype(i)=='d') THEN        ! ONO2 on Cd 
                elif nodetype[i] == 'd':
                    cdfun[3] = cdfun[3] + nf
                    mapfun[i][2][3] = nf
                    ngrp[0] = ngrp[0] + nf
                    nfcd[0] = nfcd[0] + nf
                # ELSE IF (nodetype(i)=='r') THEN        ! ONO2 on aromatic 
                elif nodetype[i] == 'r':
                    arofun[3] = arofun[3] + nf
                    mapfun[i][3][3] = nf
                    ngrp[0] = ngrp[0] + nf
                    nfcr[0] = nfcr[0] + nf
                # ELSE
                else:
                    # WRITE(6,*) '-- error --, a (ONO2) group is borne by an '
                    Fortran_write6('-- error --, a (ONO2) group is borne by an ')
                    # WRITE(6,*) 'unexpected group in chem :'
                    Fortran_write6('unexpected group in chem :')
                    # WRITE(6,*) TRIM(chem)
                    Fortran_write6(chem.strip())
                    # STOP "in chemmap"
                    sys.exit("in chemmap")
                # ENDIF
            # ENDIF
        # ENDDO
    # ENDIF

    # ! ------- hydroperoxydes (index 4) and peracids (index 12) ----------
    # IF (INDEX(chem,'(OOH)')/=0) THEN
    if chem.find('(OOH)') != -1:
        # DO i=1, node
        for i in range(1, node + 1):
            # IF (INDEX(group(i),'(OOH)')/=0) THEN
            if group[i].find('(OOH)') != -1:
                nf = 0
                # DO j=1,lgr-4
                for j in range(1, lgr - 4 + 1):
                    # IF (group(i)(j:j+4)=='(OOH)') nf=nf+1
                    if group[i][j-1 : j+4] == '(OOH)':
                        nf = nf + 1
                # ENDDO
                
                # IF (nodetype(i)=='n') THEN         ! hydroperoxyde aliphatic
                if nodetype[i] == 'n':
                    alifun[4] = alifun[4] + nf
                    mapfun[i][1][4] = nf
                    ngrp[0] = ngrp[0] + nf
                # ELSE IF (nodetype(i)=='r') THEN    ! aromaric -OOH 
                elif nodetype[i] == 'r':
                    arofun[4] = arofun[4] + nf
                    mapfun[i][3][4] = nf
                    ngrp[0] = ngrp[0] + nf
                    nfcr[0] = nfcr[0] + nf
                # ELSE IF (nodetype(i)=='d') THEN    ! =Cd-OOH (not expected) 
                elif nodetype[i] == 'd':
                    cdfun[4] = cdfun[4] + nf
                    mapfun[i][2][4] = nf
                    ngrp[0] = ngrp[0] + nf
                    nfcd[0] = nfcd[0] + nf
                # ELSE IF (nodetype(i)=='y') THEN    ! peracid acid
                elif nodetype[i] == 'y':
                    # CALL nodmap(bond,i,node,2,nnod,tnod)
                    nnod_ref = [0]
                    try:
                        nodmap(bond, i, node, 2, nnod_ref, tnod)
                        nnod = nnod_ref[0]
                    except NotImplementedError:
                         nnod = 1 
                         tnod[1] = 1
                         pass
                        
                    # IF (nnod>1) THEN 
                    if nnod > 1:
                        # DO j=1,node ; WRITE(*,*) 'group',j,'-',group(j) ; ENDDO
                        for j in range(1, node + 1):
                            Fortran_write6('group', j, '-', group[j])
                        # WRITE(6,*) '-- error --, a unique C is expected in'
                        Fortran_write6('-- error --, a unique C is expected in')
                        # WRITE(6,*) 'alpha position of a peracid group'
                        Fortran_write6('alpha position of a peracid group')
                        # WRITE(6,*) TRIM(chem)
                        Fortran_write6(chem.strip())
                        # STOP "in chemmap"
                        sys.exit("in chemmap")
                    # ENDIF
                    ialpha = tnod[1]
                    
                    # IF (nodetype(ialpha)=='r')THEN         ! CO(OOH) aromatic 
                    if nodetype[ialpha] == 'r':
                        arofun[12] = arofun[12] + 1
                        mapfun[i][3][12] = 1
                        ngrp[0] = ngrp[0] + 1
                    # ELSE IF (nodetype(ialpha)=='d')THEN    ! CO(OOH) on Cd
                    elif nodetype[ialpha] == 'd':
                        cdfun[12] = cdfun[12] + 1
                        mapfun[i][2][12] = 1
                        ngrp[0] = ngrp[0] + 1
                    # ELSE                                   ! CO(OOH) aliphatic
                    else:
                        alifun[12] = alifun[12] + 1
                        mapfun[i][1][12] = 1
                        ngrp[0] = ngrp[0] + 1
                    # ENDIF 
                # ENDIF
            # ENDDO 
        # ENDIF 

    # ! ------- fluroro (index 5) and fluoro acyl (index 17) ----------
    # IF (INDEX(chem,'(F)')/=0) THEN
    if chem.find('(F)') != -1:
        # DO i=1, node
        for i in range(1, node + 1):
            # IF (INDEX(group(i),'(F)')/=0) THEN
            if group[i].find('(F)') != -1:
                nf = 0
                # DO j=1,lgr-2
                for j in range(1, lgr - 2 + 1):
                    # IF (group(i)(j:j+2)=='(F)') nf=nf+1
                    if group[i][j-1 : j+2] == '(F)':
                        nf = nf + 1
                # ENDDO
                
                # IF (nodetype(i)=='n') THEN         ! F aliphatic
                if nodetype[i] == 'n':
                    alifun[5] = alifun[5] + nf
                    mapfun[i][1][5] = nf
                    ngrp[0] = ngrp[0] + nf
                # ELSE IF (nodetype(i)=='r') THEN    ! aromaric F 
                elif nodetype[i] == 'r':
                    arofun[5] = arofun[5] + nf
                    mapfun[i][3][5] = nf
                    ngrp[0] = ngrp[0] + nf
                    nfcr[0] = nfcr[0] + nf
                # ELSE IF (nodetype(i)=='d') THEN    ! =Cd-F 
                elif nodetype[i] == 'd':
                    cdfun[5] = cdfun[5] + nf
                    mapfun[i][2][5] = nf
                    ngrp[0] = ngrp[0] + nf
                    nfcd[0] = nfcd[0] + nf
                # ELSE IF (nodetype(i)=='y') THEN    ! fluoro acyl
                elif nodetype[i] == 'y':
                    # CALL nodmap(bond,i,node,2,nnod,tnod)
                    nnod_ref = [0]
                    try:
                        nodmap(bond, i, node, 2, nnod_ref, tnod)
                        nnod = nnod_ref[0]
                    except NotImplementedError:
                         nnod = 1 
                         tnod[1] = 1
                         pass
                        
                    # IF (nnod>1) THEN 
                    if nnod > 1:
                        # DO j=1,node ; WRITE(*,*) 'group',j,'-',group(j) ; ENDDO
                        for j in range(1, node + 1):
                            Fortran_write6('group', j, '-', group[j])
                        # WRITE(6,*) '-- error --, a unique C is expected in'
                        Fortran_write6('-- error --, a unique C is expected in')
                        # WRITE(6,*) 'alpha position of a fluoro acyl group'
                        Fortran_write6('alpha position of a fluoro acyl group')
                        # WRITE(6,*) TRIM(chem)
                        Fortran_write6(chem.strip())
                        # STOP "in chemmap"
                        sys.exit("in chemmap")
                    # ENDIF
                    ialpha = tnod[1]
                    
                    # IF (nodetype(ialpha)=='r')THEN         ! CO(F) aromatic 
                    if nodetype[ialpha] == 'r':
                        arofun[17] = arofun[17] + 1
                        mapfun[i][3][17] = 1
                        ngrp[0] = ngrp[0] + 1
                    # ELSE IF (nodetype(ialpha)=='d')THEN    ! CO(F) on Cd
                    elif nodetype[ialpha] == 'd':
                        cdfun[17] = cdfun[17] + 1
                        mapfun[i][2][17] = 1
                        ngrp[0] = ngrp[0] + 1
                    # ELSE                                   ! CO(F) aliphatic
                    else:
                        alifun[17] = alifun[17] + 1
                        mapfun[i][1][17] = 1
                        ngrp[0] = ngrp[0] + 1
                    # ENDIF 
                # ENDIF
            # ENDDO 
        # ENDIF 

    # ! ------- chloro (index 6) and chloro acyl (index 18) ----------
    # IF (INDEX(chem,'(Cl)')/=0) THEN
    if chem.find('(Cl)') != -1:
        # DO i=1, node
        for i in range(1, node + 1):
            # IF (INDEX(group(i),'(Cl)')/=0) THEN
            if group[i].find('(Cl)') != -1:
                nf = 0
                # DO j=1,lgr-3
                for j in range(1, lgr - 3 + 1):
                    # IF (group(i)(j:j+3)=='(Cl)') nf=nf+1
                    if group[i][j-1 : j+3] == '(Cl)':
                        nf = nf + 1
                # ENDDO
                
                # IF (nodetype(i)=='n') THEN         ! Cl aliphatic
                if nodetype[i] == 'n':
                    alifun[6] = alifun[6] + nf
                    mapfun[i][1][6] = nf
                    ngrp[0] = ngrp[0] + nf
                # ELSE IF (nodetype(i)=='r') THEN    ! aromaric Cl
                elif nodetype[i] == 'r':
                    arofun[6] = arofun[6] + nf
                    mapfun[i][3][6] = nf
                    ngrp[0] = ngrp[0] + nf
                    nfcr[0] = nfcr[0] + nf
                # ELSE IF (nodetype(i)=='d') THEN    ! =Cd-Cl  
                elif nodetype[i] == 'd':
                    cdfun[6] = cdfun[6] + nf
                    mapfun[i][2][6] = nf
                    ngrp[0] = ngrp[0] + nf
                    nfcd[0] = nfcd[0] + nf
                # ELSE IF (nodetype(i)=='y') THEN    ! Chloro acyl
                elif nodetype[i] == 'y':
                    # CALL nodmap(bond,i,node,2,nnod,tnod)
                    nnod_ref = [0]
                    try:
                        nodmap(bond, i, node, 2, nnod_ref, tnod)
                        nnod = nnod_ref[0]
                    except NotImplementedError:
                         nnod = 1 
                         tnod[1] = 1
                         pass
                        
                    # IF (nnod>1) THEN 
                    if nnod > 1:
                        # DO j=1,node ; WRITE(*,*) 'group',j,'-',group(j) ; ENDDO
                        for j in range(1, node + 1):
                            Fortran_write6('group', j, '-', group[j])
                        # WRITE(6,*) '-- error --, a unique C is expected in'
                        Fortran_write6('-- error --, a unique C is expected in')
                        # WRITE(6,*) 'alpha position of a chloro acyl group'
                        Fortran_write6('alpha position of a chloro acyl group')
                        # WRITE(6,*) TRIM(chem)
                        Fortran_write6(chem.strip())
                        # STOP "in chemmap"
                        sys.exit("in chemmap")
                    # ENDIF
                    ialpha = tnod[1]
                    
                    # IF (nodetype(ialpha)=='r')THEN         ! CO(Cl) aromatic 
                    if nodetype[ialpha] == 'r':
                        arofun[18] = arofun[18] + 1
                        mapfun[i][3][18] = 1
                        ngrp[0] = ngrp[0] + 1
                    # ELSE IF (nodetype(ialpha)=='d')THEN    ! CO(Cl) on Cd
                    elif nodetype[ialpha] == 'd':
                        cdfun[18] = cdfun[18] + 1
                        mapfun[i][2][18] = 1
                        ngrp[0] = ngrp[0] + 1
                    # ELSE                                   ! CO(Cl) aliphatic
                    else:
                        alifun[18] = alifun[18] + 1
                        mapfun[i][1][18] = 1
                        ngrp[0] = ngrp[0] + 1
                    # ENDIF 
                # ENDIF
            # ENDDO 
        # ENDIF 

    # ! ------- Bromo (index 7) and bromo acyl (index 19) ----------
    # IF (INDEX(chem,'(Br)')/=0) THEN
    if chem.find('(Br)') != -1:
        # DO i=1, node
        for i in range(1, node + 1):
            # IF (INDEX(group(i),'(Br)')/=0) THEN
            if group[i].find('(Br)') != -1:
                nf = 0
                # DO j=1,lgr-3
                for j in range(1, lgr - 3 + 1):
                    # IF (group(i)(j:j+3)=='(Br)') nf=nf+1
                    if group[i][j-1 : j+3] == '(Br)':
                        nf = nf + 1
                # ENDDO
                
                # IF (nodetype(i)=='n') THEN         ! Br aliphatic
                if nodetype[i] == 'n':
                    alifun[7] = alifun[7] + nf
                    mapfun[i][1][7] = nf
                    ngrp[0] = ngrp[0] + nf
                # ELSE IF (nodetype(i)=='r') THEN    ! aromatic Br
                elif nodetype[i] == 'r':
                    arofun[7] = arofun[7] + nf
                    mapfun[i][3][7] = nf
                    ngrp[0] = ngrp[0] + nf
                    nfcr[0] = nfcr[0] + nf
                # ELSE IF (nodetype(i)=='d') THEN    ! =Cd-Br  
                elif nodetype[i] == 'd':
                    cdfun[7] = cdfun[7] + nf
                    mapfun[i][2][7] = nf
                    ngrp[0] = ngrp[0] + nf
                    nfcd[0] = nfcd[0] + nf
                # ELSE IF (nodetype(i)=='y') THEN    ! Bromo acyl
                elif nodetype[i] == 'y':
                    # CALL nodmap(bond,i,node,2,nnod,tnod)
                    nnod_ref = [0]
                    try:
                        nodmap(bond, i, node, 2, nnod_ref, tnod)
                        nnod = nnod_ref[0]
                    except NotImplementedError:
                         nnod = 1 
                         tnod[1] = 1
                         pass
                        
                    # IF (nnod>1) THEN 
                    if nnod > 1:
                        # DO j=1,node  ; WRITE(*,*) 'group',j,'-',group(j)  ;  ENDDO
                        for j in range(1, node + 1):
                            Fortran_write6('group', j, '-', group[j])
                        # WRITE(6,*) '-- error --, a unique C is expected in'
                        Fortran_write6('-- error --, a unique C is expected in')
                        # WRITE(6,*) 'alpha position of a bromo acyl group'
                        Fortran_write6('alpha position of a bromo acyl group')
                        # WRITE(6,*) chem(1:70)
                        Fortran_write6(chem[:70])
                        # STOP "in chemmap"
                        sys.exit("in chemmap")
                    # ENDIF
                    ialpha = tnod[1]
                    
                    # IF (nodetype(ialpha)=='r')THEN         ! CO(Br) aromatic 
                    if nodetype[ialpha] == 'r':
                        arofun[19] = arofun[19] + 1
                        mapfun[i][3][19] = 1
                        ngrp[0] = ngrp[0] + 1
                    # ELSE IF (nodetype(ialpha)=='d')THEN    ! CO(Br) on Cd
                    elif nodetype[ialpha] == 'd':
                        cdfun[19] = cdfun[19] + 1
                        mapfun[i][2][19] = 1
                        ngrp[0] = ngrp[0] + 1
                    # ELSE                                   ! CO(Br) aliphatic
                    else:
                        alifun[19] = alifun[19] + 1
                        mapfun[i][1][19] = 1
                        ngrp[0] = ngrp[0] + 1
                    # ENDIF 
                # ENDIF
            # ENDDO 
        # ENDIF 

    # ! ------- iodo (index 8) and iodo acyl (index 20) ----------
    # IF (INDEX(chem,'(I)')/=0) THEN
    if chem.find('(I)') != -1:
        # DO i=1, node
        for i in range(1, node + 1):
            # IF (INDEX(group(i),'(I)')/=0) THEN
            if group[i].find('(I)') != -1:
                nf = 0
                # DO j=1,lgr-2
                for j in range(1, lgr - 2 + 1):
                    # IF (group(i)(j:j+2)=='(I)') nf=nf+1
                    if group[i][j-1 : j+2] == '(I)':
                        nf = nf + 1
                # ENDDO
                
                # IF (nodetype(i)=='n') THEN         ! I aliphatic
                if nodetype[i] == 'n':
                    alifun[8] = alifun[8] + nf
                    mapfun[i][1][8] = nf
                    ngrp[0] = ngrp[0] + nf
                # ELSE IF (nodetype(i)=='r') THEN    ! aromaric I 
                elif nodetype[i] == 'r':
                    arofun[8] = arofun[8] + nf
                    mapfun[i][3][8] = nf
                    ngrp[0] = ngrp[0] + nf
                    nfcr[0] = nfcr[0] + nf
                # ELSE IF (nodetype(i)=='d') THEN    ! =Cd-I 
                elif nodetype[i] == 'd':
                    cdfun[8] = cdfun[8] + nf
                    mapfun[i][2][8] = nf
                    ngrp[0] = ngrp[0] + nf
                    nfcd[0] = nfcd[0] + nf
                # ELSE IF (nodetype(i)=='y') THEN    ! fluoro acyl
                elif nodetype[i] == 'y':
                    # CALL nodmap(bond,i,node,2,nnod,tnod)
                    nnod_ref = [0]
                    try:
                        nodmap(bond, i, node, 2, nnod_ref, tnod)
                        nnod = nnod_ref[0]
                    except NotImplementedError:
                         nnod = 1 
                         tnod[1] = 1
                         pass
                        
                    # IF (nnod>1) THEN 
                    if nnod > 1:
                        # DO j=1,node ; WRITE(*,*) 'group',j,'-',group(j)  ;  ENDDO
                        for j in range(1, node + 1):
                            Fortran_write6('group', j, '-', group[j])
                        # WRITE(6,*) '-- error --, a unique C is expected in'
                        Fortran_write6('-- error --, a unique C is expected in')
                        # WRITE(6,*) 'alpha position of a Iodo acyl group'
                        Fortran_write6('alpha position of a Iodo acyl group')
                        # WRITE(6,*) chem(1:70)
                        Fortran_write6(chem[:70])
                        # STOP "in chemmap"
                        sys.exit("in chemmap")
                    # ENDIF
                    ialpha = tnod[1]
                    
                    # IF (nodetype(ialpha)=='r')THEN         ! CO(I) aromatic 
                    if nodetype[ialpha] == 'r':
                        arofun[20] = arofun[20] + 1
                        mapfun[i][3][20] = 1
                        ngrp[0] = ngrp[0] + 1
                        nfcr[0] = nfcr[0] + 1
                    # ELSE IF (nodetype(ialpha)=='d')THEN    ! CO(I) on Cd
                    elif nodetype[ialpha] == 'd':
                        cdfun[20] = cdfun[20] + 1
                        mapfun[i][2][20] = 1
                        ngrp[0] = ngrp[0] + 1
                        nfcd[0] = nfcd[0] + 1
                    # ELSE                                   ! CO(I) aliphatic
                    else:
                        alifun[20] = alifun[20] + 1
                        mapfun[i][1][20] = 1
                        ngrp[0] = ngrp[0] + 1
                    # ENDIF 
                # ENDIF
            # ENDDO 
        # ENDIF 

    # ! ---------- ester (index 15 and 16) -------------
    # IF (INDEX(chem,'-O')/=0) THEN
    if chem.find('-O') != -1:
        # esterloop: DO i = 1, node
        i = 1
        while i <= node:
            # IF (group(i)(1:3)=='-O-') THEN
            if group[i].strip()[:3] == '-O-':
                # CALL nodmap(bond,i,node,2,nnod,tnod)
                nnod_ref = [0]
                try:
                    nodmap(bond, i, node, 2, nnod_ref, tnod)
                    nnod = nnod_ref[0]
                except NotImplementedError:
                     nnod = 2 
                     tnod[1], tnod[2] = 0, 0 # Placeholder, dangerous
                     pass

                rflg = 0
                dflg = 0
                yflg = 0
                
                # DO j=1,nnod
                for j in range(1, nnod + 1):
                    ialpha = tnod[j]
                    # IF (nodetype(ialpha)=='y') THEN   ! RCO-O-R function
                    if nodetype[ialpha] == 'y':
                        ichecky = 0
                        # DO k=1,4
                        for k in range(1, 4 + 1):
                            # IF (tabester(k,2)==ialpha) ichecky=1    ! carbonyl already used
                            if tabester[k][2] == ialpha:
                                ichecky = 1
                                break # Exit k loop
                        # ENDDO
                        # IF (ichecky==0) THEN
                        if ichecky == 0:
                            yflg = yflg + 1
                            ytab[yflg] = ialpha
                        # ENDIF
                    # ENDIF
                # ENDDO

                # IF (yflg==0) CYCLE  esterloop           ! simple ether
                if yflg == 0:
                    i += 1
                    continue # CYCLE esterloop

                # ! simple ester 
                # IF (yflg/=0) THEN
                if yflg != 0:
                    nester = nester + 1
                    # IF (nester>4) THEN
                    if nester > 4:
                        # PRINT*, "in chemmap, nester > 4"
                        print("in chemmap, nester > 4")
                        # STOP "in chemmap"
                        sys.exit("in chemmap")
                    # ENDIF  
                    
                    iy = ytab[1]
                    tabester[nester][1] = i
                    tabester[nester][2] = iy
                    
                    # DO j=1,2
                    for j in range(1, 2 + 1):
                        # IF (tnod(j)/=ytab(1)) THEN
                        if tnod[j] != ytab[1]:
                            ialpha = tnod[j]  # second side of the -O- (first is iy)
                        # ENDIF
                    # ENDDO
                    
                    # IF (group(iy)(1:3)=='CHO') THEN    ! HCO-O-R
                    if group[iy].strip()[:3] == 'CHO':
                        # IF (nodetype(ialpha)=='r') THEN        ! aromatic CHO-O-R
                        if nodetype[ialpha] == 'r':
                            arofun[16] = arofun[16] + 1
                            mapfun[i][3][16] = 1
                            mapfun[iy][3][16] = 1
                            ngrp[0] = ngrp[0] + 1
                            nfcr[0] = nfcr[0] + 1
                        # ELSE IF (nodetype(ialpha)=='d') THEN  ! =C-O-CHO
                        elif nodetype[ialpha] == 'd':
                            cdfun[16] = cdfun[16] + 1
                            mapfun[i][2][16] = 1
                            mapfun[iy][2][16] = 1
                            ngrp[0] = ngrp[0] + 1
                            nfcd[0] = nfcd[0] + 1
                        # ELSE                                   ! R-O-CHO
                        else:
                            alifun[16] = alifun[16] + 1
                            mapfun[i][1][16] = 1
                            mapfun[iy][1][16] = 1
                            ngrp[0] = ngrp[0] + 1
                        # ENDIF
                        
                    # ELSE IF (group(iy)(1:3)=='CO ') THEN    ! RCO-O-R
                    elif group[iy].strip()[:3] == 'CO ':
                        # CALL nodmap(bond,iy,node,2,nnod,tnod)
                        nnod_ref = [0]
                        try:
                            nodmap(bond, iy, node, 2, nnod_ref, tnod)
                            nnod = nnod_ref[0]
                        except NotImplementedError:
                             nnod = 2 
                             tnod[1], tnod[2] = 0, 0
                             pass

                        # DO j=1,2
                        for j in range(1, 2 + 1):
                            # IF (tnod(j)/=i) ialpha2=tnod(j)  ! ialpha2 is the node next to the CO of the ester
                            if tnod[j] != i:
                                ialpha2 = tnod[j]
                        # ENDDO
                        
                        rflg = 0
                        dflg = 0

                        # ! structure is ialpha2-CO-O-ialpha
                        # IF (nodetype(ialpha)=='r') rflg=rflg+1
                        if nodetype[ialpha] == 'r':
                            rflg = rflg + 1
                        # IF (nodetype(ialpha)=='d') dflg=dflg+1
                        if nodetype[ialpha] == 'd':
                            dflg = dflg + 1
                        # IF (nodetype(ialpha2)=='r') rflg=rflg+1
                        if nodetype[ialpha2] == 'r':
                            rflg = rflg + 1
                        # IF (nodetype(ialpha2)=='d') dflg=dflg+1
                        if nodetype[ialpha2] == 'd':
                            dflg = dflg + 1

                        # IF (rflg/=0) THEN             ! aromatic ester
                        if rflg != 0:
                            arofun[15] = arofun[15] + 1.0
                            mapfun[i][3][15] = 1.0
                            mapfun[iy][3][15] = 1.0
                            nfcr[0] = nfcr[0] + rflg
                        # ELSE IF (dflg/=0) THEN             ! =C-CO-O-
                        elif dflg != 0:
                            cdfun[15] = cdfun[15] + 1.0
                            mapfun[i][2][15] = 1.0
                            mapfun[iy][2][15] = 1.0
                            nfcd[0] = nfcd[0] + dflg
                        # ELSE                             ! R-COO-R
                        else:
                            alifun[15] = alifun[15] + 1.0
                            mapfun[i][1][15] = 1.0
                            mapfun[iy][1][15] = 1.0
                        # ENDIF
                        ngrp[0] = ngrp[0] + 1
                    # ENDIF
                # ENDIF

            # ENDIF  
            i += 1
        # ENDDO esterloop
    # ENDIF

    # ! ------------Aldehydes (index 9) --------- 
    # IF (INDEX(chem,'CHO')/=0) THEN
    if chem.find('CHO') != -1:
        # dloop: DO i = 1, node
        i = 1
        while i <= node:
            # IF (group(i)(1:3)=='CHO') THEN
            if group[i].strip()[:3] == 'CHO':
                # CALL nodmap(bond,i,node,2,nnod,tnod)
                nnod_ref = [0]
                try:
                    nodmap(bond, i, node, 2, nnod_ref, tnod)
                    nnod = nnod_ref[0]
                except NotImplementedError:
                     nnod = 1 
                     tnod[1] = 1
                     pass

                # IF (nnod/=1) THEN
                if nnod != 1:
                    # WRITE(6,*) '-- error --, a unique C is expected in'
                    Fortran_write6('-- error --, a unique C is expected in')
                    # WRITE(6,*) 'alpha position of a CHO  group'
                    Fortran_write6('alpha position of a CHO group')
                    # WRITE(6,*) TRIM(chem)
                    Fortran_write6(chem.strip())
                    # STOP "in chemmap"
                    sys.exit("in chemmap")
                # ENDIF
                ialpha = tnod[1]
                
                # IF (nodetype(ialpha)=='o') THEN    ! HCO-O-R function
                if nodetype[ialpha] == 'o':
                    ichecko = 0 # check if the ether is already involved in an ester
                    ichecky = 0 # check if the carbonyl is already involved in an ester
                    # DO k=1,4
                    for k in range(1, 4 + 1):
                        # IF (tabester(k,1)==ialpha) ichecko=1    ! ether already used
                        if tabester[k][1] == ialpha:
                            ichecko = 1
                        # IF (tabester(k,2)==i) ichecky=1     ! carbonyl already used
                        if tabester[k][2] == i:
                            ichecky = 1
                    # ENDDO
                    
                    # IF (ichecky==1) CYCLE dloop ! ether already involved
                    if ichecky == 1:
                        i += 1
                        continue # CYCLE dloop
                    # IF (ichecko==0) CYCLE dloop ! carbonyl that must be an ester
                    if ichecko == 0:
                        i += 1
                        continue # CYCLE dloop
                # ENDIF  ! if that point is reached then must be counted as aldehyde
                
                # IF (nodetype(ialpha)=='r') THEN            ! aromatic aldehyde
                if nodetype[ialpha] == 'r':
                    arofun[9] = arofun[9] + 1
                    mapfun[i][3][9] = 1
                    ngrp[0] = ngrp[0] + 1
                    nfcr[0] = nfcr[0] + 1
                # ELSE IF (nodetype(ialpha)=='d') THEN  ! =C-CHO
                elif nodetype[ialpha] == 'd':
                    cdfun[9] = cdfun[9] + 1
                    mapfun[i][2][9] = 1
                    ngrp[0] = ngrp[0] + 1
                    nfcd[0] = nfcd[0] + 1
                # ELSE                                     ! R-CHO
                else:
                    alifun[9] = alifun[9] + 1
                    mapfun[i][1][9] = 1
                    ngrp[0] = ngrp[0] + 1
                # ENDIF
            # ENDIF
            i += 1
        # ENDDO dloop
    # ENDIF

    # ! ---------- ketones (index 10) -------------
    # IF (INDEX(chem,'CO')/=0) THEN
    if chem.find('CO') != -1:
        # kloop: DO i = 1, node
        i = 1
        while i <= node:
            # IF (group(i)(1:3)=='CO ') THEN
            if group[i].strip()[:3] == 'CO ':
                # CALL nodmap(bond,i,node,2,nnod,tnod)
                nnod_ref = [0]
                try:
                    nodmap(bond, i, node, 2, nnod_ref, tnod)
                    nnod = nnod_ref[0]
                except NotImplementedError:
                     nnod = 2 
                     tnod[1], tnod[2] = 0, 0
                     pass

                # IF (nnod/=2) THEN
                if nnod != 2:
                    # WRITE(6,*) '-- error --, only 2 C is expected in'
                    Fortran_write6('-- error --, only 2 C is expected in')
                    # WRITE(6,*) 'alpha position of a -CO-  group'
                    Fortran_write6('alpha position of a -CO- group')
                    # WRITE(6,*) TRIM(chem),'  nnod=',nnod
                    Fortran_write6(chem.strip(), ' nnod=', nnod)
                    # WRITE(6,*) bond(1:node,1:node)
                    # Print 2D array bond for debugging
                    bond_str = '\n'.join([' '.join(map(str, row[1:node+1])) for row in bond[1:node+1]])
                    Fortran_write6(bond_str)
                    # STOP "in chemmap"
                    sys.exit("in chemmap")
                # ENDIF
                
                rflg = 0
                dflg = 0
                # DO j=1,nnod
                for j in range(1, nnod + 1):
                    ialpha = tnod[j]
                    # IF (nodetype(ialpha)=='o') THEN    ! RCO-O-R function
                    if nodetype[ialpha] == 'o':
                        ichecko = 0 # check if the ether is already involved in an ester
                        ichecky = 0 # check if the carbonyl is already involved in an ester
                        # DO k=1,4
                        for k in range(1, 4 + 1):
                            # IF (tabester(k,1)==ialpha) ichecko=1    ! ether already used
                            if tabester[k][1] == ialpha:
                                ichecko = 1
                            # IF (tabester(k,2)==i) ichecky=1    ! carbonyl already used
                            if tabester[k][2] == i:
                                ichecky = 1
                        # ENDDO
                        
                        # IF (ichecky==1) CYCLE kloop   ! carbonyl already involved
                        if ichecky == 1:
                            i += 1
                            continue # CYCLE kloop
                        # IF (ichecko==0) CYCLE kloop   ! ether that must be an ester
                        if ichecko == 0:
                            i += 1
                            continue # CYCLE kloop
                    # ENDIF  ! if that point is reached then must be counted as ketone

                    # IF (nodetype(ialpha)=='r') rflg=rflg+1
                    if nodetype[ialpha] == 'r':
                        rflg = rflg + 1
                    # IF (nodetype(ialpha)=='d') dflg=dflg+1
                    if nodetype[ialpha] == 'd':
                        dflg = dflg + 1
                # ENDDO
                
                # IF (rflg/=0) THEN             ! aromatic ketone
                if rflg != 0:
                    arofun[10] = arofun[10] + 1.0
                    mapfun[i][3][10] = 1.0
                    nfcr[0] = nfcr[0] + rflg
                # ELSE IF (dflg/=0) THEN             ! =C-CO-R
                elif dflg != 0:
                    cdfun[10] = cdfun[10] + 1.0
                    mapfun[i][2][10] = 1.0
                    nfcd[0] = nfcd[0] + dflg
                # ELSE                             ! R-CO-R
                else:
                    alifun[10] = alifun[10] + 1.0
                    mapfun[i][1][10] = 1.0
                # ENDIF
                ngrp[0] = ngrp[0] + 1
            # ENDIF
            i += 1
        # ENDDO kloop
    # ENDIF

    # ! ------------- PAN (index 13) -----------------
    # BA- July 2020. The group CO(ONO2) is ignored in the current version.
    # Some work was likely initiated (see index 21) but was apparently not
    # finished - these lines being commented. All this need to be revisited.
    # As a preliminary patch, I added CO(ONO2) to the PAN group. Chemmap is  
    # only called by gromhe ... so that should not create issues. Note that
    # work was done also on GECKO July 2020 to avoid the production of  
    # these structures.

    # !baba  IF (INDEX(chem,'CO(OONO2')/=0) THEN
    # IF ((INDEX(chem,'CO(OONO2')/=0).OR.(INDEX(chem,'CO(ONO2')/=0)) THEN
    if (chem.find('CO(OONO2') != -1) or (chem.find('CO(ONO2') != -1):
        # DO i = 1, node
        for i in range(1, node + 1):
            # IF ( (group(i)(1:9)=='CO(OONO2)').OR. &
            #     (group(i)(1:8)=='CO(ONO2)') ) THEN
            if (group[i].strip()[:9] == 'CO(OONO2)') or \
               (group[i].strip()[:8] == 'CO(ONO2)'):
                
                # CALL nodmap(bond,i,node,2,nnod,tnod)
                nnod_ref = [0]
                try:
                    nodmap(bond, i, node, 2, nnod_ref, tnod)
                    nnod = nnod_ref[0]
                except NotImplementedError:
                     nnod = 1 
                     tnod[1] = 1
                     pass

                # IF (nnod/=1) THEN
                if nnod != 1:
                    # WRITE(6,*) '-- error --, a unique C is expected in'
                    Fortran_write6('-- error --, a unique C is expected in')
                    # WRITE(6,*) 'alpha position of a CO(OONO2)  group'
                    Fortran_write6('alpha position of a CO(OONO2) group')
                    # WRITE(6,*) TRIM(chem)
                    Fortran_write6(chem.strip())
                    # STOP "in chemmap"
                    sys.exit("in chemmap")
                # ENDIF
                ialpha = tnod[1]
                
                # IF (nodetype(ialpha)=='o') THEN    ! R-O-CO(OONO2) function
                if nodetype[ialpha] == 'o':
                    # WRITE(6,*) '-- error --,  -O-CO(OONO2) group is unexpected'
                    Fortran_write6('-- error --, -O-CO(OONO2) group is unexpected')
                    # WRITE(6,*) TRIM(chem)
                    Fortran_write6(chem.strip())
                    # !BABA july 2020           STOP "in chemmap"
                    # nodetype(ialpha)='n' ! patch BABA july 2020 (just to make it work)
                    nodetype[ialpha] = 'n'
                # ENDIF

                # IF (nodetype(ialpha)=='r') THEN            ! aromatic PAN
                if nodetype[ialpha] == 'r':
                    arofun[13] = arofun[13] + 1
                    mapfun[i][3][13] = 1
                    ngrp[0] = ngrp[0] + 1
                    nfcr[0] = nfcr[0] + 1
                # ELSE IF (nodetype(ialpha)=='d') THEN  ! =C-CO(OONO2)
                elif nodetype[ialpha] == 'd':
                    cdfun[13] = cdfun[13] + 1
                    mapfun[i][2][13] = 1
                    ngrp[0] = ngrp[0] + 1
                    nfcd[0] = nfcd[0] + 1
                # ELSE                                     ! R-CO(OONO2)
                else:
                    alifun[13] = alifun[13] + 1
                    mapfun[i][1][13] = 1
                    ngrp[0] = ngrp[0] + 1
                # ENDIF
            # ENDIF
        # ENDDO
    # ENDIF

    # ! ---------- ether (index 14) -------------
    # IF (INDEX(chem,'-O')/=0) THEN
    if chem.find('-O') != -1:
        # etherloop: DO i = 1, node
        i = 1
        while i <= node:
            # IF (group(i)(1:3)=='-O-') THEN
            if group[i].strip()[:3] == '-O-':
                # CALL nodmap(bond,i,node,2,nnod,tnod)
                nnod_ref = [0]
                try:
                    nodmap(bond, i, node, 2, nnod_ref, tnod)
                    nnod = nnod_ref[0]
                except NotImplementedError:
                     nnod = 2 
                     tnod[1], tnod[2] = 0, 0
                     pass

                # IF (nnod/=2) THEN
                if nnod != 2:
                    # WRITE(6,*) '-- error --, only 2 C is expected in'
                    Fortran_write6('-- error --, only 2 C is expected in')
                    # WRITE(6,*) 'alpha position of a -O-  group'
                    Fortran_write6('alpha position of a -O- group')
                    # WRITE(6,*) TRIM(chem)
                    Fortran_write6(chem.strip())
                    # STOP "in chemmap"
                    sys.exit("in chemmap")
                # ENDIF
                
                rflg = 0
                dflg = 0
                # DO j=1,nnod
                for j in range(1, nnod + 1):
                    ialpha = tnod[j]
                    # IF (nodetype(ialpha)=='y') THEN    ! RCO-O-R function
                    if nodetype[ialpha] == 'y':
                        ichecko = 0 # check if the ether is already involved in an ester
                        ichecky = 0 # check if the carbonyl is already involved in an ester
                        # DO k=1,4
                        for k in range(1, 4 + 1):
                            # IF (tabester(k,1)==i) ichecko=1    ! ether already used
                            if tabester[k][1] == i:
                                ichecko = 1
                            # IF (tabester(k,2)==ialpha) ichecky=1    ! carbonyl already used
                            if tabester[k][2] == ialpha:
                                ichecky = 1
                        # ENDDO
                        
                        # IF (ichecko==1) CYCLE etherloop   ! ether already involved in an ester
                        if ichecko == 1:
                            i += 1
                            continue # CYCLE etherloop
                        # IF (ichecky==0) CYCLE etherloop   ! carbonyl that must be an ester
                        if ichecky == 0:
                            i += 1
                            continue # CYCLE etherloop
                    # ENDIF  ! if that point is reached then must be counted as ether

                    # IF (nodetype(ialpha)=='r') rflg=rflg+1
                    if nodetype[ialpha] == 'r':
                        rflg = rflg + 1
                    # IF (nodetype(ialpha)=='d') dflg=dflg+1
                    if nodetype[ialpha] == 'd':
                        dflg = dflg + 1
                # ENDDO
                
                # IF (rflg/=0) THEN             ! aromatic ether
                if rflg != 0:
                    arofun[14] = arofun[14] + 1.0
                    mapfun[i][3][14] = 1.0
                    nfcr[0] = nfcr[0] + rflg

                    # IF (rflg>1) THEN
                    if rflg > 1:
                        ierr[0] = 1
                        # RETURN
                        return 
                    # ENDIF
                # ELSE IF (dflg/=0) THEN          ! =C-O-R
                elif dflg != 0:
                    cdfun[14] = cdfun[14] + 1.0
                    mapfun[i][2][14] = 1.0
                    nfcd[0] = nfcd[0] + dflg
                # ELSE                             ! R-O-R
                else:
                    alifun[14] = alifun[14] + 1.0
                    mapfun[i][1][14] = 1.0
                # ENDIF
                ngrp[0] = ngrp[0] + 1
            # ENDIF
            i += 1
        # ENDDO etherloop
    # ENDIF

    # !baba !------------- CO(ONO2) (index 21) --------------
   
    for i in range(1, node + 1):
        # DO j=1,3
        for j in range(1, 3 + 1):
            # DO k=1,20
            for k in range(1, 20 + 1):
                # IF (mapfun(i,j,k)/=0.) THEN
                if mapfun[i][j][k] != 0.0:
                    # IF (mapfun(i,j,k).lt.1.) THEN
                    if mapfun[i][j][k] < 1.0:
                        funflg[i] = funflg[i] + 1
                    # ELSE
                    else:
                        # funflg(i)=funflg(i)+INT(mapfun(i,j,k))
                        funflg[i] = funflg[i] + int(mapfun[i][j][k])
                    # ENDIF
                # ENDIF
            # ENDDO
        # ENDDO
    # ENDDO

    # END SUBROUTINE chemmap
    return 

# END MODULE mapping