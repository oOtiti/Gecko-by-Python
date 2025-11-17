#---------------------------------------------------------------------
# Add ring-joining numerical characters to group strings     
#          at nodes given by rjg.                                     
#---------------------------------------------------------------------
import keyparameter
def rjgadd(nring,group,rjg) :
    # nring     ! # of rings in chem 
    # rjg(:,:)  ! group # of ring-joining pairs for each ring
    # group(:)  ! groups at position (node) i
# loop in reverse order => chars in numerical order if >1 exist at any node
    j=0
    for n in range(nring,0,-1) :
        for ii in range(1,3) :
            i=rjg[n-1][ii-1]
            if i!=0 :
                if(group[i-1][0:2] =='-O') or (group[i-1][1]=='d'):
                    j=3
                else :
                    j=2
                temp=group[i-1][:j-1] + keyparameter.digit[n-1] + group[i-1][j-1:]
                group[i-1]=temp

#-------------------------------------------------------------
# Strip ring-joining numerical characters from group strings. 
# Subroutine rjgadd puts them back.                           
#-------------------------------------------------------------
def rjgrm(nring,group,rjg) :
    # nring      ! # of rings in chem
    # group(:)   ! groups at position (node) i
    #I rjg(:,:)  ! group # of ring-joining pairs for each ring
    j=0
    for i in range(len(rjg)):
        for j in range(len(rjg[i])):
            rjg[i][j]=0
    for n in range(1,nring+1) :
        for i in range(1,len(group)+1) :
            if (group[i-1][0:2]=='-O') or (group[i-1][1]=='d'):
                j=3
            else :
                j=2
            if group[i-1][j-1]==keyparameter.digit[n-1]:
                temp=group[i-1][:(j-1)] + group[i-1][j:]
                group[i-1]=temp
                if rjg[n-1][0]==0 :
                    rjg[n-1][0]=i
                else :
                    rjg[n-1][1]=i

#-----------------------------------------------------------------
# PURPOSE: Add ring-joining numerical characters to group strings 
#          at nodes given by rj.                                  
#-----------------------------------------------------------------
def rjsadd(nring,chem,rjs) :
    # nring
    # rjs(:,:)  ! character # of ring-joining pairs for each ring
    # chem  ! chemical formula
# loop in reverse order => chars in numerical order if >1 exist at any node
    for n in range(nring, 0, -1):
        for j in range(2, 0, -1):
            i = rjs[n-1][j-1]  
            if i == 0:
                continue
            chem = chem[:i-1] + keyparameter.digit[n-1] + chem[i-1:]
    return chem

#----------------------------------------------------------------------
# PURPOSE: Strip ring-joining numerical characters from chem strings.  
#          Subroutine rjsadd puts them back.                           
#----------------------------------------------------------------------

def rjsrm(nring, chem, rjs):
    # nring 
    # chem ! chemical name 
    # rjs(:,:)  ! character # of ring-joining pairs for each ring
    for i in range(len(rjs)):
        for j in range(len(rjs[i])):
            rjs[i][j] = 0
    
    space_pos = chem.find(' ')
    lenchem = space_pos if space_pos != -1 else -1 
    
    for n in range(1, nring + 1):
        for j in range(1, 3):
            ptr = 1  
            for i in range(1, lenchem + 1):
                if chem[i-1] in ('C', 'c', '-'):
                    ptr = i + 1 
                    if ptr <= len(chem) and chem[ptr-1] == 'd':
                        ptr += 1
                    if ptr >= 2 and ptr <= len(chem) and chem[ptr-2:ptr] == '-O':
                        ptr += 1

                    if ptr <= len(chem) and chem[ptr-1] == keyparameter.digit[n-1]:
                        rjs[n-1][j-1] = ptr
                        chem = chem[:ptr-1] + chem[ptr:]
                        lenchem -= 1
                        break  
    return chem  