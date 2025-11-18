import keyparameter
import sys
#=======================================================================
# PURPOSE: increment and add a reference to the reference list provided
# as input.
#=======================================================================
def addref(progname,reference,nref,reflist,chem=None) :
    species=''
    i=0
    # check if chem is provided
    if chem is not None:
        species=chem
    else :
        species=' '*keyparameter.mxlfo
    for i in range(1,nref+1) :
        if(reflist[i-1]==reference) :
            return
    nref=nref+1

    if nref > len(reflist):
        print("--error-- too many references added for the list")
        print("          in subroutine:", progname.rstrip()) 
        print("size of reflist is: ", len(reflist)) 
        if species[0] != ' ':  
            print("          for the species:", species.rstrip())
        print("          reference was:", reference.rstrip())
        print(reflist)
        sys.exit("in addref")
    reflist[nref-1]=reference

#=======================================================================
# PURPOSE: increment number of reactions.
#=======================================================================
def addrx(progname,chem,nr,flag) :
    mesg=' '*70
    nr=nr+1
    if nr >len(flag) :
        mesg = '-error- too many reactions created for species: ' + chem.rstrip() + ' in addrx'
        stoperr(progname,mesg,chem)
    flag[nr-1]=1

#=======================================================================
# PURPOSE: increment number of products in a reaction.
#=======================================================================
def add1tonp(progname,chem,np) :
    np=np+1
    if np>keyparameter.mxpd :
        print("--error-- too many products created, np > mxpd")
        print("          for species: ",chem.rstrip())
        print("          in subroutine: ",progname.rstrip())
        sys.exit("in add1tonp")

#=======================================================================
# PURPOSE: generates error message output
#=======================================================================
def stoperr(prog,mesg,chem) :
    kill=True
    print(' ')
    print("--error-- in: " + prog.rstrip())
    print(mesg.rstrip()) 
    print("for species: " + chem.rstrip()) 
    if kill:
        sys.exit("in stoperr") 

#=======================================================================
# PURPOSE: set bond value
#=======================================================================
def setbond(bondtb,x,y,bondvalue) :
    bondtb[x-1][y-1]=bondvalue
    bondtb[y-1][x-1]=bondvalue

#===========================================================
# PURPOSE: remove empty nodes in group and bond
#===========================================================
def erase_blank(bond, group):
    icheck = 0
    i = 0
    j = 0
    k = 0
    last = 0
    n = 0
    nca = 0
    
    tbond = [row.copy() for row in bond]
    elem_len = len(group[0]) if group else 0
    tgroup = [elem.ljust(elem_len) for elem in group]
    
    size_group = len(tgroup)
    for n in range(size_group):
        if len(tgroup[n]) > 0 and tgroup[n][0] != ' ':
            nca += 1
            last = n
    
    if nca < 2:
        return
    if last == nca - 1:
        return

    icheck = 0
    i = 0
    
    while True:
        if i >= nca:
            break
        
        if tgroup[i].strip() == '':
            for j in range(i, last):
                tgroup[j] = tgroup[j+1].ljust(elem_len)
            tgroup[last] = ' ' * elem_len
        
            for j in range(len(tbond[0])):
                for k in range(i, last):
                    tbond[k][j] = tbond[k+1][j]
                tbond[last][j] = 0
        
            for j in range(len(tbond)):
                for k in range(i, last):
                    tbond[j][k] = tbond[j][k+1]
                tbond[j][last] = 0
    
        if tgroup[i].strip() != '':
            i += 1
    
        icheck += 1
        if icheck > len(bond) * len(bond[0]):
            print("--error-- in erase_blank")
            print("infinite loop when erasing blanks lines")
            sys.exit("in erase_blank")
    
    for i in range(len(bond)):
        bond[i] = tbond[i].copy()
    for i in range(len(group)):
        group[i] = tgroup[i]

# ======================================================================
# PURPOSE: count the number of occurence of a substring in a string.
# ======================================================================
def countstring(line, str_sub):
    n = 0
    if len(str_sub) == 0:
        return n
    p = 0
    while True:
        ipos = line.find(str_sub, p)
        if ipos == -1:
            return n
        n += 1
        p = ipos + len(str_sub)

def kval(arrh, T):
    import math
    return arrh[0] * (T ** arrh[1]) * math.exp(-arrh[2] / T)