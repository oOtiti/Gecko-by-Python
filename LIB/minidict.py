from keyparameter import mxlfo, mxlco  # string length in the dictionary
import sys
mxminid = 3000  # INTEGER, PARAMETER :: mxminid=3000
nmini = 0       # INTEGER :: nmini
nam_minid = [' ' * mxlco for _ in range(mxminid)] 
fo_minid = [' ' * mxlfo for _ in range(mxminid)] 
# ======================================================================
def get_fo(nam, formula):
    formula = ' ' * mxlfo
    for i in range(1, nmini + 1):
        if nam_minid[i - 1] == nam:
            formula = fo_minid[i - 1]
            return formula 
    formula = nam
    return formula
# ======================================================================
def add_fo(nam, formula):
    global nmini, nam_minid, fo_minid
    toadd = True  
    for i in range(1, nmini + 1):
        if formula == fo_minid[i - 1]:
            toadd = False
            break
    
    if toadd:
        nmini = nmini + 1
        if nmini > mxminid:
            print("maximum formula in minid reached - see module minidict")
            sys.exit("in minid")
        fo_minid[nmini - 1] = formula
        nam_minid[nmini - 1] = nam

# ======================================================================
def clean_minid():
    global mxminid,nmini,nam_minid,fo_minid
    nmini=0
    nam_minid = [' ' * mxlco for _ in range(mxminid)]
    fo_minid = [' ' * mxlfo for _ in range(mxminid)]