# MODULE dictstacktool
from toolbox import stoperr, addref
from keyparameter import mxlfl, mxlco, mxlfo, mxpd, mecu
from references import mxlcod
from keyflag import enolflg, isomerfg
from cdtool import switchenol
from atomtool import cnum, onum
from searching import srch, srh5
from database import nspsp, dictsp
from dictstackdb import (nrec, ninorg, inorglst, dict, namlst, dbrch, diccri,
                         stabl, lotopstack, nhldvoc, holdvoc, nhldrad, holdrad,
                         level)
from namingtool import naming
from switchisomer import isomer
from tempflag import iflost
from normchem import stdchm
from rxwrttool import rxwrit, rxinit

# CONTAINS
# ! SUBROUTINE bratio(pchem,brtio,pname)
# ! SUBROUTINE add_topvocstack(inchem,brtio,pname,ncom,com)
# ! SUBROUTINE loader(chem,name)
# ! SUBROUTINE scrollc1stack(c1stack)
# ! SUBROUTINE addc1chem(pchem,pname,nstck,c1stack,ierr)
# ! SUBROUTINE addc1dict(pchem,pname,fgrp)

# !=======================================================================
# ! PURPOSE : gives the short name of the species (pchem) given as input 
# ! and updates the stack and dictionary arrays if necessary.
# !  If the species is already known, the routine just returns its short 
# !  name. If the species is new, then:                                   
# !    1- a short name is given to the species (see naming)          
# !    2- the "dictionary" arrays are updated (dict, namlst, dbrch).            
# !       Species is added in such a way that the tables remain sorted   
# !    3- the species is added to stack for future reactions (see loader)
# !=======================================================================
def bratio(inchem, brtio, pname, ncom, com):
    from toolbox import stoperr, addref
    from keyparameter import mxlfl, mxlco, mxlfo, mxpd, mecu
    from references import mxlcod
    from keyflag import enolflg, isomerfg
    from cdtool import switchenol
    from atomtool import cnum, onum
    from searching import srch, srh5
    from database import nspsp, dictsp
    from dictstackdb import (nrec, ninorg, inorglst, dict, namlst, dbrch, diccri,
                            stabl, lotopstack, nhldvoc, holdvoc, nhldrad, holdrad,
                            level)
    from namingtool import naming
    from switchisomer import isomer
    from tempflag import iflost
    from normchem import stdchm
    from rxwrttool import rxwrit, rxinit
    # 初始化变量
    progname = 'bratio'
    mesg = ''
    fgrp = ''
    pchem = inchem  # make a working copy of inchem
    dicptr = 0
    namptr = 0
    nca = 0
    chg = 0
    i = 0
    ipos = 0
    ierr = 0
    tabinfo = [0] * len(diccri[0]) if diccri else []  # SIZE(diccri,2)
    loswitch = False  # switch for keto/enol change
    
    mx1c = 15  # size of the C1 stack
    c1stack = [''] * mx1c  # C1 stack 
    nstck = 0  # # of element in C1 stack 
    
    if inchem == ' ':
        return  # return if no species 

    # check keyword (can be given by special mechanisms)
    if inchem[:5] == 'EXTRA':
        pname = 'EXTRA '
        return
    elif inchem[:2] == 'HV':
        pname = 'HV    '
        return
    elif inchem[:6] == 'OXYGEN':
        pname = 'OXYGEN'
        return
    elif inchem[:6] == 'ISOM  ':
        pname = 'ISOM  '
        return
    elif inchem[:4] == '(+M)':
        pname = '(+M)  '
        return  # !! "FALLOFF"=7
    
    pname = ' '
    nca = cnum(pchem) + onum(pchem)
    tabinfo = [0] * len(tabinfo)

    # special name (inorganics, formulae that cannot be held, and C1)
    # =============================================================

    # no carbon - check if known in the list of inorganics or keywords. 
    if nca == 0:
        for i in range(1, ninorg + 1):
            # inorglst(i)(10:mxlfo+11) -> Python: inorglst[i-1][9:mxlfo+11]
            # pchem(1:LEN_TRIM(pchem)+1) -> pchem[:len(pchem.strip())+1]
            substr_inorg = inorglst[i-1][9 : mxlfo+11]
            substr_pchem = pchem[: len(pchem.strip()) + 1]
            if substr_pchem in substr_inorg and substr_inorg.index(substr_pchem) == 0:
                pname = inorglst[i-1][:mxlco]
                return

    # catcher for "COO"
    if inchem[:4] == 'COO ':
        pchem = 'CO2 ' + pchem[4:]

    # Check if "special species" (formula starts with '#')
    # =============================================================
    if pchem[0] == '#':
        ipos = 0
        dicptr = srch(nrec, pchem, dict)

        # if already recorded ...
        if dicptr > 0:
            pname = dict[dicptr-1][:mxlco]  # dict(dicptr) -> dict[dicptr-1]
            dbrch[dicptr-1] = max(dbrch[dicptr-1], brtio)
            return

        # If not recorded, then search in the special list 
        ipos = 0
        for i in range(1, nspsp + 1):
            # dictsp(i)(10:mxlfo+11) -> dictsp[i-1][9:mxlfo+11]
            # pchem(1:INDEX(pchem," ")) -> pchem[:pchem.index(' ')]
            space_idx = pchem.index(' ') if ' ' in pchem else len(pchem)
            substr_pchem = pchem[:space_idx]
            substr_dictsp = dictsp[i-1][9 : mxlfo+11]
            if substr_pchem in substr_dictsp and substr_dictsp.index(substr_pchem) == 0:
                ipos = i
                pname = dictsp[i-1][:mxlco]
                namptr = srh5(pname, namlst, nrec)
                namptr = -namptr
                # if species is already in the namlist, we don't need to add it.
                if namptr < 0:
                    return
                fgrp = dictsp[i-1][mxlfo+12:]
                break

        if ipos == 0:
            mesg = "Special species cannot be managed. Check formula."
            stoperr(progname, mesg, pchem)

    # C1 species, must already be known !
    # =============================================================
    elif nca == 1:
        dicptr = srch(nrec, pchem, dict)
        if dicptr > 0:  # species already included
            pname = dict[dicptr-1][:mxlco]
            return

        # not yet known ...
        elif iflost == 0:
            # check if species & chemistry is "available" in C1 routines (add all then & return pname) 
            ierr = 1
            nstck = 0
            c1stack = [''] * mx1c
            addc1chem(pchem, pname, nstck, c1stack, ierr)

            # species unknown - stop
            if ierr != 0:
                mesg = "Following C1 species not in the dictionary"
                stoperr(progname, mesg, pchem)

            # If additional species was produced by pchem, then add also those species      
            if nstck != 0:
                scrollc1stack(c1stack)

            return  # all clear - return with pname

    # Regular species 
    # =============================================================
    else:
        # Regular formula and must start with a "C,c, or O"  
        if (pchem[0] != 'C') and (pchem[:2] != '-O') and (pchem[0] != 'c'):
            mesg = "Species cannot be managed. Check formula."
            stoperr(progname, mesg, pchem)

        # check enols and switch to keto form
        if enolflg:
            if 'Cd' in pchem:
                switchenol(pchem, loswitch)  # pchem is returned as std keto if enol
                if loswitch:
                    addref(progname, 'KETOENOL', ncom, com, inchem)

        # Search if pchem is already recorded. If yes, return the short name
        dicptr = srch(nrec, pchem, dict)
        if dicptr > 0:
            pname = dict[dicptr-1][:mxlco]
            dbrch[dicptr-1] = max(dbrch[dicptr-1], brtio)
            return

        # the species is unknown 
        # --------------------------------------

        # LUMP SECONDARY SPECIES - not available. See old version of geckoa 
        # including subroutine "lump_sec" and the corresponding bratio version. 

        # ISOMER SWITCH - search if an isomer is already known and replace formula 
        if isomerfg and nca > 3 and '.' not in pchem:
            isomer(pchem, brtio, chg, tabinfo)
            if chg == 1:  # pchem was switched - get corresponding name
                dicptr = srch(nrec, pchem, dict)
                if dicptr <= 0:
                    mesg = "expected species after isomer switch not found"
                    stoperr(progname, mesg, pchem)
                dbrch[dicptr-1] = max(dbrch[dicptr-1], brtio)
                pname = dict[dicptr-1][:mxlco]
                addref(progname, 'SWAPISOMER', ncom, com, inchem)
                return 

        # Get name (pname) for pchem and position (namptr) in namlst for addition
        naming(pchem, namptr, pname, fgrp)

    # update stack and dictionary array (species may come be a '#' species)
    # ====================================================================

    # if the flag to stop the chemistry is raised, then return
    if iflost == 1:
        return

    # raise the counters
    nrec += 1
    if nrec > len(dict):
        mesg = "Too many species added to dictionary"
        stoperr(progname, mesg, inchem)

    # add new name and raise the name array
    namptr += 1
    # namlst(namptr+1:nrec+1)=namlst(namptr:nrec) -> Python索引偏移-1
    namlst[namptr : nrec+1] = namlst[namptr-1 : nrec]
    namlst[namptr-1] = pname  # namlst(namptr) -> namlst[namptr-1]

    # raise dictionary and branching table. Add new line for the new species 
    dicptr = abs(dicptr) + 1
    # dict(dicptr+1:nrec+1)=dict(dicptr:nrec)
    dict[dicptr : nrec+1] = dict[dicptr-1 : nrec]
    # dbrch(dicptr+1:nrec+1)=dbrch(dicptr:nrec)
    dbrch[dicptr : nrec+1] = dbrch[dicptr-1 : nrec]
    # WRITE(dict(dicptr),'(a6,3x,a120,2x,a15)') pname, pchem, fgrp
    dict[dicptr-1] = f"{pname:<6}   {pchem:<120}  {fgrp:<15}"
    dbrch[dicptr-1] = brtio

    # store information required to search for an isomer (info saved in tabinfo)
    if isomerfg:
        # diccri(dicptr+1:nrec+1,:) = diccri(dicptr:nrec,:)
        diccri[dicptr : nrec+1] = diccri[dicptr-1 : nrec]
        diccri[dicptr-1] = tabinfo.copy()

    loader(pchem, pname)

# !=======================================================================
# ! PURPOSE: raise flag for an addition on top of the VOC stack before
# ! calling bratio. If the species is already known in the mechanism, then
# ! nothing should happen: the routine simply return the short name 
# ! (pname) corresponding to the input species (inchem). If the species 
# ! is new, then it will be added on top of the VOC stack, with a # of 
# ! generation identical to the current value.
# !=======================================================================
def add_topvocstack(inchem, brtio, pname, ncom, com):
    from toolbox import stoperr, addref
    from keyparameter import mxlfl, mxlco, mxlfo, mxpd, mecu
    from references import mxlcod
    from keyflag import enolflg, isomerfg
    from cdtool import switchenol
    from atomtool import cnum, onum
    from searching import srch, srh5
    from database import nspsp, dictsp
    from dictstackdb import (nrec, ninorg, inorglst, dict, namlst, dbrch, diccri,
                            stabl, lotopstack, nhldvoc, holdvoc, nhldrad, holdrad,
                            level)
    from namingtool import naming
    from switchisomer import isomer
    from tempflag import iflost
    from normchem import stdchm
    from rxwrttool import rxwrit, rxinit
    # 初始化变量
    savstabl = stabl 

    # change data in dictstackdb
    stabl -= 1  # because loader will add 1 to stabl once called
    lotopstack = True  # raise flag for top addition
    
    # return short name and add species to the stack (if new)
    bratio(inchem, brtio, pname, ncom, com)

    # restore data in dictstackdb 
    stabl = savstabl  # curretnt value if stabl
    lotopstack = False  # default value for stack management

# !=======================================================================
# ! PURPOSE: Load molecules in the stack for further reactions. 
# ! Two stacks are considered, one for the non radical VOC and one for 
# ! the radicals. Reactions for radical are written first (see main.f).
# ! New species are added to the bottom of the stack (first in first out)
# ! as the default situation.
# !                                                                     
# ! Information loaded in the stack comprises (holdvoc & holdrad):
# !   holdvoc(1:6) = short name of the chemical                       
# !   holdvoc(7:126) = chemical formula                               
# !   holdvoc(127:129) = # of stable generations
# !   holdvoc(130:132) = # of levels (including radicals)
# !                                                                     
# ! INPUT from module dictstackdb:
# !  - level      : # of level (stable+radical) needed to produce "chem"
# !  - stabl      : # of generation (stable) needed to produce "chem"
# !  - lotopstack : if true, then VOC is added on top on the stack (lifo)  
# ! INPUT/OUTPUT from module dictstackdb:
# !  - nhldvoc    : number of (stable) VOC in the stack 
# !  - holdvoc(i) : list of the VOC in the stack
# !  - nhldrad    : number of radical in the stack
# !  - holdrad(i) : list of the radicals in the stack
# !=======================================================================
def loader(chem, idnam):
    from toolbox import stoperr, addref
    from keyparameter import mxlfl, mxlco, mxlfo, mxpd, mecu
    from references import mxlcod
    from keyflag import enolflg, isomerfg
    from cdtool import switchenol
    from atomtool import cnum, onum
    from searching import srch, srh5
    from database import nspsp, dictsp
    from dictstackdb import (nrec, ninorg, inorglst, dict, namlst, dbrch, diccri,
                            stabl, lotopstack, nhldvoc, holdvoc, nhldrad, holdrad,
                            level)
    from namingtool import naming
    from switchisomer import isomer
    from tempflag import iflost
    from normchem import stdchm
    from rxwrttool import rxwrit, rxinit
    # 初始化变量
    progname = 'loader'
    mesg = ''
    tplevel = 0
    tpstabl = 0

    # raise the counters:
    tplevel = level + 1 
    if '.' not in chem:
        tpstabl = stabl + 1
    else:
        tpstabl = stabl

    # VOC (non radical) species
    # --------------------------
    if '.' not in chem:  # usual VOC
        nhldvoc += 1
        if nhldvoc > len(holdvoc):
            mesg = "Stack full. # of VOC species > VOC stack size"
            stoperr(progname, mesg, chem)
        if lotopstack:  # special case - add species on top of the stack (lifo)
            # holdvoc(2:nhldvoc)=holdvoc(1:nhldvoc-1)
            holdvoc[1:nhldvoc] = holdvoc[:nhldvoc-1]
            # WRITE(holdvoc(1),'(a6,a120,i3,i3)') idnam,chem,tpstabl,tplevel
            holdvoc[0] = f"{idnam:<6}{chem:<120}{tpstabl:3d}{tplevel:3d}"
        else:  # default (fifo)
            # WRITE(holdvoc(nhldvoc),'(a6,a120,i3,i3)') idnam,chem,tpstabl,tplevel 
            holdvoc[nhldvoc-1] = f"{idnam:<6}{chem:<120}{tpstabl:3d}{tplevel:3d}"
      
    # radical species
    # ---------------
    else:  # usual radical
        nhldrad += 1
        if nhldrad > len(holdrad):
            mesg = "Stack full. # of radical species > radical stack size"
            stoperr(progname, mesg, chem)
        # WRITE(holdrad(nhldrad),'(a6,a120,i3,i3)') idnam,chem,tpstabl,tplevel
        holdrad[nhldrad-1] = f"{idnam:<6}{chem:<120}{tpstabl:3d}{tplevel:3d}"

# ! ======================================================================
# ! PURPOSE: Write the species in the "C1stack" to the dictionaries and 
# ! the corresponding chemistry to the mechanism. That chemistry may lead
# ! to the addition of more C1 species to the stack. The routine scrolls
# ! all species in the stack and stops when the stack is empty.
# !
# ! For future development: all this C1 chemistry/species is "hard coded".
# ! Some future development might be useful to make it "soft" (reading
# ! input files)
# ! ======================================================================
def scrollc1stack(c1stack):
    # 初始化变量
    progname = 'scrollc1stack'
    mesg = ''
    chem = ''
    icount = 0
    nstck = 0
    ierr = 0
    dicptr = 0
    dumnam = ''

    nstck = sum(1 for x in c1stack if x != '')
    icount = 0

    # loop over C1 species in stack
    while True:  # just in case, count the ...   
        icount += 1  # ... number of iteration to ...
        if icount == 100:  # ... infinite loop.
            mesg = "infinite loop when adding C1 species"
            stoperr(progname, mesg, c1stack[0])

        # get next species from stack (exit if empty)
        chem = c1stack[0]
        # PRINT*, 'next from stack:',chem 
        if chem == ' ':
            break  # no more species in stack, exit

        # check if the species already added to dictionary 
        dicptr = srch(nrec, chem, dict)
        # PRINT*, 'dictpr=',dicptr
        if dicptr > 0:
            if nstck == 1:
                break  # exit - nothing left in stack 
            # c1stack(1:nstck-1)=c1stack(2:nstck)   ! move -1 stack list
            c1stack[:nstck-1] = c1stack[1:nstck]
            nstck -= 1

        # add new species
        else:
            ierr = 1
            addc1chem(chem, dumnam, nstck, c1stack, ierr)
            if ierr != 0:  # species not available
                mesg = "Following C1 species not in the dictionary"
                stoperr(progname, mesg, c1stack[0])

            if nstck == 1:
                break  # exit - nothing left in stack 
            # c1stack(1:nstck-1)=c1stack(2:nstck)   ! move -1 stack list
            c1stack[:nstck-1] = c1stack[1:nstck]
            nstck -= 1

# ! ======================================================================
# ! PURPOSE: check if the formula provided as input (pchem) is "known" as
# ! a "regular" species. If so, then the species is added in the
# ! dictionaries and the chemistry (if any) is added in the mechanism.
# ! The subroutine returns the "regular" short name (pname) for the input
# ! species. 
# ! New C1 species (provided from the chemistry of pchem) that may not
# ! already be considered are added to a C1 stack for further addition in
# ! dictionaries and mechanisms.
# ! NB: if new C1 species are added here, they MUST ALSO be added to
# !     array "optional_c1" in spsptool.f90
# ! For future development: all this C1 chemistry/species is "hard coded".
# ! Some further development might be useful to make it "soft" (reading
# ! input files).
# ! The list of allowed formulae and corresponding names is:
# ! CH2OO  : CH2.(OO.)
# ! CH3O3H : CH3(OOOH)
# ! COHOOH : CH2(OH)(OOH)
# ! COHO3H : CH2(OH)(OOOH)
# ! CHONO2 : CHO(NO2)
# ! C10001 : CH2(ONO2)(OO.)
# ! C10002 : CO(NO2)(ONO2)
# ! C10003 : CO(OH)(NO2)
# ! C10004 : CO(OH)(ONO2)
# ! C10005 : CO(OH)(OOH)
# ! C10006 : CO(ONO2)(ONO2)
# ! C10007 : CO(ONO2)(OOH)
# ! C10008 : CO(OOH)(OO.)
# ! C10009 : CO(NO2)(OOOH)
# ! C10010 : CO(OOOH)(OOOH)
# ! H2CO3  : CO(OH)(OH)
# ! HCOO3H : CHO(OOOH)                   
# ! NO1001 : CH2(OH)(ONO2)
# ! NH1001 : CH2(ONO2)(OOH)
# ! 2V1001 : CH2(NO2)(OO.)
# ! C1ND00 : CHO(ONO2)
# ! N01003 : CH3(ONO)        ! still in dic_c1.dat and mch_singlec.dat
# !
# ! to be continued ...
# ! ======================================================================
def addc1chem(pchem, pname, nstck, c1stack, ierr):
    from toolbox import add1tonp

    # 初始化变量
    fgrp = ''
    tchem = ''
    np = 0
    namptr = 0
    p = [''] * mxpd
    r = [''] * 3
    s = [0.0] * mxpd
    arrh = [0.0] * 3
    idreac = 0
    nlabel = 0
    xlabel = 0.0
    folow = [0.0] * 3
    fotroe = [0.0] * 4
    brtio = 0.0
    namacetic = ' '  # intialize at 1st call only

    mxcom = 10
    nref = 0
    ref = [''] * mxcom
    progname = 'addc1chem'
    mesg = ''

    # rate constants for CH2.(OO). criegee radical
    kuni_ch2oo = [1.66E+01, 4.02, 8024.]
    kmh2o_ch2oo = [8.13E-18, 1.23, 698.]
    kdh2o_ch2oo = [7.95E-18, 1.24, -1510.]
    kso2_ch2oo = [3.70E-11, 0., 0.]
    kno2_ch2oo = [3.00E-12, 0., 0.]
    khno3_ch2oo = [5.40E-10, 0., 0.]
    krcooh_ch2oo = [1.20E-10, 0., 0.]
    
    pname = ' '
    ierr = 1  # raise flag as default 

    # ----------
    # CH2.(OO).
    # ----------
    if pchem == 'CH2.(OO.)':
        # add the species to the dictionary
        pname = 'CH2OO'
        fgrp = '4 '
        addc1dict(pchem, pname, fgrp)
        ierr = 0  # flag down, species found

        # add chemistry

        # unimolecular decomposition (CO2+H2, CO+H2O, 50% each)
        np = 0
        rxinit(r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)
        nref = 0
        ref = [''] * mxcom
        addref(progname, 'MN22SAR2', nref, ref, pchem)

        add1tonp(progname, pchem, np)
        s[np-1] = 0.5
        p[np-1] = 'CO2  '
        add1tonp(progname, pchem, np)
        s[np-1] = 0.5
        p[np-1] = 'H2   '
        add1tonp(progname, pchem, np)
        s[np-1] = 0.5
        p[np-1] = 'CO   '
        add1tonp(progname, pchem, np)
        s[np-1] = 0.5
        p[np-1] = 'H2O  '

        r[0] = pname
        arrh = kuni_ch2oo.copy()
        rxwrit(mecu, r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe, ref)

        # reaction with water
        np = 0
        rxinit(r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)
        nref = 0
        ref = [''] * mxcom
        addref(progname, 'MN22SAR2', nref, ref, pchem)

        add1tonp(progname, pchem, np)
        s[np-1] = 0.55
        p[np-1] = 'COHOOH'
        add1tonp(progname, pchem, np)
        s[np-1] = 0.40
        p[np-1] = 'CH2O  '
        add1tonp(progname, pchem, np)
        s[np-1] = 0.40
        p[np-1] = 'H2O2  '
        add1tonp(progname, pchem, np)
        s[np-1] = 0.05
        p[np-1] = 'HCOOH '

        r[0] = pname
        r[1] = 'EXTRA'
        idreac = 2
        nlabel = 500
        arrh = kmh2o_ch2oo.copy()
        rxwrit(mecu, r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe, ref)

        # reaction with water dimer
        np = 0
        rxinit(r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)
        nref = 0
        ref = [''] * mxcom
        addref(progname, 'MN22SAR2', nref, ref, pchem)

        add1tonp(progname, pchem, np)
        s[np-1] = 0.55
        p[np-1] = 'COHOOH'
        add1tonp(progname, pchem, np)
        s[np-1] = 0.40
        p[np-1] = 'CH2O  '
        add1tonp(progname, pchem, np)
        s[np-1] = 0.40
        p[np-1] = 'H2O2  '
        add1tonp(progname, pchem, np)
        s[np-1] = 0.05
        p[np-1] = 'HCOOH '

        r[0] = pname
        r[1] = 'EXTRA'
        idreac = 2
        nlabel = 502
        arrh = kdh2o_ch2oo.copy()
        rxwrit(mecu, r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe, ref)

        # new C1 species that may need to be added    
        namptr = srh5('COHOOH', namlst, nrec)
        if namptr < 0:
            nstck += 1
            if nstck > len(c1stack):
                mesg = "too many species added to C1 stack"
                stoperr(progname, mesg, pchem)
            c1stack[nstck-1] = 'CH2(OH)(OOH) '

        # reaction with SO2
        np = 0
        rxinit(r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)
        nref = 0
        ref = [''] * mxcom
        addref(progname, 'MN22SAR2', nref, ref, pchem)

        add1tonp(progname, pchem, np)
        s[np-1] = 1.
        p[np-1] = 'CH2O  '
        add1tonp(progname, pchem, np)
        s[np-1] = 1.
        p[np-1] = 'SULF  '

        r[0] = pname
        r[1] = 'SO2 '
        arrh = kso2_ch2oo.copy()
        rxwrit(mecu, r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe, ref)

        # reaction with NO2 (temporary version, product expected is CH2(NO2)(OO.)
        np = 0
        rxinit(r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)
        nref = 0
        ref = [''] * mxcom
        addref(progname, 'MN22SAR2', nref, ref, pchem)

        add1tonp(progname, pchem, np)
        s[np-1] = 1.
        p[np-1] = '2V1001'

        r[0] = pname
        r[1] = 'NO2 '
        arrh = kno2_ch2oo.copy()
        rxwrit(mecu, r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe, ref)

        # new C1 species that may need to be added    
        namptr = srh5('2V1001', namlst, nrec)
        if namptr < 0:
            nstck += 1
            if nstck > len(c1stack):
                mesg = "too many species added to C1 stack"
                stoperr(progname, mesg, pchem)
            c1stack[nstck-1] = 'CH2(NO2)(OO.) '

        # reaction with HNO3 (temporary version, product expected is CH2(ONO2)(OOH)
        np = 0
        rxinit(r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)
        nref = 0
        ref = [''] * mxcom
        addref(progname, 'MN22SAR2', nref, ref, pchem)

        add1tonp(progname, pchem, np)
        s[np-1] = 1.
        p[np-1] = 'NH1001'

        r[0] = pname
        r[1] = 'HNO3 '
        arrh = khno3_ch2oo.copy()
        rxwrit(mecu, r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe, ref)

        # new C1 species that may need to be added    
        namptr = srh5('NH1001', namlst, nrec)
        if namptr < 0:
            nstck += 1
            if nstck > len(c1stack):
                mesg = "too many species added to C1 stack"
                stoperr(progname, mesg, pchem)
            c1stack[nstck-1] = 'CH2(ONO2)(OOH) '

        # reaction with HCOOH
        np = 0
        rxinit(r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)
        nref = 0
        ref = [''] * mxcom
        addref(progname, 'MN22SAR2', nref, ref, pchem)

        add1tonp(progname, pchem, np)
        s[np-1] = 1.
        tchem = 'CHO-O-CH2(OOH)'
        stdchm(tchem)
        brtio = 1.
        bratio(tchem, brtio, p[np-1], nref, ref)

        r[0] = pname
        r[1] = 'HCOOH '
        arrh = krcooh_ch2oo.copy()
        rxwrit(mecu, r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe, ref)

        # reaction with CH3COOH
        np = 0
        rxinit(r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)
        nref = 0
        ref = [''] * mxcom
        addref(progname, 'MN22SAR2', nref, ref, pchem)

        add1tonp(progname, pchem, np)
        s[np-1] = 1.
        tchem = 'CH3CO-O-CH2(OOH)'
        stdchm(tchem)
        brtio = 1.
        bratio(tchem, brtio, p[np-1], nref, ref)

        r[0] = pname
        if namacetic == ' ':
            tchem = 'CH3CO(OH) '
            brtio = 1.  # CH3COOH assumed as a primary species
            add_topvocstack(tchem, brtio, namacetic, nref, ref)
        r[1] = namacetic
        arrh = krcooh_ch2oo.copy()
        rxwrit(mecu, r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe, ref)

    # ------------
    # CH3(OOOH)
    # ------------
    elif pchem == 'CH3(OOOH)':
        # these species are produced from the RO2 + OH reaction (see Jenkin et al., 2019)
        # use rate for CH3OOOH + HO => CH3O. + H2O + O2 from Anglada and Solé (2018)
        # use decompostion rate from Assaf et al., 2018 (80% path 1, 20% path 2)

        # add the species to the dictionary
        pname = 'CH3O3H'
        fgrp = '   '
        addc1dict(pchem, pname, fgrp)
        ierr = 0  # flag down, species found

        # add chemistry
        # reaction with OH
        np = 0
        rxinit(r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)
        nref = 0
        ref = [''] * mxcom
        addref(progname, 'JA18KMV000', nref, ref, pchem)

        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'CH3O  '

        r[0] = pname
        r[1] = 'HO '
        arrh[0] = 1.46E-12
        arrh[1] = 0.
        arrh[2] = -1037
        rxwrit(mecu, r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe, ref)

        # decomposition
        np = 0
        rxinit(r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)
        nref = 0
        ref = [''] * mxcom
        addref(progname, 'EA18KMR000', nref, ref, pchem)

        add1tonp(progname, pchem, np)
        s[np-1] = 0.8
        p[np-1] = 'CH3O  '
        add1tonp(progname, pchem, np)
        s[np-1] = 0.8
        p[np-1] = 'HO2   '
        add1tonp(progname, pchem, np)
        s[np-1] = 0.2
        p[np-1] = 'CH3OH '

        r[0] = pname
        arrh[0] = 1.52E10
        arrh[1] = 1.35
        arrh[2] = 12000
        rxwrit(mecu, r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe, ref)

    # ------------
    # CHO(OOOH)
    # ------------
    elif pchem == 'CHO(OOOH)':
        # these species are produced from the RO2 + OH reaction (see Jenkin et al., 2019)
        # use rate for CH3OOOH + HO => CH3O. + H2O + O2 from Anglada and Solé (2018)
        # use decompostion rate from Assaf et al., 2018 (80% path 1, 20% path 2)

        # add the species to the dictionary
        pname = 'HCOO3H'
        fgrp = '   '
        addc1dict(pchem, pname, fgrp)
        ierr = 0  # flag down, species found

        # add chemistry
        # reaction with OH
        np = 0
        rxinit(r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)
        nref = 0
        ref = [''] * mxcom
        addref(progname, 'JA18KMV000', nref, ref, pchem)

        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'CO    '
        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'HO2   '

        r[0] = pname
        r[1] = 'HO '
        arrh[0] = 1.46E-12
        arrh[1] = 0.
        arrh[2] = -1037
        rxwrit(mecu, r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe, ref)

        # decomposition
        np = 0
        rxinit(r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)
        nref = 0
        ref = [''] * mxcom
        addref(progname, 'EA18KMR000', nref, ref, pchem)

        add1tonp(progname, pchem, np)
        s[np-1] = 0.8
        p[np-1] = 'CO    '
        add1tonp(progname, pchem, np)
        s[np-1] = 1.6
        p[np-1] = 'HO2   '
        add1tonp(progname, pchem, np)
        s[np-1] = 0.2
        p[np-1] = 'HCOOH '

        r[0] = pname
        arrh[0] = 1.52E10
        arrh[1] = 1.35
        arrh[2] = 12000
        rxwrit(mecu, r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe, ref)

    # ------------
    # CH2(OH)(OOH)
    # ------------
    elif pchem == 'CH2(OH)(OOH)':
        # add the species to the dictionary
        pname = 'COHOOH'
        fgrp = 'HO '
        add1tonp(pchem, pname, fgrp)
        ierr = 0  # flag down, species found

        # add chemistry
        # reaction with OH
        np = 0
        rxinit(r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)
        nref = 0
        ref = [''] * mxcom
        addref(progname, 'HA18KMV000', nref, ref, pchem)

        add1tonp(progname, pchem, np)
        s[np-1] = 0.55
        p[np-1] = 'CH2O  '
        add1tonp(progname, pchem, np)
        s[np-1] = 0.55
        p[np-1] = 'HO2   '
        add1tonp(progname, pchem, np)
        s[np-1] = 0.45
        p[np-1] = 'HCOOH '
        add1tonp(progname, pchem, np)
        s[np-1] = 0.45
        p[np-1] = 'HO    '

        r[0] = pname
        r[1] = 'HO '
        arrh[0] = 7.1E-12
        arrh[1] = 0.
        arrh[2] = 0.
        rxwrit(mecu, r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe, ref)

    # ------------
    # CH2(OH)(OOOH)
    # ------------
    elif pchem == 'CH2(OH)(OOOH)':
        # add the species to the dictionary
        pname = 'COHO3H'
        fgrp = 'O '
        addc1dict(pchem, pname, fgrp)
        ierr = 0  # flag down, species found

        # add chemistry
        #       ..... WRITE THE REACTIONS (INCLUDING PHASE TRANSFER) ...
        np = 0
        rxinit(r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)
        nref = 0
        ref = [''] * mxcom
        addref(progname, 'JO23PRSCOM', nref, ref, pchem)

        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'CH2O  '
        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'HO2   '

        r[0] = pname
        arrh[0] = 0.1
        arrh[1] = 0.
        arrh[2] = 0.
        r[1] = 'HV   '
        idreac = 1
        nlabel = 4
        xlabel = arrh[0]
        rxwrit(mecu, r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe, ref)

    # ------------
    # CHO(NO2)
    # ------------
    elif pchem == 'CHO(NO2)':
        # add the species to the dictionary
        pname = 'CHONO2'
        fgrp = 'DV '
        addc1dict(pchem, pname, fgrp)
        ierr = 0  # flag down, species found

        # add chemistry
        #       ..... WRITE THE REACTIONS (INCLUDING PHASE TRANSFER) ...
        np = 0
        rxinit(r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)
        nref = 0
        ref = [''] * mxcom
        addref(progname, 'JO23PRSCOM', nref, ref, pchem)

        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'CO    '
        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'HO2   '
        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'NO2   '

        r[0] = pname
        arrh[0] = 0.1
        arrh[1] = 0.
        arrh[2] = 0.
        r[1] = 'HV   '
        idreac = 1
        nlabel = 4
        xlabel = arrh[0]
        rxwrit(mecu, r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe, ref)

    # ------------
    # CH2(ONO2)(OO.)
    # ------------
    elif pchem == 'CH2(ONO2)(OO.)':
        # add the species to the dictionary
        pname = 'C10001'
        fgrp = '2N  '
        addc1dict(pchem, pname, fgrp)
        ierr = 0  # flag down, species found

        # add chemistry
        # decomposition
        np = 0
        rxinit(r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)
        nref = 0
        ref = [''] * mxcom  # ; CALL addref(progname,'XXXXXXXXXX',nref,ref,pchem)

        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'CH2O  '
        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'NO2   '

        r[0] = pname
        arrh[0] = 1.00E6
        arrh[1] = 0.
        arrh[2] = 0.
        rxwrit(mecu, r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe, ref)

    # ------------
    # CO(NO2)(ONO2)
    # ------------
    elif pchem == 'CO(NO2)(ONO2)':
        # add the species to the dictionary
        pname = 'C10002'
        fgrp = 'KNV '
        addc1dict(pchem, pname, fgrp)
        ierr = 0  # flag down, species found

        # add chemistry
        #       ..... WRITE THE REACTIONS (INCLUDING PHASE TRANSFER) ...

        np = 0
        rxinit(r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)
        nref = 0
        ref = [''] * mxcom
        addref(progname, 'JO23PRSCOM', nref, ref, pchem)

        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'CO2  '
        add1tonp(progname, pchem, np)
        s[np-1] = 2.0
        p[np-1] = 'NO2   '

        r[0] = pname
        arrh[0] = 0.1
        arrh[1] = 0.
        arrh[2] = 0.
        r[1] = 'HV   '
        idreac = 1
        nlabel = 4
        xlabel = arrh[0]
        rxwrit(mecu, r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe, ref)

    # ------------
    # CO(OH)(NO2)
    # ------------
    elif pchem == 'CO(OH)(NO2)':
        # add the species to the dictionary
        pname = 'C10003'
        fgrp = 'AV '
        addc1dict(pchem, pname, fgrp)
        ierr = 0  # flag down, species found

        # add chemistry
        #       ..... WRITE THE REACTIONS (INCLUDING PHASE TRANSFER) ...

        np = 0
        rxinit(r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)
        nref = 0
        ref = [''] * mxcom
        addref(progname, 'JO23PRSCOM', nref, ref, pchem)

        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'HO2   '
        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'CO2   '
        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'NO2   '

        r[0] = pname
        arrh[0] = 0.1
        arrh[1] = 0.
        arrh[2] = 0.
        r[1] = 'HV   '
        idreac = 1
        nlabel = 4
        xlabel = arrh[0]
        rxwrit(mecu, r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe, ref)

    # ------------
    # CO(OH)(ONO2)
    # ------------
    elif pchem == 'CO(OH)(ONO2)':
        # add the species to the dictionary
        pname = 'C10004'
        fgrp = 'AN '
        addc1dict(pchem, pname, fgrp)
        ierr = 0  # flag down, species found

        # add chemistry
        #       ..... WRITE THE REACTIONS (INCLUDING PHASE TRANSFER) ...

        np = 0
        rxinit(r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)
        nref = 0
        ref = [''] * mxcom
        addref(progname, 'JO23PRSCOM', nref, ref, pchem)

        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'HO2   '
        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'CO2   '
        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'NO3   '

        r[0] = pname
        arrh[0] = 0.1
        arrh[1] = 0.
        arrh[2] = 0.
        r[1] = 'HV   '
        idreac = 1
        nlabel = 4
        xlabel = arrh[0]
        rxwrit(mecu, r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe, ref)

    # ------------
    # CO(OH)(OOH)
    # ------------
    elif pchem == 'CO(OH)(OOH)':
        # add the species to the dictionary
        pname = 'C10005'
        fgrp = 'AH '
        addc1dict(pchem, pname, fgrp)
        ierr = 0  # flag down, species found

        # add chemistry
        #       ..... WRITE THE REACTIONS (INCLUDING PHASE TRANSFER) ...

        np = 0
        rxinit(r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)
        nref = 0
        ref = [''] * mxcom
        addref(progname, 'JO23PRSCOM', nref, ref, pchem)

        add1tonp(progname, pchem, np)
        s[np-1] = 2.0
        p[np-1] = 'HO    '
        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'CO2   '

        r[0] = pname
        arrh[0] = 0.1
        arrh[1] = 0.
        arrh[2] = 0.
        r[1] = 'HV   '
        idreac = 1
        nlabel = 4
        xlabel = arrh[0]
        rxwrit(mecu, r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe, ref)

    # ------------
    # CO(ONO2)(ONO2)
    # ------------
    elif pchem == 'CO(ONO2)(ONO2)':
        # add the species to the dictionary
        pname = 'C10006'
        fgrp = 'NN '
        addc1dict(pchem, pname, fgrp)
        ierr = 0  # flag down, species found

        # add chemistry
        #       ..... WRITE THE REACTIONS (INCLUDING PHASE TRANSFER) ...

        np = 0
        rxinit(r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)
        nref = 0
        ref = [''] * mxcom
        addref(progname, 'JO23PRSCOM', nref, ref, pchem)

        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'NO3   '
        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'CO2   '
        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'NO2   '

        r[0] = pname
        arrh[0] = 0.1
        arrh[1] = 0.
        arrh[2] = 0.
        r[1] = 'HV   '
        idreac = 1
        nlabel = 4
        xlabel = arrh[0]
        rxwrit(mecu, r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe, ref)

    # ------------
    # CO(ONO2)(OOH)
    # ------------
    elif pchem == 'CO(ONO2)(OOH)':
        # add the species to the dictionary
        pname = 'C10007'
        fgrp = 'NH '
        addc1dict(pchem, pname, fgrp)
        ierr = 0  # flag down, species found

        # add chemistry
        #       ..... WRITE THE REACTIONS (INCLUDING PHASE TRANSFER) ...

        np = 0
        rxinit(r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)
        nref = 0
        ref = [''] * mxcom
        addref(progname, 'JO23PRSCOM', nref, ref, pchem)

        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'HO    '
        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'CO2   '
        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'NO3   '

        r[0] = pname
        arrh[0] = 0.1
        arrh[1] = 0.
        arrh[2] = 0.
        r[1] = 'HV   '
        idreac = 1
        nlabel = 4
        xlabel = arrh[0]
        rxwrit(mecu, r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe, ref)

    # ------------
    # CO(OOH)(OO.)
    # ------------
    elif pchem == 'CO(OOH)(OO.)':
        # add the species to the dictionary
        pname = 'C10008'
        fgrp = '3H  '
        addc1dict(pchem, pname, fgrp)
        ierr = 0  # flag down, species found

        # add chemistry
        # decomposition
        np = 0
        rxinit(r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)
        nref = 0
        ref = [''] * mxcom  # ; CALL addref(progname,'XXXXXXXXXX',nref,ref,pchem)

        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'CO2   '
        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'HO    '

        r[0] = pname
        arrh[0] = 1.00E6
        arrh[1] = 0.
        arrh[2] = 0.
        rxwrit(mecu, r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe, ref)

    # ------------
    # CO(NO2)(OOOH)
    # ------------
    elif pchem == 'CO(NO2)(OOOH)':
        # add the species to the dictionary
        pname = 'C10009'
        fgrp = 'KNV '
        addc1dict(pchem, pname, fgrp)
        ierr = 0  # flag down, species found

        # add chemistry
        #       ..... WRITE THE REACTIONS (INCLUDING PHASE TRANSFER) ...
        np = 0
        rxinit(r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)
        nref = 0
        ref = [''] * mxcom
        addref(progname, 'JO23PRSCOM', nref, ref, pchem)

        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'HO2   '
        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'CO2   '
        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'NO2   '

        # assign photolysis @ 0.1 x jNO2
        r[0] = pname
        arrh[0] = 0.1
        arrh[1] = 0.
        arrh[2] = 0.
        r[1] = 'HV   '
        idreac = 1
        nlabel = 4
        xlabel = arrh[0]
        rxwrit(mecu, r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe, ref)

    # --------------
    # CO(OOOH)(OOOH)
    # --------------
    elif pchem == 'CO(OOOH)(OOOH)':
        # add the species to the dictionary
        pname = 'C10010'
        fgrp = 'KNV '
        addc1dict(pchem, pname, fgrp)
        ierr = 0  # flag down, species found
        #       ..... WRITE THE REACTIONS (INCLUDING PHASE TRANSFER) ...
        np = 0
        rxinit(r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)
        nref = 0
        ref = [''] * mxcom
        addref(progname, 'JO23PRSCOM', nref, ref, pchem)

        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'HO    '
        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'CO2   '
        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'HO2   '

        # assign photolysis @ 0.1 x jNO2
        r[0] = pname
        arrh[0] = 0.1
        arrh[1] = 0.
        arrh[2] = 0.
        r[1] = 'HV   '
        idreac = 1
        nlabel = 4
        xlabel = arrh[0]
        rxwrit(mecu, r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe, ref)
              
        # ------------
        # CO(OH)(OH)
        # ------------
    elif pchem == 'CO(OH)(OH)':
        # add the species to the dictionary
        pname = 'H2CO3'
        fgrp = 'D  '
        addc1dict(pchem, pname, fgrp)
        ierr = 0  # flag down, species found

        # add chemistry
        #       ..... WRITE THE REACTIONS (INCLUDING PHASE TRANSFER) ...

        np = 0
        rxinit(r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)
        nref = 0
        ref = [''] * mxcom
        addref(progname, 'JO23PRSCOM', nref, ref, pchem)

        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'HO    '
        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'CO2   '
        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'HO2   '

        r[0] = pname
        arrh[0] = 0.1
        arrh[1] = 0.
        arrh[2] = 0.
        r[1] = 'HV   '
        idreac = 1
        nlabel = 4
        xlabel = arrh[0]
        rxwrit(mecu, r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe, ref)

    # ------------
    # CH2(OH)(ONO2)
    # ------------
    elif pchem == 'CH2(OH)(ONO2)':
        # The chemistry of this compound has to be checked. It's produced during 
        # the oxidation of alkanes. For the moment, the constants are the same as CH3OH.
        # It's assumed that only a reaction with OH occurs.

        # add the species to the dictionary
        pname = 'NO1001'
        fgrp = 'NO  '
        addc1dict(pchem, pname, fgrp)
        ierr = 0  # flag down, species found

        # add chemistry
        # decomposition
        np = 0
        rxinit(r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)
        nref = 0
        ref = [''] * mxcom  # ; CALL addref(progname,'XXXXXXXXXX',nref,ref,pchem)

        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'HCOOH '
        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'NO2   '

        r[0] = pname
        r[1] = 'HO    '
        arrh[0] = 3.10E-12
        arrh[1] = 0.
        arrh[2] = 0.
        rxwrit(mecu, r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe, ref)

    # ------------
    # CH2(ONO2)(OOH)
    # ------------
    elif pchem == 'CH2(ONO2)(OOH)':
        # The chemistry of this compound has to be checked. It's produced by 
        # the CH2(.OO.)+HNO3 reaction. The rate constant is based on Jenkin et
        # al., 2018 SAR and assume H abstraction from the CH2 group only. 
        # Reaction products are preliminary. 

        # add the species to the dictionary
        pname = 'NH1001'
        fgrp = 'NH  '
        addc1dict(pchem, pname, fgrp)
        ierr = 0  # flag down, species found

        # add chemistry
        np = 0
        rxinit(r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)
        nref = 0
        ref = [''] * mxcom
        addref(progname, 'MJ18KMV000', nref, ref, pchem)
        addref(progname, 'NH1001COM1', nref, ref, pchem)

        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'HCOO2H '
        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'NO3   '
        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'H2O   '

        r[0] = pname
        r[1] = 'HO    '
        arrh[0] = 4.95E-12
        arrh[1] = 0.
        arrh[2] = 719.
        rxwrit(mecu, r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe, ref)

    # ------------
    # CH2(NO2)(OO.)
    # ------------
    # The chemistry of this compound still need to be developped. The 
    # species is produced by the CH2(.OO.)+NO2 reaction. To avoid the 
    # accumulation of the radical, decomposition is performed with an 
    # arbitrary high rate constant.  
    elif pchem == 'CH2(NO2)(OO.)':
        # add the species to the dictionary
        pname = '2V1001'
        fgrp = '2V  '
        addc1dict(pchem, pname, fgrp)
        ierr = 0  # flag down, species found

        # add chemistry
        # decomposition
        np = 0
        rxinit(r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)
        nref = 0
        ref = [''] * mxcom
        addref(progname, 'RO2TEMPO', nref, ref, pchem)

        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'CH2O  '
        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'NO3   '

        r[0] = pname
        arrh[0] = 1.00E6
        arrh[1] = 0.
        arrh[2] = 0.
        rxwrit(mecu, r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe, ref)

    # ------------
    # CHO(ONO2)
    # ------------
    elif pchem == 'CHO(ONO2)':
        # this species is produced as a minor product during the oxidation of 
        # some cyclic HC. We assume that this species behaves like R-CO(ONO2)

        # add the species to the dictionary
        pname = 'C1ND00'
        fgrp = 'N  '
        addc1dict(pchem, pname, fgrp)
        ierr = 0  # flag down, species found

        # add chemistry
        # decomposition
        np = 0
        rxinit(r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)
        nref = 0
        ref = [''] * mxcom  # ; CALL addref(progname,'XXXXXXXXXX',nref,ref,pchem)

        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'CO2   '
        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'HO2   '
        add1tonp(progname, pchem, np)
        s[np-1] = 1.0
        p[np-1] = 'NO2   '

        r[0] = pname
        arrh[0] = 1.00E1
        arrh[1] = 0.
        arrh[2] = 0.
        rxwrit(mecu, r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe, ref)

    # ------------
    # CH3(ONO)
    # ------------
    #  elif pchem == 'CH3(ONO)':
    #
    # !!-- add the species to the dictionary
    #    pname = 'N01003'
    #    fgrp = '   '
    #    addc1dict(pchem, pname, fgrp)
    #    ierr = 0  # flag down, species found
    #
    # !!-- add chemistry
    # !! reaction with OH
    #    np = 0
    #    rxinit(r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)
    #    nref = 0
    #    ref = [''] * mxcom  # ; CALL addref(progname,'XXXXXXXXX',nref,ref,pchem)
    #
    #    add1tonp(progname, pchem, np)
    #    s[np-1] = 0.55
    #    p[np-1] = 'CH2O  '
    #    add1tonp(progname, pchem, np)
    #    s[np-1] = 0.45
    #    p[np-1] = 'NO    '
    #
    #    r[0] = pname
    #    r[1] = 'HO '
    #    arrh[0] = 0.301E-12
    #    arrh[1] = 0.
    #    arrh[2] = 0.
    #    rxwrit(mecu, r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe, ref)
    #    
    # !! photolysis
    #    np = 0
    #    rxinit(r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe)
    #    nref = 0
    #    ref = [''] * mxcom  # ; CALL addref(progname,'XXXXXXXXX',nref,ref,pchem)
    #
    #    add1tonp(progname, pchem, np)
    #    s[np-1] = 0.55
    #    p[np-1] = 'CH3O  '
    #    add1tonp(progname, pchem, np)
    #    s[np-1] = 0.45
    #    p[np-1] = 'NO    '
    #
    #    r[0] = pname
    #    idreac = 1
    #    nlabel = 40000   
    #    arrh[0] = 1.00E00
    #    arrh[1] = 0.
    #    arrh[2] = 0.
    #    rxwrit(mecu, r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe, ref)

# ! ======================================================================
# ! PURPOSE: Add the formula (pchem) and name (pname) provided as input 
# ! in the dictionaries.
# ! ======================================================================
def addc1dict(pchem, pname, fgrp):
    from toolbox import stoperr, addref
    from keyparameter import mxlfl, mxlco, mxlfo, mxpd, mecu
    from references import mxlcod
    from keyflag import enolflg, isomerfg
    from cdtool import switchenol
    from atomtool import cnum, onum
    from searching import srch, srh5
    from database import nspsp, dictsp
    from dictstackdb import (nrec, ninorg, inorglst, dict, namlst, dbrch, diccri,
                            stabl, lotopstack, nhldvoc, holdvoc, nhldrad, holdrad,
                            level)
    from namingtool import naming
    from switchisomer import isomer
    from tempflag import iflost
    from normchem import stdchm
    from rxwrttool import rxwrit, rxinit
    # 初始化变量
    progname = 'addc1dict'
    mesg = ''
    dicptr = 0
    namptr = 0

    # get pointer to add the species in namlst and dict
    dicptr = srch(nrec, pchem, dict)
    namptr = srh5(pname, namlst, nrec)
    if dicptr > 0:
        mesg = "formula for a C1 species to be added to dict already exists"
        stoperr(progname, mesg, pchem)
    if namptr > 0:
        mesg = "short name assigned to a C1 species already exists"
        stoperr(progname, mesg, pchem)

    # raise the counters
    nrec += 1
    if nrec > len(dict):
        mesg = "Too many species added to dictionary"
        stoperr(progname, mesg, pchem)

    # add new name and raise the name array
    namptr = abs(namptr) + 1
    # namlst(namptr+1:nrec+1)=namlst(namptr:nrec)
    namlst[namptr : nrec+1] = namlst[namptr-1 : nrec]
    namlst[namptr-1] = pname  # namlst(namptr) -> namlst[namptr-1]

    # raise dictionary and branching table. Add new line for the new species 
    dicptr = abs(dicptr) + 1
    # dict(dicptr+1:nrec+1)=dict(dicptr:nrec)
    dict[dicptr : nrec+1] = dict[dicptr-1 : nrec]
    # dbrch(dicptr+1:nrec+1)=dbrch(dicptr:nrec)
    dbrch[dicptr : nrec+1] = dbrch[dicptr-1 : nrec]
    # WRITE(dict(dicptr),'(a6,3x,a120,2x,a15)') pname, pchem, fgrp
    dict[dicptr-1] = f"{pname:<6}   {pchem:<120}  {fgrp:<15}"
    dbrch[dicptr-1] = 1.0  # set max yield to 1 (C1 species, number not needed)

    # raise index of isomer table
    if isomerfg:
        # diccri(dicptr+1:nrec+1,:) = diccri(dicptr:nrec,:)
        diccri[dicptr : nrec+1] = diccri[dicptr-1 : nrec]
        diccri[dicptr-1] = [0] * len(diccri[0]) if diccri else []