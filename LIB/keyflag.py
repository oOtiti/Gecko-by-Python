critvp = None        # log(Pvap) below which chemistry is ignored
# various cut-offs
brcut = None         # rx cut off branching ratio below threshold
yldcut = None        # (branching ratio * yield) cut off threshold
rxloss = None        # species yield below which is considered lost carbon
maxgen = None        # maximum # of generations allowed
TK = None            # user-set temperature
Tref = 298.0         # reference temperature（不可修改，对应 Fortran 的 PARAMETER）
Mref = 2.5e19        # default 3rd body M (molec cm-3)（不可修改）


# -- CHEMISTRY FLAG - SAR SELECTOR
kisom_sar = None     # R(O.) isomerisation rate selector (1=Atkinson 2007; 2=Vereecken 2009)
kdiss_sar = None     # R(O.) decomposition rate selector (1=Atkinson 2007; 2=Vereecken 2009)
kohadd_sar = None    # OH addition to alkene selector (1=Peeters 1997; 2=Ziemann 2009; 3=Jenkin 2016)
pvap_sar = None      # vapor pressure SAR selector: 1=JR-MY; 2=Nannoolal; 3=SIMPOL-1


# -- FLAGS
g2pfg = None                  # gas/particle partitioning flag
g2wfg = None                  # gas/wall partitioning flag (chamber simulations)
isomerfg = None               # allow isomer substitution
enolflg = None                # allow switching enols to ketones
pamfg = None                  # "PAM" simulations (H chemistry considered)
highnoxfg = None              # only NO chemistry for RO2
rx_ro_no2 = None              # add RO+NO2 reaction
C_NO2 = None                  # NO2 concentration (molec.cm-3) for RO+NO2 reaction
rx_ro2_no2 = None             # add RO2+NO2 reversible reaction
dhffg = None                  # activate DHF formation
chafg = None                  # consider cyclic hemi acetal as non-volatile
bimolecrx4criegee = None      # allow bimolecular reaction for stabilized criegee
rx_ro2_multiclass = None      # treat RO2+RO2 with 9 classes (True) or only CH3O2 (False)
rx_ro2_oh = None              # add RO2+OH reaction


# -- INFOS
sar_only_fg = None            # write kOH/kNO3/kO3 for SAR assessment
wrtsarinfo = None             # write SAR info (groups, ...) in files
wrtopeinfo = False            # write "operator info" (OLD STUFF, Fortran 中是 PARAMETER)


# -- OUTPUT FLAGS
screenfg = None               # write additional info on screen
wrtdhf = None                 # write formation enthalpies file
wrtpvap = None                # write vapor pressures file
wrthenry = None               # write Henry's law coefs file（原 Fortran 中是 wrthenry，修正拼写）
wrtref = None                 # write reaction references
wrttg = None                  # write Tg data file
wrtdepo = None                # write deposition data file
wrtkivoci = None              # write kivoci and kjvocj file
wrtmaxyield = None            # write maximum yields file
#-----------------------------------
# default environmental parameters for mechanism generation
# => values to use if not supplied by user
#-----------------------------------
import keyparameter

def define_defaults() :
    global critvp, brcut, yldcut, rxloss, maxgen, TK, C_NO2, \
           rx_ro_no2, rx_ro2_no2, rx_ro2_multiclass, \
           pvap_sar, kisom_sar, kdiss_sar, kohadd_sar, \
           g2pfg, g2wfg, isomerfg, highnoxfg, dhffg, chafg, \
           enolflg, pamfg, bimolecrx4criegee, rx_ro2_oh, \
           wrtsarinfo, screenfg, wrtdhf, wrtpvap, wrthenry, wrtref, \
           wrttg, wrtdepo, wrtkivoci, wrtmaxyield, sar_only_fg
# various cut-offs
    maxgen = 20              # maximum # of generations allowed
    TK     = 298.0            # temperature for which to generate the mechanism
    critvp = -13             # log(Pvap) below which chemistry is ignored
    brcut  = 0.05            # rx cut off branching ratio below which a reaction pathway is ignored
    yldcut = 1e-3            # rx cut off if: (rx branching ratio)*(yield) is below threshold
    rxloss = 1e-10           # send to lost carbon if species yield below threshold
    C_NO2  = 2.5e10          # NO2 concentration (molec.cm-3) used in roselector (rochem) for RO+NO2 reaction

# reaction switches
    rx_ro_no2 = False      # add the RO+NO2 reaction to the RO chemistry
    rx_ro2_no2= False      # add the RO2+NO2 reversible reaction to the RO2 chemistry
    rx_ro2_multiclass=True # treat RO2+RO2 with the 9 classes (true) of ro2 or only CH3O2 (false)

# SAR 
    pvap_sar   = 2           # SAR for vapor pressure: 1=JR-MY; 2=Nannoolal; 3=SIMPOL-1
    kisom_sar  = 2           # R(O.) selector for isomerisation reaction rate (1=Atkinson 2007; 2=Vereecken 2009)
    kdiss_sar  = 2           # R(O.) selector for decomposition reaction rate (1=Atkinson 2007; 2=Vereecken 2009)
    kohadd_sar = 3           # OH addition to alkene (1=Peeters 1997; 2=Ziemann 2009; 3=Jenkin 2016)
    
# flags
    g2pfg      = True      # flag to write for gas/particle partitioning "reactions"
    g2wfg      = False     # flag to write for gas/wall partitioning "reactions" (chamber simulations)
    isomerfg   = True      # allow isomer substitution
    highnoxfg  = False     # only NO chemistry for RO2 considered 
    dhffg      = False     # flag to activate DHF formation 
    chafg      = False     # consider cyclic hemi acetal (CHA) species as non volatile
    pamfg      = False     # "PAM" simulations: H chemistry could be considered
    bimolecrx4criegee = False # allow bimolecular reaction for stabilized criegee (keep H2O)
    enolflg    = True      # allow to switch enols to ketones -------------- WARNING --------- /#\ GECKO-A does NOT treat enol for now /#\
    rx_ro2_oh  = False     # add the RO2+OH reaction in RO2 chemistry  ----- WARNING --------- /#\ if true, leads to ROOOH product for which GECKO has a really simple chemistry so far /#\

# Debug/infos
    wrtsarinfo = False     # write info about SAR (groups, ...) in dedicated files  
    sar_only_fg= False     # write kOH, kNO3, kO3 for SAR assessment (rate in database not longer used) - if true, run the reactions for the input species only

# default directories
    keyparameter.dirgecko  = '../'
    keyparameter.dirout = 'OUT/'

# write output files
    screenfg   = False     # write additional info on screen during generation
    wrtdhf     = False     # write output file with formation enthalpies information
    wrtpvap    = True      # write output file with the vapor pressures
    wrthenry   = True      # write output file with the Henry's law coefs
    wrtref     = True      # write the references for the reactions in the mechanism
    wrttg      = True      # write output file with Tg data
    wrtdepo    = False     # write output file with data fo deposition
    wrtkivoci  = False     # write output file with kivoci and kjvocj
    wrtmaxyield= False     # write output file with maximum yields

define_defaults()