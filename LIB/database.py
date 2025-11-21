import keyparameter
import references
# 常量定义
mxkdb = 1000  # max # of rate constant in a given database
mxkwr = 110   # max # of known reaction in a given database
mxbg = 500    # max # of benson group
lbg = 24      # max length of a benson group (string)
mxspsp = 3000 # max # of special species
mxog = 10000  # max # of entry (reaction) in the oge database
mxfn = 3000   # max # of species with fixed names (i.e. not set by gecko)

# ==============================================================================
# Rate constants for VOC+oxidant
# ==============================================================================

# VOC+OH database
nkohdb = 0  # # of species in the database
# formula of the species (长度为 keyparameter.mxlfo)
kohdb_chem = [' ' * keyparameter.mxlfo for _ in range(mxkdb)]
kohdb_298 = [0.0 for _ in range(mxkdb)]  # rate constant @ 298 K
# arrhenius parameter for the rate constant (A, n, Ea)
kohdb_arr = [[0.0 for _ in range(3)] for _ in range(mxkdb)]
# comment's code for the rate (长度为 references.mxlcod)
kohdb_com = [[' ' * references.mxlcod for _ in range(3)] for _ in range(mxkdb)]

# VOC+O3 database
nko3db = 0  # # of species in the database
# formula of the species (长度为 keyparameter.mxlfo)
ko3db_chem = [' ' * keyparameter.mxlfo for _ in range(mxkdb)]
ko3db_298 = [0.0 for _ in range(mxkdb)]  # rate constant @ 298 K
# arrhenius parameter for the rate constant (A, n, Ea)
ko3db_arr = [[0.0 for _ in range(3)] for _ in range(mxkdb)]
# comment's code for the rate (长度为 references.mxlcod)
ko3db_com = [[' ' * references.mxlcod for _ in range(3)] for _ in range(mxkdb)]

# VOC+NO3 database
nkno3db = 0  # # of species in the database
# formula of the species (长度为 keyparameter.mxlfo)
kno3db_chem = [' ' * keyparameter.mxlfo for _ in range(mxkdb)]
kno3db_298 = [0.0 for _ in range(mxkdb)]  # rate constant @ 298 K
# arrhenius parameter for the rate constant (A, n, Ea)
kno3db_arr = [[0.0 for _ in range(3)] for _ in range(mxkdb)]
# comment's code for the rate (长度为 references.mxlcod)
kno3db_com = [[' ' * references.mxlcod for _ in range(3)] for _ in range(mxkdb)]

# ==============================================================================
# Known mechanism for VOC+oxidant ("kw" data base)
# ==============================================================================

# VOC+OH mechanism
nkwoh = 0  # # of species in the database
# formula of the VOC reacting with OH (长度为 keyparameter.mxlfo)
kwoh_rct = [' ' * keyparameter.mxlfo for _ in range(mxkwr)]
nkwoh_pd = [0 for _ in range(mxkwr)]  # # of channel involved in mechanism
# yield of each channel (mxkwr x keyparameter.mxnr)
kwoh_yld = [[0.0 for _ in range(keyparameter.mxnr)] for _ in range(mxkwr)]
# formula of the main product per channel (长度为 keyparameter.mxlfo, mxkwr x keyparameter.mxnr)
kwoh_pd = [[' ' * keyparameter.mxlfo for _ in range(keyparameter.mxnr)] for _ in range(mxkwr)]
# coproducts per channel (长度为 keyparameter.mxlco, mxkwr x keyparameter.mxnr x keyparameter.mxcopd)
kwoh_copd = [[[' ' * keyparameter.mxlco for _ in range(keyparameter.mxcopd)] 
              for __ in range(keyparameter.mxnr)] for ___ in range(mxkwr)]
# comment's code for the rate (长度为 references.mxlcod)
kwoh_com = [[' ' * references.mxlcod for _ in range(3)] for _ in range(mxkwr)]

# VOC+O3 mechanism
nkwo3 = 0  # # of species in the database
# formula of the VOC reacting with OH (长度为 keyparameter.mxlfo)
kwo3_rct = [' ' * keyparameter.mxlfo for _ in range(mxkwr)]
nkwo3_pd = [0 for _ in range(mxkwr)]  # # of channel involved in mechanism
# yield of each channel (mxkwr x keyparameter.mxnr)
kwo3_yld = [[0.0 for _ in range(keyparameter.mxnr)] for _ in range(mxkwr)]
# formula of the main product per channel 3 (长度为 keyparameter.mxlfo, mxkwr x keyparameter.mxnr)
kwo3_pd = [[' ' * keyparameter.mxlfo for _ in range(keyparameter.mxnr)] for _ in range(mxkwr)]
# coproducts per channel (长度为 keyparameter.mxlco, mxkwr x keyparameter.mxnr x keyparameter.mxcopd)
kwo3_copd = [[[' ' * keyparameter.mxlco for _ in range(keyparameter.mxcopd)] 
              for __ in range(keyparameter.mxnr)] for ___ in range(mxkwr)]
# comment's code for the rate (长度为 references.mxlcod)
kwo3_com = [[' ' * references.mxlcod for _ in range(3)] for _ in range(mxkwr)]

# VOC+NO3 mechanism
nkwno3 = 0  # # of species in the database
# formula of the VOC reacting with OH (长度为 keyparameter.mxlfo)
kwno3_rct = [' ' * keyparameter.mxlfo for _ in range(mxkwr)]
nkwno3_pd = [0 for _ in range(mxkwr)]  # # of channel involved in mechanism
# yield of each channel (mxkwr x keyparameter.mxnr)
kwno3_yld = [[0.0 for _ in range(keyparameter.mxnr)] for _ in range(mxkwr)]
# formula of the main product per channel 3 (长度为 keyparameter.mxlfo, mxkwr x keyparameter.mxnr)
kwno3_pd = [[' ' * keyparameter.mxlfo for _ in range(keyparameter.mxnr)] for _ in range(mxkwr)]
# coproducts per channel (长度为 keyparameter.mxlco, mxkwr x keyparameter.mxnr x keyparameter.mxcopd)
kwno3_copd = [[[' ' * keyparameter.mxlco for _ in range(keyparameter.mxcopd)] 
               for __ in range(keyparameter.mxnr)] for ___ in range(mxkwr)]
# comment's code for the rate (长度为 references.mxlcod)
kwno3_com = [[' ' * references.mxlcod for _ in range(3)] for _ in range(mxkwr)]

# ==============================================================================
# Known mechanism for radicals (RO2, RCOO2, RO and criegee)
# ==============================================================================

# RO2 chemistry
nkwro2 = 0  # # of entry (reaction) in the database
# arrhenius parameter for the rate constant (A, n, Ea)
kwro2_arrh = [[0.0 for _ in range(3)] for _ in range(mxkwr)]
# stoi. coef. of the products
kwro2_stoi = [[0.0 for _ in range(4)] for _ in range(mxkwr)]
# formula of the reactants (长度为 keyparameter.mxlfo, mxkwr x 2)
kwro2_rct = [[' ' * keyparameter.mxlfo for _ in range(2)] for _ in range(mxkwr)]
# formula of the products (长度为 keyparameter.mxlfo, mxkwr x 4)
kwro2_prd = [[' ' * keyparameter.mxlfo for _ in range(4)] for _ in range(mxkwr)]
# comment's code for the ro2 reaction (长度为 references.mxlcod, mxkwr x 4)
kwro2_com = [[' ' * references.mxlcod for _ in range(4)] for _ in range(mxkwr)]

# RCOO2 chemistry
nkwrco3 = 0  # # of entry (reaction) in the database
# arrhenius parameter for the rate constant (A, n, Ea)
kwrco3_arrh = [[0.0 for _ in range(3)] for _ in range(mxkwr)]
# stoi. coef. of the products
kwrco3_stoi = [[0.0 for _ in range(4)] for _ in range(mxkwr)]
# formula of the reactants (长度为 keyparameter.mxlfo, mxkwr x 2)
kwrco3_rct = [[' ' * keyparameter.mxlfo for _ in range(2)] for _ in range(mxkwr)]
# formula of the products (长度为 keyparameter.mxlfo, mxkwr x 4)
kwrco3_prd = [[' ' * keyparameter.mxlfo for _ in range(4)] for _ in range(mxkwr)]
# comment's code for the rcoo2 reaction (长度为 references.mxlcod, mxkwr x 4)
kwrco3_com = [[' ' * references.mxlcod for _ in range(4)] for _ in range(mxkwr)]

# RO chemistry
nkwro = 0  # # of entry (reaction) in the database
# arrhenius parameter for the rate constant (A, n, Ea)
kwro_arrh = [[0.0 for _ in range(3)] for _ in range(mxkwr)]
# stoi. coef. of the products
kwro_stoi = [[0.0 for _ in range(4)] for _ in range(mxkwr)]
# formula of the reactants (长度为 keyparameter.mxlfo, mxkwr x 2)
kwro_rct = [[' ' * keyparameter.mxlfo for _ in range(2)] for _ in range(mxkwr)]
# formula of the products (长度为 keyparameter.mxlfo, mxkwr x 4)
kwro_prd = [[' ' * keyparameter.mxlfo for _ in range(4)] for _ in range(mxkwr)]
# comment's code for the ro reaction (长度为 references.mxlcod, mxkwr x 4)
kwro_com = [[' ' * references.mxlcod for _ in range(4)] for _ in range(mxkwr)]

# Criegee chemistry
nkwcri = 0  # # of entry (reaction) in the database
# arrhenius parameter for the rate constant (A, n, Ea)
kwcri_arrh = [[0.0 for _ in range(3)] for _ in range(mxkwr)]
# stoi. coef. of the products
kwcri_stoi = [[0.0 for _ in range(4)] for _ in range(mxkwr)]
# formula of the reactants (长度为 keyparameter.mxlfo, mxkwr x 2)
kwcri_rct = [[' ' * keyparameter.mxlfo for _ in range(2)] for _ in range(mxkwr)]
# formula of the products (长度为 keyparameter.mxlfo, mxkwr x 4)
kwcri_prd = [[' ' * keyparameter.mxlfo for _ in range(4)] for _ in range(mxkwr)]
# comment's code for the ro reaction (长度为 references.mxlcod, mxkwr x 3)
kwcri_com = [[' ' * references.mxlcod for _ in range(3)] for _ in range(mxkwr)]

# ==============================================================================
# photolysis labels & known photolysis reactions
# ==============================================================================
njdat = 0  # # of known photo. reactions
# formula of the species being photolyzed (长度为 keyparameter.mxlfo)
jchem = [' ' * keyparameter.mxlfo for _ in range(mxkwr)]
# the two main photodissociation fragments (长度为 keyparameter.mxlfo, mxkwr x 2)
jprod = [[' ' * keyparameter.mxlfo for _ in range(2)] for _ in range(mxkwr)]
jlabel = [0 for _ in range(mxkwr)]  # ID number of photolytic reaction i
# additional inorganic coproduct (长度为 keyparameter.mxlco)
coprodj = [' ' * keyparameter.mxlco for _ in range(mxkwr)]
# comment's code for the photolysis reactions (长度为 references.mxlcod, mxkwr x 3)
j_com = [[' ' * references.mxlcod for _ in range(3)] for _ in range(mxkwr)]
nj40 = 0  # # of J40 data in database
jlab40 = [0 for _ in range(mxkwr)]  # label for which J4O is provided
j40 = [0.0 for _ in range(mxkwr)]  # J40 values for various labels

# ==============================================================================
# Benson groups database
# ==============================================================================
nbson = 0  # # of group in the database
bsonval = [0.0 for _ in range(mxbg)]  # values of Benson group
bsongrp = [' ' * lbg for _ in range(mxbg)]  # Benson group string (长度为 lbg)

# ==============================================================================
# Simpol SAR for vapor pressure
# ==============================================================================
p0simpgrp298 = 0.0  # constant term for psat @ 298K
psimpgrp298 = [0.0 for _ in range(30)]  # grp contribution for psat @ 298K
h0simpgrp298 = 0.0  # constant term for heat of vaporisation @ 298K
hsimpgrp298 = [0.0 for _ in range(30)]  # grp contribution for heat of vaporisation @ 298K

# ==============================================================================
# Nannoolal SAR for vapor pressure and Tb
# ==============================================================================
# grp contribution. 1st column is Tb weight, 2nd is Psat (219 x 2)
nanweight = [[0.0 for _ in range(2)] for _ in range(219)]

# ==============================================================================
# Henry's law coefficients (HLC)
# ==============================================================================
nhlcdb = 0  # # number of species in the hlc database
# formula of the species (长度为 keyparameter.mxlfo)
hlcdb_chem = [' ' * keyparameter.mxlfo for _ in range(mxkdb)]
hlcdb_dat = [[0.0 for _ in range(3)] for _ in range(mxkdb)]  # data
# comment's code for the data (长度为 references.mxlcod, mxkdb x 3)
hlcdb_com = [[' ' * references.mxlcod for _ in range(3)] for _ in range(mxkdb)]

# ==============================================================================
# Hydration constants (for effective HLC)
# ==============================================================================
nkhydb = 0  # # number of species in the khyd database
# formula of the species (长度为 keyparameter.mxlfo)
khydb_chem = [' ' * keyparameter.mxlfo for _ in range(mxkdb)]
khydb_dat = [[0.0 for _ in range(1)] for _ in range(mxkdb)]  # data
# comment's code for the data (长度为 references.mxlcod, mxkdb x 3)
khydb_com = [[' ' * references.mxlcod for _ in range(3)] for _ in range(mxkdb)]

# ==============================================================================
# special species list (#species)
# ==============================================================================
nspsp = 0  # # of "special" species (e.g. furane)
lospsp = [False for _ in range(mxspsp)]  # logical, "true" if the species has been used.
# dictionary line of the special species (长度为 keyparameter.mxldi)
dictsp = [' ' * keyparameter.mxldi for _ in range(mxspsp)]

# ==============================================================================
# special species chemistry (oge = Out GEnerator)
# ==============================================================================
noge = 0  # # of reactions given (oge = "out generator")
ogelab = [0 for _ in range(mxog)]  # label for the reaction (if EXTRA or HV)
# reactants for reaction i (长度为 keyparameter.mxlfo, mxog x 3)
ogertve = [[' ' * keyparameter.mxlfo for _ in range(3)] for _ in range(mxog)]
# products for reaction i (长度为 keyparameter.mxlfo, mxog x keyparameter.mxpd)
ogeprod = [[' ' * keyparameter.mxlfo for _ in range(keyparameter.mxpd)] for _ in range(mxog)]
# arrehnius coefficients (A, n, Ea)
ogearh = [[0.0 for _ in range(3)] for _ in range(mxog)]
# stochiometric coefficients for product j
ogestoe = [[0.0 for _ in range(keyparameter.mxpd)] for _ in range(mxog)]
# aux. info. for reaction i (e.g. falloff react.)
ogeaux = [[0.0 for _ in range(7)] for _ in range(mxog)]

# ==============================================================================
# fixed names data (fn database)
# ==============================================================================
nfn = 0  # # of species with fixed names
# names of fixed species (长度为 keyparameter.mxlco)
namfn = [' ' * keyparameter.mxlco for _ in range(mxfn)]
# standardized formula of species with fixed name (长度为 keyparameter.mxlfo)
chemfn = [' ' * keyparameter.mxlfo for _ in range(mxfn)]