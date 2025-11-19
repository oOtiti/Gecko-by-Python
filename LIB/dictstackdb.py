from keyparameter import mxldi, mxlco  # string length in the dictionary

# dictionaries 
mxspe=7000000  # max # of species allowed in the mechanism
mxc1=3000      # max # of inorganic and C1 species allowed in the mechanism

nrec=0               # # of species recorded in dict      
ninorg=0             # # of inorganic species
nwpspe=0             # # of species recorded in the condensed (wall,part.) phase
dict=[' ' * mxldi for _ in range(mxspe)]  # dictionary line (code+formula+fg)
namlst=[' ' * mxlco for _ in range(mxspe)]  # name (lco=6 characters) of the species used
inorglst=[' ' * mxldi for _ in range(mxc1)]  # list of inorganic species (code+formula+fg)
dbrch=[0.0 for _ in range(mxspe)]  # yield attach to a formula in dict
# to be introduced in a "type dict"
dctmw=[0.0 for _ in range(mxspe)]     # molecular weight of species in dict
dctatom=[[0 for _ in range(9)] for _ in range(mxspe)]  # atoms & radical for species in dict
dcthenry=[0.0 for _ in range(mxspe)]  # Henry's law coeff. for species in dict
dctnan=[[0.0 for _ in range(2)] for _ in range(mxspe)]  # Pvap(:,1) & heat(:,2) (Nannoonal) for species in dict
dctsim=[[0.0 for _ in range(2)] for _ in range(mxspe)]  # Pvap(:,1) & heat(:,2) (Simpol) for species in dict
dctmyr=[[0.0 for _ in range(2)] for _ in range(mxspe)]  # Pvap(:,1) & heat(:,2) (Myrdal&Yalkowski) for species in dict

# stack                                         
mxsvoc=1500000  # max # of species in the voc stack
mxsrad=200      # max # of species in the radical stack
lenss=132       # length of a stack string (code+formula+i3+i3)

nhldvoc=0             # # of (stable) VOC in the stack
holdvoc=[' ' * lenss for _ in range(mxsvoc)]  # VOCs in the stack (name[a6]+formula[a120]+stabl[i3]+level[i3])
nhldrad=0             # # of radical in the stack
holdrad=[' ' * lenss for _ in range(mxsrad)]  # radicals in the stack(name[a6]+formula[a120]+stabl[i3]+level[i3])
stabl=0               # # of generations needed to produce the current species
level=0               # # of intermediates (rad + stable) needed to produce the current species
lotopstack=False      # if true, add species on top of VOC stack (default is false)

# isomer stack
mxiso=2000      # max # of isomers for a given molecule 
mxcri=21        # max # of criteria used to discriminate isomers
diccri=[[0 for _ in range(mxcri)] for _ in range(mxspe)]  #

# tetrahydrofuran (CHA) stack (to manage phase partitioning)
mxcha=10000      # max # of tetrahydrofuran (CHA) allowed 
ncha=0               # # of species recorded in dict      
chatab=[' ' * mxlco for _ in range(mxcha)]    # idnam of the CHA 