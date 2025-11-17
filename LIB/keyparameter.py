#Set the list of key parameters for gecko
from math import sqrt
mxlco=6         # maximum length of species names (code)
mxlfo=120       # maximum length of a formula
mxlfl=15        # maximum length of functionality list (string)
mxldi=146       # string length in the dictionary
mxnode=35       # maximum number of nodes allowed
mxlgr=21        # maximum length of a string in a group
mxring=4        # maximum rings allowed
mxcp=99         # maximum # of copies of formula allowed & max # of "tracks" 
mxhyd=10        # maximum # of distinct position to add H2O (i.e. # of carbonyls)
mxhiso=5000     # maximum # of hydrate isomer a molecule can have
mxpd=38         # maximum # of products per reaction
mxnr=35         # maximum # of reactions per species
mxcopd=10       # maximum # of coproducts per generated species
mxnp=4          # maximum # of product per reaction in the output mechanism
mxps=600        # maximum # of "primary" species that can be given as input 
mxtrk=5         # maximum # of distinct Cd track in a molecule
mxlcd=4         # maximum length of a Cd track (current max is 4 for C=C-C=C)
mxlest=6        # maximum length of an ester track (current max is 6: CO-O-CO-O-CO-O)
  
# file unit
dctu=7          # dictionary file unit
prmu=14         # list of primary species used (for findname)
logu=15         # log file unit
refu=16         # mechanism with reference file unit
mecu=17         # mechanism file unit
scru=18         # file unit to redirect info from screen to file
kohu=19         # kivoci file unit (record k_OH for the various species) 
kno3u=59        # kjvocj file unit (record k_NO3 for the various species) 
waru=20         # warning file unit 
gasu=21         # gas phase species file unit
prtu=22         # particle phase species file unit
walu=23         # wall phase species file unit
ohu =26         # rate constant OH file unit (e.g. SAR assessment)
o3u =27         # rate constant O3 file unit (e.g. SAR assessment)
no3u=28         # rate constant NO3 file unit (e.g. SAR assessment)
dhfu=29         # dummy info file for dhf reactions
saru=30         # SAR info about group ...
tfu1=31         # temporary file unit 1
tfu2=32         # temporary file unit 2
tfu3=33         # temporary file unit 3
tfu4=34         # temporary file unit 4
mcou=37         # mechanism file unit with comment#
pfu1=41         # peroxy file unit 1
pfu2=42         # peroxy file unit 2
pfu3=43         # peroxy file unit 3
pfu4=44         # peroxy file unit 4
pfu5=45         # peroxy file unit 5
pfu6=46         # peroxy file unit 6
pfu7=47         # peroxy file unit 7
pfu8=48         # peroxy file unit 8
pfu9=49         # peroxy file unit 9
nmlu=70         # namelist file unit
rszu=789        # record stack size unit (information file)

# directory for the input/output files
dirgecko  =""
dirout    ="" 

# digit (to rate groups and ring joining char)     
digit =('1','2','3','4')                                                 

# priorities to rate/sort functionalities in groups
pri='CHOCHC CH(CH2 CH3O.)OO. OOH ONOF  Br2Br)Cl2Cl)NO2NO)OH)'

def generate_prime(x):
    temp=[]
    for i in range(2,x+1):
        is_prime=True
        for j in range(2,sqrt(i)+1):
            if i%j==0: 
                is_prime=False
                break
        if is_prime :
            temp.append(i)
    return tuple(temp)
        
    
# prime numbers (less than 1000)
prim=generate_prime(1000);
