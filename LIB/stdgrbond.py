import keyparameter,toolbox
#=======================================================================
# PURPOSE: Contruction of the group matrix and the carbon-carbon bond   
# matrix of a molecule (chem) given as input. Characters (/, \) used to 
# mark cis/trans C=C bonds in chem are removed in groups. The zebond 
# matrix (optional argument#) handle the cis/trans configuration.  
#                                                                   
# First groups are constructed. A group starts with a carbon "C", an
# aromatic carbon "c", or an ether bond "-O-", and ends with either 
# a next "C" or "c", or a double-bond "=", or an ether bond "-O-", 
# or  blank "NUL". In each group the trailing parentheses are deleted 
# if it is an opening with no closing and vice-versa. Next the carbon 
# skeleton is built with carbon, oxygen (if ether) and remaining 
# parentheses and double-bonds. The bond matrix is then built based on 
# the skeleton. The program checks the valence for each carbon center.            
#=======================================================================
def grbond (chem,group,bond,dbflg,nring,zebind) :
