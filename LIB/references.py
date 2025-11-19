
mxreac_info=1000 # max # of comment/references 
mxlcod=10    # max length of comment's code in database
mxlreac_info=400 # max length of a comment/reference
  
nreac_info=0         # # of available "reac_info" (references)  
code = [' ' * mxlcod for _ in range(mxreac_info)]  # Code for each "reac_info" 
reac_info = [' ' * mxlreac_info for _ in range(mxreac_info)]# Full text correxponding to the reac_info 