import numpy
# FLASH refinement calculator

nyb = 8 # number of cells per block
nblocky = 2  # number of blocks in y
refine = 5  # amr refinement level 
domain = 0.1 # your domain length, cm
resolve = 5e-04 # what you want to resolve cm
#x = 2.5e-04
x = domain/(nyb*nblocky*(2**(refine-1)))
#domain = x*(nyb*nblocky*(2**(refine-1)))
ncells = resolve/x

print("x = ",x*(10**4),"um")
print("Number of cells in ",resolve*1e4,"um = ", ncells)
#print(domain)