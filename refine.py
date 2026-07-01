import numpy
# FLASH refinement calculator

nyb = 8 # number of cells per block
nblocky = 16  # number of blocks in y
refine = 4  # amr refinement level 
domain = 1500e-4 # your domain length, cm
resolve = 10e-4 # what you want to resolve cm
x = domain/(nyb*nblocky*(2**(refine-1)))
ncells = resolve/x

print("Number of cells in ",resolve*1e4,"um = ", ncells)