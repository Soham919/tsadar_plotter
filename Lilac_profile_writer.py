import numpy as np
from matplotlib import pyplot as plt
from scipy import constants
from GJS_analyze import GJS_analyze
from pathlib import Path

field = 'density'
fn = "D-GJ-C-232_H2N2_1500psi.h5"
full_data = GJS_analyze(fn, field)
r, z, q = full_data
q = q*(80/1500)  # scale the density to match whatever psi you want

# ------ Convert data to lower resolution for writing to input deck ------ #

step = 10  # resolution of the profile in index space
z1 = np.zeros((1,len(z)//step))
r1 = np.zeros((len(r)//step,1))
values = np.zeros((len(r)//step, len(z)//step))

# sample values at desired steps
for i in range(0, len(r), step):
    r1[i//step] = r[i]
    for j in range(0, len(z), step):
        z1[0,j//step] = z[j]
        values[i//step, j//step] = q[i, j] 

R1, Z1 = np.meshgrid(r1, z1)
values = values.T  # transpose to match the orientation of R1, Z1
d = (r1[1]-r1[0]).squeeze()*1000 # um, thickness of the cell in the r direction

# plot to see resolution 
zh = 5 # mm, height at which to take the profile
rd = 8 # mm, radius at which you want to put the foil
idr = np.abs(r1-rd).argmin()  # index of r closest to rd
idz = np.abs(z1-zh).argmin()  # index of z closest to zh
print(f"{z1[0,idz]}")
dens_profile = values[0:idr, idz]  # density profile at z = zh, for r <= rd
dens_profile = np.concatenate((dens_profile[::-1], dens_profile,)) # make it symmetric about r=0

fig, ax = plt.subplots(1,2, figsize=(12,5))
im = ax[0].pcolormesh(R1, Z1, values.T, cmap='viridis', shading='auto')
ax[0].axvline(x=r1[idr], color='red', linestyle='--', label=f'r = {rd} mm')
ax[0].axhline(y=z1[0,idz], color='red', linestyle='--', label=f'z = {zh} mm')
ax[0].set_xlabel('R (mm)')
ax[0].set_ylabel('Z (mm)')

ax[1].plot(dens_profile/1000, color='blue', label=f'z = {z1[0,idz]:.2f} mm')
ax[1].set_xlabel('R (mm)')
ax[1].set_ylabel('Density (g/cm^3)')
ax[1].legend()
plt.show()


# ----------- write to input deck ------------ #

# -- template for the input deck, with placeholders for the cell thickness and density -- #
t1 = """$rhydro
  igeom=-1
  istart=0,
  job_name='sban_kinshock26A',
  tcut=10,
  toutf=0.,
  iprint(6)=0,iprint(4)=0,iprint(15)=0,iprint(16)=0,
  ipfav=0,ipfs=0,
  idtime=1,fvol=0.05,
  igrwth=0,fflux=-0.06,
  ivoidc=2,
  iradt=3,itntran=0,ntemp=2,
  irpost=0,
  ippost(13)=3,  
  ipostf=0,tpostf=25.0,
  irayt=-5,zoning='user',
  ifreq(1)=3,
  ilasar=8,
  ps_file = '600ps_flat.txt',
  dtmax=1.e-13,
 $end  

 $prof 
  matcod=110,
  ncell=250,
  thick=500.0,
  itnuc=0,
  ieos=4,
  izion=1,
  dens_mat=1.0,
  temp=10.0,
 $end     
"""
template = """
$prof 
 matcod = 1, 7,
 ncell = 10,
 thick = {v1},
 itnuc = 0,
 ieos = 8, 4,
 mixfrac = 0.9, 0.1,
 dens_mat = {v2},
 temp = 10.0,
 $end
"""
t3 = """  
 $prof 
  matcod= 14, 7,
  ncell=300,
  thick=1.0,
  itnuc=0,
  ieos=4, 4,
  mixfrac = 0.428, 0.571,
  izion=1,
  dens_mat=3.2,
  temp=10.0,
 $end   
 
 $prof ncell=-1, 
 $end
 
 &rbeam
    beam_file = 'profile_sg5_sshv_SSD_415.dat',
 &end
"""
 # write the input deck with the profiles for each cell
with open("lilac_data_input_80psi.txt", "w") as f:
    
    f.write(t1) # write the first part of the input deck
    f.write("\n")
    
    # write the profiles for each cell
    for i in range(dens_profile.size-1):  
        v1 = d # thickness of the cell in um
        v2 = dens_profile[i]/1000  # density of the cell in g/cm^3
        f.write(template.format(v1=v1, v2=v2))
        f.write("\n")
    
    f.write(t3) # write the last part of the input deck
