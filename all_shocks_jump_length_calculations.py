import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
from Plasma_func import cs_ij, cs, RH, mfp_ii, f_rho

eps = constants.epsilon_0 # epislon naught SI
kb = constants.Boltzmann # Boltzmann constant SI
e = constants.e # charge of electron(Coulomb)
me = constants.electron_mass # mass electron(kg)
mp = constants.proton_mass # mass of H(kg)
pi = constants.pi # pi


############ 92535-92536-92537 #############
# H_rho = 4/30  # 10^20 cm^-3
# He_rho = 1/30 # 10^20 cm^-3
# Te = 50 # eV
# Te_s = 400 # eV
# Ti_s = 700/1000 # KeV

# us = (4.5/6)*10**6
# cs_H_He2 = cs_ij(H_rho, He_rho, 1, 4, 1, 2,Te)
# cs_H_He1 = cs_ij(H_rho, He_rho, 1, 4, 1, 1,Te)
# M2 = us/cs_H_He2
# M1 = us/cs_H_He1
# print("Mach number for He 2+ =", M2)
# print("Mach number for He 1+ =", M1)

# [rho2,u2,T2] = RH(3.74,5/30,40,50)  # actual Mach number seems to be 3.74
# print("rho2 = ",rho2,"x 10^20 cm^-3")
# print("u2 = ",u2,"x 10^6 cm/s")
# print("T2 = ",T2/1000,"KeV")

# ## PEDESTAL LENGTH 
# hconc = np.linspace(0.01, 0.99, 100)
# Temp = [0.5, 0.6, 0.7]  # KeV
# df = [0.35,0.25,0.15]  # ne
# mfp_12 = np.zeros((len(df),len(Temp),len(hconc)))
# mfp_11 = np.zeros((len(df),len(Temp),len(hconc)))
# mfp_1 = np.zeros((len(df),len(Temp),len(hconc)))

# ##### FOR H-HE PLASMA PLOT ######
# for j in range(len(df)):
#     for i in range(len(Temp)):
#         mfp_12[j,i,:] = mfp_ii(1,4,1,2,Temp[i],f_rho(hconc,1,2,df[j])[1])
#         mfp_11[j,i,:] = mfp_ii(1,1,1,1,Temp[i],f_rho(hconc,1,2,df[j])[0])
#         mfp_1[j,i,:] = ((1/mfp_12[j,i,:])+(1/mfp_11[j,i,:]))**(-1)

# p_length = ((mp/me)**(0.5))*mfp_1

# ##### FOR H-N PLASMA PLOT ######
# # for j in range(len(df)):
# #     for i in range(len(Temp)):
# #         mfp_12[j,i,:] = mfp_ii(1,14,1,4,Temp[i],f_rho(hconc,1,4,df[j])[1])
# #         mfp_11[j,i,:] = mfp_ii(1,1,1,1,Temp[i],f_rho(hconc,1,4,df[j])[0])
# #         mfp_1[j,i,:] = ((1/mfp_12[j,i,:])+(1/mfp_11[j,i,:]))**(-1)
        
# # p_length = ((mp/me)**(0.5))*mfp_1

# ## PEDESTAL LENGTH PLOTS
# fig, ax = plt.subplots()
# #plt.plot(hconc,mfp_hhe,color= 'r',label='H-He')
# #plt.plot(hconc,mfp_hh2,color = 'b',label='H-H mfp')

# ax.plot(hconc,p_length[0,0,:],color='red',label=r"ne = 0.35 x $10^{20} cm^{-3}$",linestyle='-')
# ax.plot(hconc,p_length[0,1,:],color='red')
# ax.plot(hconc,p_length[0,2,:],color='red')

# ax.plot(hconc,p_length[1,0,:],color='blue',label=r"ne = 0.2 x $10^{20} cm^{-3}$",linestyle='-')
# ax.plot(hconc,p_length[1,1,:],color='blue')
# ax.plot(hconc,p_length[1,2,:],color='blue')

# ax.plot(hconc,p_length[2,0,:],color='mediumseagreen',label=r"ne = 0.15 x $10^{20} cm^{-3}$",linestyle='-')
# ax.plot(hconc,p_length[2,1,:],color='mediumseagreen')
# ax.plot(hconc,p_length[2,2,:],color='mediumseagreen')
#     #ax.scatter(Temp[i],p_length[i,50],label=f"T_{i} = {Temp[i]} KeV")
# ax.set_xlabel("H conc.(/1)")
# ax.set_ylabel(r"$\sqrt{\frac{m_i}{m_e}}\lambda_{ij}$ (mm)")
# ax.set_title(r"Pedestal length for $H^{+}-He^{+2}$ plasma shock")
# ax.text(0.6, 0.8, "Temp = 0.5, 0.6, 0.7 KeV",
#         transform=ax.transAxes,
#         fontsize=10,
#         color='black',
#         bbox=dict(facecolor='white', edgecolor='black', boxstyle='square,pad=0.5'))

# ax.legend()
# plt.show()

############ 92527-92530 ############
ne_s = 0.4
H_rho = 0.4  # 10^20 cm^-3
Te = 50 # eV
Te_s = 500 # eV
Ti_s = 700/1000 # KeV

us = (5/6)*10**6
cs_H_H = cs(1,1,50)
M = us/cs_H_H

print("Mach number for 92527-92530 shock =", M)


[rho2,u2,T2] = RH(M,H_rho,30,50)  # actual Mach number seems to be 3.74
print("rho2 = ",rho2,"x 10^20 cm^-3")
print("u2 = ",u2,"x 10^6 cm/s")
print("T2 = ",T2/1000,"KeV")

## PEDESTAL LENGTH 
Temp = [0.5, 0.6, 0.75]  # KeV
df = [0.4,0.35]  # ne
mfp_1 = np.zeros([len(df),len(Temp)])
for j in range(len(df)):
    for i in range(len(Temp)):
        mfp_1[j,i] = mfp_ii(1,1,1,1,Temp[i],df[j])

p_length = ((mp/me)**(0.5))*mfp_1

## PEDESTAL LENGTH PLOTS
fig, ax = plt.subplots()
#plt.plot(hconc,mfp_hhe,color= 'r',label='H-He')
#plt.plot(hconc,mfp_hh2,color = 'b',label='H-H mfp')

ax.plot(Temp,p_length[0,:],color='red',label=r"ne = 0.4 x $10^{20} cm^{-3}$",linestyle='-')

ax.plot(Temp,p_length[1,:],color='blue',label=r"ne = 0.35 x $10^{20} cm^{-3}$",linestyle='-')

ax.set_xlabel("Temperature(KeV)")
ax.set_ylabel(r"$\sqrt{\frac{m_i}{m_e}}\lambda_{ij}$ (mm)")
ax.set_title(r"Pedestal length for $H^{+}-H^{+}$ plasma shock")

ax.legend()
plt.show()