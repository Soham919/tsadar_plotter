import numpy as np
import matplotlib.pyplot as plt
import sys
from pathlib import Path
sys.path.append(r"\\profiles\Users$\sban\Documents\Scripts")
from fitter_1D import fit_1d_curve


# t = [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9]
# tw4 = []
# tw3_1d_1um_resolve = np.array([3320, 4100, 4830, 5360, 6090, 6620, 7160, 7710, 8240, 8610, 9140, 9680, 10060, 10600, 10960, 11500, 11860, 12400])/1000
# tw3 = np.array([3340, 4100, 4780, 5400, 6000, 6550, 7070, 7600, 8090, 8600, 9070, 9550, 10020, 10460, 10890, 11340, 11750, 12170])/1000
# tw2 = np.array([3290, 3960, 4570, 5140, 5680, 6200, 6690, 7140, 7620, 8070, 8500, 8930, 9340, 9750, 10150, 10540, 10950, 11310])/1000
# tw1 = np.array([3140, 3680, 4190, 4650, 5090, 5480, 5880, 6270, 6650, 7000, 7370, 7690, 8030, 8360, 8680, 8980, 9300, 9600])/1000
# tw05 = np.array([3070, 3310, 3680, 3960, 4230, 4510, 4750, 5000, 5220, 5460, 5690, 5900, 6130, 6330, 6540, 6740, 6940, 7130])/1000
# u3 = np.gradient(tw3, t)
# u2 = np.gradient(tw2, t)
# u1 = np.gradient(tw1, t)
# u05 = np.gradient(tw05, t)

# ------------------------------ 2D 3TW data (1um resolve) for comparison ------------------------------------------------------ #
t2D = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]
e3 = np.array([200, 220, 240, 270, 340, 490, 630, 790, 940, 1090, 1220, 1840, 2350, 2820, 3230, 3620, np.nan, np.nan, np.nan]) 
e2 = np.array([200, 205, 218, 245, 295, 389, 502, 608, 702, 800, 890, 1230, 1520, 1770, 1990, 2190, 2370, 2550, 2730])
e1 = np.array([200, 210, 230, 240, 330, 400, 450, 510, 560, 600, 650, 830, 980, 1120, 1240, 1370, 1480, 1580, 1700 ]) 
compare_1d_2d = np.array([200, 220, 240, 320, 410, 490, 670, 870, 950, 1140, 1330, 1990, 2630, 3210, 3860, 4420, 4890, 5460, 6030])
ue3 = np.gradient(e3, t2D)
ue2 = np.gradient(e2, t2D)
ue1 = np.gradient(e1, t2D)
u1d2d = np.gradient(compare_1d_2d, t2D)
# ------------------------------------------------------------------------------------------------------------------------------ #
# ------------------------------- Time deresolved 2D 3TW data (1um resolve) ---------------------------------------------------- #
# t2D_td = [0, 0.2, 0.4, 0.6, 0.8, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]
# e3_td = np.array([200, 250, 370, 650, 960, 1240, 1850, 2380, 2840, 3250, 3640, 3990, np.nan, np.nan]) 
# e2_td = np.array([200, 250, 320, 460, 740, 910, 1250, 1540, 1730, 2020, 2240, 2440, 2630, 2810])
# e1_td = np.array([200, 250, 330, 450, 550, 630, 810, 950, 1100, 1240, 1370, 1490, 1610, 1720]) 

# ue3_td = np.gradient(e3_td, t2D_td)
# ue2_td = np.gradient(e2_td, t2D_td)
# ue1_td = np.gradient(e1_td, t2D_td)
# ------------------------------------------------------------------------------------------------------------------------------ #

# --------------------------------- 2D offset runs ----------------------------------------------------------------------------- #
t2D_offset = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9]

e3_8mm_offset = np.array([200, 240, 270, 350, 500, 770, 990, 1170, 1380, 1560, 1740, 2500, 3090, 3610, 4070, 4480, 4820, 5180, 5490, 5780, 6060, 6350, 6610, 6880, 7150, 7400, 7650])
e3_6_5mm_offset = np.array([200, 240, 270, 380, 570, 730, 890, 1060, 1200, 1360, 1500, 2090, 2590, 3050, 3430, 3800, 4120, 4440, 4740, 5030, 5310, 5580, 5850, 6100, 6370, 6600, 6850])
e3_5mm_offset = np.array([200, 250, 260, 320, 460, 610, 770, 920, 1090, 1240, 1380, 1980, 2460, 2910, 3320, 3690, 4030, 4350, 4670, 4990, 5260, 5560, 5830, 6100, np.nan, np.nan, np.nan])

ue3_8mm_offset = np.gradient(e3_8mm_offset, t2D_offset)
ue3_6_5mm_offset = np.gradient(e3_6_5mm_offset, t2D_offset)
ue3_5mm_offset = np.gradient(e3_5mm_offset, t2D_offset)
# ------------------------------------------------------------------------------------------------------------------------------ #
# --------------------------------- 2D offset density variation runs ----------------------------------------------------------- #
t2D_offset = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9]
fp = Path(r"C:\Simulation_data\FLASH\2D\offset\dens_vary")
half_dens_8mm = np.load(fp / "1mmSpot_0.5dens_8mm_offset" / "shock_positions.npy")
half_dens_8mm[5] = 750
half_dens_8mm[6] = 1030
half_dens_8mm[7] = 1240
half_dens_8mm[8] = 1450
half_dens_8mm[9] = 1650
half_dens_8mm[10] = 1830
quart_dens_8mm = np.load(fp / "1mmSpot_0.25dens_8mm_offset" / "shock_positions.npy")
quart_dens_8mm[5] = 720
quart_dens_8mm[6] = 1030
quart_dens_8mm[7] = 1310
quart_dens_8mm[8] = 1500
quart_dens_8mm[9] = 1740
quart_dens_8mm[10] = 1950
half_dens_6_5mm = np.load(fp / "1mmSpot_0.5dens_6.5mm_offset" / "shock_positions.npy")
half_dens_6_5mm[-1] = np.nan
half_dens_6_5mm[-2] = np.nan
half_dens_6_5mm[-3] = np.nan
quart_dens_6_5mm = np.load(fp / "1mmSpot_0.25dens_6.5mm_offset" / "shock_positions.npy")

u1 = np.gradient(half_dens_8mm, t2D_offset)
u2 = np.gradient(quart_dens_8mm, t2D_offset)
u3 = np.gradient(half_dens_6_5mm, t2D_offset)
u4 = np.gradient(quart_dens_6_5mm, t2D_offset)
# ------------------------------------------------------------------------------------------------------------------------------ #
# --------------------------------- 2D spot sizes same power ----------------------------------------------------------- #
t2D_spot = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9]
fp = Path(r"C:\Simulation_data\FLASH\2D\Spot_size_vary")
#half_dens_8mm = np.load(fp / "100um_spot_Si" / "shock_positions.npy")
e2_power = np.load(fp / "100umSpot_samePower" / "shock_positions.npy")
u2_power = np.gradient(e2_power, t2D_spot)
e1_power = np.load(fp / "10umSpot_samePower" / "shock_positions.npy")
u1_power = np.gradient(e1_power, t2D_spot)

def sedov_velocity(t, C):
    """
    Sedov-Taylor blast wave solution for shock speed as a function of time.

    Parameters
    ----------
    t : array-like
        Time array (seconds).
    B : float
        Coefficient B.
    C : float
        Coefficient C.

    Returns
    -------
    R : array-like
        Shock speed as a function of time (meters/second).
    """
    return (t/C)**(-1/2)

fit = fit_1d_curve(t2D_offset[10:-1], ue3_5mm_offset[10:-1], sedov_velocity, plot=True, xlabel="Time (ns)", ylabel="Shock speed (mm/ns)", title="Sedov-Taylor Fit for 3TW 5mm offset", data_label="3TW 5mm offset", fit_label="Sedov-Taylor Fit")
# ------------------------------------------------------------------------------------------------------------------------------ #

#fig, ax = plt.subplots(figsize=(8, 5))

# ------------ 1D different powers ----------------- #
# ax.plot(t, tw05, marker="o", label="0.5 TW")
# ax.plot(t, tw1, marker="o", label="1 TW")
# ax.plot(t, tw2, marker="o", label="2 TW")
# ax.plot(t, tw3, marker="o", label="3 TW")
# ax.plot(t, tw3_1d_1um_resolve, marker="o", label="3 TW - 1um resolve")   

# ax.plot(t, u05, marker="o", label="0.5 TW")
# ax.plot(t, u1, marker="o", label="1 TW")
# ax.plot(t, u2, marker="o", label="2 TW")
# ax.plot(t, u3, marker="o", label="3 TW")   

# --------- 2D data spot sizes same intensity ---------------- #
# ax.plot(t2D, e1, marker="o", label="10um 3TW")
# ax.plot(t2D, e2, marker="o", label="100um 3TW")
# ax.plot(t2D, e3, marker="o", label="1000um 3TW")
# ax.plot(t2D, compare_1d_2d, marker="o", label="1D - 3TW")

# ax.plot(t2D, ue1, marker="o", label="10um 3TW - 2D")
# ax.plot(t2D, ue2, marker="o", label="100um 3TW - 2D")
# ax.plot(t2D, ue3, marker="o", label="1000um 3TW - 2D")
# ax.plot(t2D, u1d2d, marker="o", label="3TW - 1D")

# --------- 2D data spot sizes same power---------------- #
# ax.plot(t2D, e2, marker="o", label="100um same Intensity")
# ax.plot(t2D_spot, e2_power, marker="o", label="100um same Power")
#ax.plot(t2D_spot, e1_power, marker="o", label="10um same Power")

#ax.plot(t2D, ue3, marker="o", label="1000um same power - 2D")
#ax.plot(t2D_spot, u2_power, marker="o", label="100um same Power - 2D")
# ax.plot(t2D_spot, u1_power, marker="o", label="10um same Power - 2D")


# ax.plot(t2D, ue1, marker="o", label="10um 3TW - 2D")
# ax.plot(t2D, ue2, marker="o", label="100um 3TW - 2D")
# ax.plot(t2D, ue3, marker="o", label="1000um 3TW - 2D")
# ax.plot(t2D, u1d2d, marker="o", label="3TW - 1D")

# --------- 2D time deresolved data spot sizes ---------------- #
# ax.plot(t2D_td, e1_td, marker="o", label="10um time deresolved")
# ax.plot(t2D_td, e2_td, marker="o", label="100um time deresolved")
# ax.plot(t2D_td, e3_td, marker="o", label="1000um time deresolved")




# ------------ 2D offset data ------------------------- #
# ax.plot(t2D_offset, e3_8mm_offset, marker="o", label="8mm offset 2D")
# ax.plot(t2D_offset, e3_6_5mm_offset, marker="o", label="6.5mm offset 2D")
# ax.plot(t2D_offset, e3_5mm_offset, marker="o", label="5mm offset 2D")

# ax.plot(t2D_offset, ue3_8mm_offset, marker="o", label="8mm offset 2D")
# ax.plot(t2D_offset, ue3_6_5mm_offset, marker="o", label="6.5mm offset 2D")
# ax.plot(t2D_offset, ue3_5mm_offset, marker="o", label="5mm offset 2D")

# ------------- 2D offset density variation data ------------------------- #

# ax.plot(t2D_offset, e3_8mm_offset, color='blue', marker="o", markerfacecolor='white', markeredgecolor='blue', markersize=7, label=r"8mm offset $\rho$")
# ax.plot(t2D_offset, e3_6_5mm_offset, color='blue', marker="^", markerfacecolor='white', markeredgecolor='blue', markersize=7, label=r"6.5mm $\rho$")
# ax.plot(t2D_offset, half_dens_6_5mm, color='red', marker="^", markerfacecolor='white', markeredgecolor='red', markersize=7, label=r"6.5mm offset $\rho/2$")
# ax.plot(t2D_offset, quart_dens_6_5mm, color='green', marker="^", markerfacecolor='white', markeredgecolor='green', markersize=7, label=r"6.5mm offset $\rho/4$")   
# ax.plot(t2D_offset, half_dens_8mm, color='red', marker="o", markerfacecolor='white', markeredgecolor='red', markersize=7, label=r"8mm offset $\rho/2$")
# ax.plot(t2D_offset, quart_dens_8mm, color='green', marker="o", markerfacecolor='white', markeredgecolor='green', markersize=7, label=r"8mm offset $\rho/4$")

# ax.plot(t2D_offset, ue3_8mm_offset, color='blue', marker="o", markerfacecolor='white', markeredgecolor='blue', markersize=7, label=r"8mm offset $\rho$")
# ax.plot(t2D_offset, ue3_6_5mm_offset, color='blue', marker="^", markerfacecolor='white', markeredgecolor='blue', markersize=7, label=r"6.5mm $\rho$")
# ax.plot(t2D_offset, u3, color='red', marker="^", markerfacecolor='white', markeredgecolor='red', markersize=7, label=r"6.5mm offset $\rho/2$")
# ax.plot(t2D_offset, u4, color='green', marker="^", markerfacecolor='white', markeredgecolor='green', markersize=7, label=r"6.5mm offset $\rho/4$")   
# ax.plot(t2D_offset, u1, color='red', marker="o", markerfacecolor='white', markeredgecolor='red', markersize=7, label=r"8mm offset $\rho/2$")
# ax.plot(t2D_offset, u2, color='green', marker="o", markerfacecolor='white', markeredgecolor='green', markersize=7, label=r"8mm offset $\rho/4$")

# ax.legend()
# ax.set_xlabel("Time (ns)")
# ax.set_ylabel("Shock speed (mm/ns)")
# #ax.set_ylim(0, 3500)
# plt.title("Shock speed vs time for different laser spot sizes")
# ax.grid(True, alpha=0.3)
# plt.show()