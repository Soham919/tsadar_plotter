import numpy as np
import matplotlib.pyplot as plt

t = [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9]
tw4 = []
tw3_1d_1um_resolve = np.array([3320, 4100, 4830, 5360, 6090, 6620, 7160, 7710, 8240, 8610, 9140, 9680, 10060, 10600, 10960, 11500, 11860, 12400])/1000
tw3 = np.array([3340, 4100, 4780, 5400, 6000, 6550, 7070, 7600, 8090, 8600, 9070, 9550, 10020, 10460, 10890, 11340, 11750, 12170])/1000
tw2 = np.array([3290, 3960, 4570, 5140, 5680, 6200, 6690, 7140, 7620, 8070, 8500, 8930, 9340, 9750, 10150, 10540, 10950, 11310])/1000
tw1 = np.array([3140, 3680, 4190, 4650, 5090, 5480, 5880, 6270, 6650, 7000, 7370, 7690, 8030, 8360, 8680, 8980, 9300, 9600])/1000
tw05 = np.array([3070, 3310, 3680, 3960, 4230, 4510, 4750, 5000, 5220, 5460, 5690, 5900, 6130, 6330, 6540, 6740, 6940, 7130])/1000
u3 = np.gradient(tw3, t)
u2 = np.gradient(tw2, t)
u1 = np.gradient(tw1, t)
u05 = np.gradient(tw05, t)

# 2D 3TW data (1um resolve) for comparison
t2D = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]
e3 = np.array([200, 220, 240, 270, 340, 490, 630, 790, 940, 1090, 1220, 1840, 2350, 2820, 3230, 3620, np.nan, np.nan, np.nan]) 
e2 = np.array([200, 205, 218, 245, 295, 389, 502, 608, 702, 800, 890, 1230, 1520, 1770, 1990, 2190, 2370, 2550, 2730])
e1 = np.array([200, 210, 230, 240, 330, 400, 450, 510, 560, 600, 650, 830, 980, 1120, 1240, 1370, 1480, 1580, 1700 ]) 
compare_1d_2d = np.array([200, 220, 240, 320, 410, 490, 670, 870, 950, 1140, 1330, 1990, 2630, 3210, 3860, 4420, 4890, 5460, 6030])
ue3 = np.gradient(e3, t2D)
ue2 = np.gradient(e2, t2D)
ue1 = np.gradient(e1, t2D)
u1d2d = np.gradient(compare_1d_2d, t2D)


fig, ax = plt.subplots(figsize=(8, 5))
# ax.plot(t, tw05, marker="o", label="0.5 TW")
# ax.plot(t, tw1, marker="o", label="1 TW")
# ax.plot(t, tw2, marker="o", label="2 TW")
# ax.plot(t, tw3, marker="o", label="3 TW")
# ax.plot(t, tw3_1d_1um_resolve, marker="o", label="3 TW - 1um resolve")   

# ax.plot(t, u05, marker="o", label="0.5 TW")
# ax.plot(t, u1, marker="o", label="1 TW")
# ax.plot(t, u2, marker="o", label="2 TW")
# ax.plot(t, u3, marker="o", label="3 TW")   

# 2D data
# ax.plot(t2D, e1, marker="o", label="10um 3TW - 2D")
# ax.plot(t2D, e2, marker="o", label="100um 3TW - 2D")
# ax.plot(t2D, e3, marker="o", label="1000um 3TW - 2D")
# ax.plot(t2D, compare_1d_2d, marker="o", label="3TW - 1D")

ax.plot(t2D, ue1, marker="o", label="10um 3TW - 2D")
ax.plot(t2D, ue2, marker="o", label="100um 3TW - 2D")
ax.plot(t2D, ue3, marker="o", label="1000um 3TW - 2D")
ax.plot(t2D, u1d2d, marker="o", label="3TW - 1D")

ax.legend()
ax.set_xlabel("Time (ns)")
ax.set_ylabel("Shock speed (mm)")
plt.title("Shock front speed vs time for different laser spot sizes")
ax.grid(True, alpha=0.3)
plt.show()