import numpy as np
import matplotlib.pyplot as plt

t = [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9]
tw4 = []
tw3 = np.array([3340, 4100, 4780, 5400, 6000, 6550, 7070, 7600, 8090, 8600, 9070, 9550, 10020, 10460, 10890, 11340, 11750, 12170])/1000
tw2 = np.array([3290, 3960, 4570, 5140, 5680, 6200, 6690, 7140, 7620, 8070, 8500, 8930, 9340, 9750, 10150, 10540, 10950, 11310])/1000
tw1 = np.array([3140, 3680, 4190, 4650, 5090, 5480, 5880, 6270, 6650, 7000, 7370, 7690, 8030, 8360, 8680, 8980, 9300, 9600])/1000
tw05 = np.array([3070, 3310, 3680, 3960, 4230, 4510, 4750, 5000, 5220, 5460, 5690, 5900, 6130, 6330, 6540, 6740, 6940, 7130])/1000
u3 = np.gradient(tw3, t)
u2 = np.gradient(tw2, t)
u1 = np.gradient(tw1, t)
u05 = np.gradient(tw05, t)

fig, ax = plt.subplots(figsize=(8, 5))
# ax.plot(t, tw05, marker="o", label="0.5 TW")
# ax.plot(t, tw1, marker="o", label="1 TW")
# ax.plot(t, tw2, marker="o", label="2 TW")
# ax.plot(t, tw3, marker="o", label="3 TW")   

ax.plot(t, u05, marker="o", label="0.5 TW")
ax.plot(t, u1, marker="o", label="1 TW")
ax.plot(t, u2, marker="o", label="2 TW")
ax.plot(t, u3, marker="o", label="3 TW")   

ax.legend()
ax.set_xlabel("Time (ns)")
ax.set_ylabel("Shock front speed (mm/ns)")
plt.title("Shock front speed vs time for different laser powers")
ax.grid(True, alpha=0.3)
plt.show()