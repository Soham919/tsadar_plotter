import numpy as np
import matplotlib.pyplot as plt

sg1 = np.loadtxt(
    "//profiles/Users$/sban/Documents/Kinshock/Kinshock-27A/Re__1_ns_pulse/sg10v001.txt",
    delimiter=","
)

sg_600_t = np.array([0, 0.177, 0.634, 0.818])
sg_600_P = np.array([0, 0.5, 0.5, 0])

sg_1_t = np.array([0, 0.177, 1.00, 1.177])
sg_1_P = np.array([0, 0.5, 0.5, 0])

# Integrate power over time
E_sg600 = np.trapezoid(sg_600_P, sg_600_t)
E_sg1   = np.trapezoid(sg_1_P, sg_1_t)

# Loaded sg1 curve (with same 0.461 ns time shift as plot)
E_sg1_loaded = np.trapezoid(sg1[:, 1], sg1[:, 0] - 0.461)

print("600 ps pulse integral:", E_sg600, "TW ns")
print("1 ns pulse integral:", E_sg1, "TW ns")
print("Loaded pulse integral:", E_sg1_loaded, "TW ns")

plt.plot(sg1[:, 0] - 0.461, sg1[:, 1], label='1 ns')
plt.plot(sg_600_t, sg_600_P, label='600 ps')
plt.plot(sg_1_t, sg_1_P, label='1 ns', color='green')

plt.xlabel('Time (ns)')
plt.ylabel('Power (TW)')
plt.legend()
plt.show()