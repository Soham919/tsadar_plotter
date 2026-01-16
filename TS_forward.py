import numpy as np
import matplotlib.pyplot as plt
from Plasma_func import mfp_ii

l = np.linspace(0.01,0.99,100)
x = 1/(1-l)
y = 10/l

z = (1/x + 1/y)**(-1)
plt.plot(l,z)
plt.plot(l,x)
plt.plot(l,y)
plt.show()