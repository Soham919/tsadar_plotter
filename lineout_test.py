import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append(r"\\profiles\Users$\sban\Documents\Scripts")
from lineout_tool import *

x = np.linspace(0, 10, 100)
y = np.linspace(0, 10, 100)

def f(x, y):
    return np.sin(x) * np.cos(y)

z = f(x[:, None], y[None, :])

fig, ax  = plt.subplots(figsize=(8, 6))
ax.pcolormesh(x, y, z, shading='auto', cmap='viridis')
#ax.colorbar(label='f(x, y) = sin(x) * cos(y)') 
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('2D Color Plot of f(x, y)')


lineout_tool = InteractiveLineout(
    ax=ax,
    x=x,
    y=y,
    data=z,
    width=0.1,
    lineout_label="random function",
)

plt.show()
