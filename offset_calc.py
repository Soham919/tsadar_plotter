import numpy as np
import matplotlib.pyplot as plt

e_disp = 0.27093
i_disp = 0.0057

ep_px = 565  #electron pointing pixel
ecw_px = 523  #electron cw pixel

ip_px = 625  #ion pointing pixel
icw_px = 493  #ion cw pixel

epw_offset = 526.5 - (ecw_px*e_disp)
iaw_offset = 526.5 - (icw_px*i_disp)

print(f"epw offset = {epw_offset}")
print(f"iaw_offset = {iaw_offset}")

x = np.linspace(0,10,10)



