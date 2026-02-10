import numpy
import matplotlib.pyplot as plt


def GJS_target_disp(P,A,t,m):
    F = P*A
    a = F/m
    d = 0.5*a*(t**2)
    return d

P = 2e5  # Pressure in Pa
A1 = 4e-6   # Area of the Si3N4 window in m^2
l1 = 10**(-6) # thickness of the Si3N4 window in m
rho1 = (3.2)*(10**(3)) # mass density of Si3N4 in kg/m^3
m1 = rho1*A1*l1 # mass of window in kg
A2 = (25-4)*(10**(-6))   # Area of the Si frame in m^2
l2 = 200*(10**(-6)) # thickness of the Si frame in m
rho2 = (2.33)*(10**(3)) # mass density of Si in kg/m^3
m2 = A2*l2*rho2 # mass of frame in kg
m = m1 + m2 # ~ total mass of target in kg

print(m)
t = (160)*(10**(-6)) # time in seconds

print(f"Target is displaced by {GJS_target_disp(P,(A1+A2),t,m)*1000} mm")