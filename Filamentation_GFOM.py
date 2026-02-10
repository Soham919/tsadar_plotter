"""
Resources 

IAW damping ratio:
https://www.lle.rochester.edu/media/publications/lle_review/documents/v74/4_74accurate.pdf 

"""

import numpy as np
import astropy.units as u
import astropy.constants as const
import matplotlib.pyplot as plt

from plasmapy.formulary.densities import critical_density
from plasmapy.formulary.collisions import MaxwellianCollisionFrequencies
from plasmapy.formulary.speeds import thermal_speed
from plasmapy.particles import Particle
import warnings


def iaw_ratio_eq31_32(T_e=None, T_i=None, Z=None, mu=None):
    """
    The ratio of the ion acoustic wave frequency to the damping rate, or 
    
    \omega_i  / \omega_r 
    
    At the point on the dispersion curve relevant to the provided plasma
    
    https://www.lle.rochester.edu/media/publications/lle_review/documents/v74/4_74accurate.pdf
    """
    m_e = const.m_e.si
    m_i = const.m_p.si * mu
    

    
    omega_r = np.sqrt(1 + 3*T_i/Z/T_e)
    
    omega_i = np.sqrt(np.pi/8) * (
                    np.sqrt(Z*m_e/m_i) + 
                    (Z*T_e/T_i)**(3/2)*
                    np.exp( 
                        -Z*T_e/2/T_i - 3/2
                        )
                                )
    return (omega_i/omega_r).to(u.dimensionless_unscaled).value



def iaw_ratio(T_e=None, T_i=None, Z=None, mu=None):
    """
    The ratio of the ion acoustic wave frequency to the damping rate, or 
    
    \omega_i  / \omega_r 
    
    At the point on the dispersion curve relevant to the provided plasma
    
    https://www.lle.rochester.edu/media/publications/lle_review/documents/v74/4_74accurate.pdf
    """
    m_e = const.m_e.si
    m_i = const.m_p.si * mu
    
    
    T= (T_i/Z/T_e).to(u.dimensionless_unscaled)
    M= (Z*m_e/m_i).to(u.dimensionless_unscaled)
    Γ = np.sqrt(np.pi/8)*np.sqrt(M)/T**2
    Δ = np.sqrt(np.pi/8)*T**(-7/2)*np.exp(-1/2/T - 3/2)
    
    if T>0.1:
        warnings.warn(f"T={T}, but this model assumes T<<1")
    
    

    
    
    Ωr = 1 + 3*T/2 + 15*T**2/8 + (147/16 +  Δ*(Γ+Δ) )*T**3
    
    Ωi = np.sqrt(np.pi/8) * (
                        np.sqrt(M) +
                        T**(3/2)*
                        np.exp(-1/2/T - 3/2 - 3*T)
                        )
    
    return (Ωi/Ωr).to(u.dimensionless_unscaled).value


def gamma_t(ion=None, rho0=None, n_e=None, T_e=None, T_i=None):
    
    col = MaxwellianCollisionFrequencies(test_particle='e-',
                                         field_particle=ion,
                                         T_a=T_e,
                                         T_b=T_i,
                                         n_a=n_e,
                                         n_b=n_e,
                                         Coulomb_log=10*u.dimensionless_unscaled,
                                         )
    nu_ei = col.Maxwellian_avg_ei_collision_freq
    vte = thermal_speed(T=T_e, particle='e-')
    mfp = (vte/nu_ei).to(u.um)
    Z = ion.charge_number
    
    return 1 + 1.76*Z**(5/7)*(rho0/mfp)**(4/7)



log_ne = np.linspace(18, 20.5, num=100)
ne = 10**log_ne * u.cm**-3
T_i = 1*u.eV
temperatures = [10*u.eV, 25*u.eV, 50*u.eV, 100*u.eV]
ion = Particle('p+')

probe_wavelengths = np.array([526.5, 263.25])*u.nm
probe_energy= 45*u.J
probe_pulse = 3.7*u.ns
probe_fnum = 6.7 # Hansen 2019 mitigation
probe_area = np.pi*(100*u.um/2)**2
probe_intensity = (probe_energy/probe_pulse/probe_area).to(u.W/u.cm**2)
print(probe_intensity/(10**(14)))




"""
# Turnbull s101402 
ne = 3.91e20*u.cm**-3
T_i = 10*u.eV
temperatures = np.array([1130])*u.eV
ion = Particle('N-14 4+')

probe_wavelengths = np.array([526.5, 263.25])*u.nm
probe_energy= 30*u.J
probe_pulse = 3.7*u.ns
probe_fnum = 6.7 # Hansen 2019 mitigation
probe_area = np.pi*(100*u.um/2)**2
probe_intensity = 7.4e14 * u.W/u.cm**2
"""





fig, ax = plt.subplots()
ax.set_xlabel("Density (cm$^{-3}$)")
ax.set_xscale('log')
ax.set_ylim(0,5)
ax.axhline(1, color='gray')

colors = ['lime', 'blueviolet']
linestyles = ['solid', 'dashed','dotted', 'dashdot']

ffom = 'Pg.1'

for w, wavelength in enumerate(probe_wavelengths):
    probe_freq = (2*np.pi*u.rad*const.c.si/wavelength).to(u.rad/u.s) 
    nc = critical_density(probe_freq)
    
    transverse_speckle_width = (2/np.pi)*probe_fnum*wavelength

    for t,Te in enumerate(temperatures):
        
        
        
        
        if ffom == 'Pg.1':
            ax.set_ylabel("Filamentation FOM ")
            #Pg. 1 FFOM
            FFOM = (10*
                   (probe_intensity).to(u.W/u.cm**2).value*1e-14*
                   ((wavelength).to(u.um).value)**2*
                   (ne/nc).to(u.dimensionless_unscaled).value*
                   3/Te.to(u.keV).value*
                   (probe_fnum/8)**2 )
        elif ffom == 'Eq.1':
            
            ax.set_ylabel("Filamentation FOM (Grech) Time resolved TS ")
            
            ir = iaw_ratio(T_e=Te, T_i=T_i, Z=ion.charge_number, 
                           mu=ion.mass_number)
            gt = gamma_t(ion=ion, rho0=transverse_speckle_width,
                         T_e=Te, T_i=T_i, n_e=ne,
                         )
            
            # Note: ir is defined per Table 1 in Turnbull's PRL, but
            # the ratio is inverted in the actual FFOM
            
            FFOM = (gt*(1/ir)*
                    (probe_intensity).to(u.W/u.cm**2).value*1e-14*
                    ((wavelength).to(u.um).value)**2*
                    (ne/nc).to(u.dimensionless_unscaled).value*
                    3/Te.to(u.keV).value*
                    (probe_fnum/8)**2 
                    )
                    
        
        ax.plot(ne.value, FFOM, linestyle=linestyles[t],
                color=colors[w])
        
        
for w, wavelength in enumerate(probe_wavelengths):
    ax.plot([],[], color=colors[w], linestyle='solid',  label=f"{wavelength.value} nm")
    
for t,Te in enumerate(temperatures):
    ax.plot([],[], color='k', linestyle=linestyles[t], label=f"{Te.value} eV")
    
ax.legend(loc='upper right')
plt.show()