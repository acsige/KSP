#%%
import numpy as np
import matplotlib.pyplot as plt
plt.rc("figure", figsize=[12,8])
plt.rc("font", size=8)
import KSP_module as ksp
from KSP_module import orbit, planetary_body, Kerbol
# %% 
# Testing during development
if __name__ == "__main__":
    #Moho = orbit(a=5263138304, e=0.2, omega=70)
    Eve = planetary_body('Eve', orbit(Kerbol, a=9832684544, e=0.01, omega=15, i=2.1),
                        GM=8.1717302e12, radius=7e5, atmo_height=9e4)
    Kerbin = planetary_body('Kerbin', orbit(Kerbol, a=13599840256, e=0),
                            GM=3.5316e12, radius=6e5, atmo_height=7e4)
    Duna = planetary_body('Duna', orbit(Kerbol, a=20726155264, e=0.051, omega=135.5, i=0.06),
                          GM=3.0136321e11, radius=3.2e5, atmo_height=5e4)
    Mun = planetary_body('Mun', orbit(Kerbin, a=1.2e6, e=0, nu0=1.7),
                         GM=6.5138398e10, radius=2e5, atmo_height=0)
    Minmus = planetary_body('Minmus', orbit(Kerbin, a=4.7e7, e=0, nu0=0.9),
                         GM=1.7658e9, radius=6e4, atmo_height=0)
# %% 
    
    ksp.plot_hohmann_orbits([Kerbin.orbit, Eve.orbit, ksp.calc_window(Kerbin.orbit, Eve.orbit, 0)])
    ksp.plot_hohmann_orbits([Kerbin.orbit, Duna.orbit, ksp.calc_window(Kerbin.orbit, Duna.orbit, 0)])

# %%
lko = ksp.calc_orbit(Kerbin,e=0,T=1840)
lko.max_alt,lko.min_alt
# %%
lko.do_maneuver(0,10)
# %%
