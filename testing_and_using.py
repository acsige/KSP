#%%
import matplotlib.pyplot as plt
plt.rc("figure", figsize=[12,8])
plt.rc("font", size=8)
from math import pi
import KSP_module as ksp
from KSP_module import orbit, planetary_body
from KSP_module import Kerbol, Kerbin, Duna, Eve, Mun, Minmus

# Testing during development
if __name__ == "__main__":
    # Plotting Hohmann transfer orbits
    transfer1 = ksp.calc_window(Kerbin.orbit, Duna.orbit, 0)
    fig,ax = ksp.plot_hohmann_orbit(Kerbin, Duna, transfer1)
    plt.show()
    print([ksp.pretty_date(transfer1.t_launch), transfer1.t_launch])
    
    transfer2 = ksp.calc_window(Kerbin.orbit, Eve.orbit, 0)
    fig,ax = ksp.plot_hohmann_orbit(Kerbin, Eve, transfer2)
    plt.show()
    print([ksp.pretty_date(transfer2.t_launch), transfer2.t_launch])

    LKO = orbit(Kerbin, min_alt = 70000.1, e=0)
    transfer3 = ksp.calc_window(LKO, Mun.orbit, 0)
    fig,ax = ksp.plot_hohmann_orbit(LKO, Mun, transfer3)
    plt.show()
    print([ksp.pretty_date(transfer3.t_launch), transfer3.t_launch])
# %%
LKO.do_maneuver(0,10)
    
# %%
