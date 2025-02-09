#%%
import matplotlib.pyplot as plt
plt.rc("figure", figsize=[12,8])
plt.rc("font", size=8)
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


#%%
lko = ksp.calc_orbit(Kerbin,e=0,T=1840)
lko.max_alt,lko.min_alt
# %%
lko.do_maneuver(0,10)
# %%
