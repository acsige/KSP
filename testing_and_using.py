#%%
import matplotlib.pyplot as plt
plt.rc("figure", figsize=[6,4])
plt.rc("font", size=8)
from math import pi, sqrt
import KSP_module as ksp
from KSP_module import Kerbol, Kerbin, Duna, Eve, Mun, Minmus, Moho
#%%
class orbit_devel(ksp.orbit):
    def __init__(self, primary, **kwargs):
        super().__init__(primary, **kwargs)
        print("orbit_devel init")
#%%
# Testing during development
if __name__ == "__main__":
    LKO = orbit_devel(Kerbin, min_alt = 70000.1, e=0)
    transfer3 = ksp.calc_window(LKO, Mun.orbit, 0)
    transfer3_devel = orbit_devel(Kerbin, **transfer3.orbit_kwargs)
    t_soi = transfer3.calc_soi_change(0)[0]

    fig,ax = ksp.plot_hohmann_orbit(LKO, Mun, transfer3)
    ax = ksp.add_soi_to_plot(ax, Mun, t_soi)
    ax = ksp.add_planetary_body_to_plot(ax, Mun, t_soi, label=None)
    ax = ksp.add_orbit_point_to_plot(ax, transfer3, t_soi, label="SOI enter")
    plt.show()
    print([ksp.pretty_date(transfer3.t_launch), transfer3.t_launch])

# %%