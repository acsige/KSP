#%%
import matplotlib.pyplot as plt
plt.rc("figure", figsize=[12,8])
plt.rc("font", size=8)
from math import pi
import KSP_module as ksp
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

    LKO = ksp.orbit(Kerbin, min_alt = 70000.1, e=0)
    transfer3 = ksp.calc_window(LKO, Mun.orbit, 0)
    fig,ax = ksp.plot_hohmann_orbit(LKO, Mun, transfer3)
    plt.show()
    print([ksp.pretty_date(transfer3.t_launch), transfer3.t_launch])
# %%
LKOe = ksp.orbit(Kerbin, min_alt = 70000.1, e=0.5)
t1,d1 = LKOe.calc_min_distance_to(Kerbin, 0, LKOe.T/2)
t2,d2 = LKOe.calc_min_distance_to(Kerbin, LKOe.T/2, LKOe.T)
print(t1,d1)
fig,ax = ksp.initialize_plot()
ax = ksp.add_orbit_to_plot(ax, LKOe)
ax = ksp.add_orbit_point_to_plot(ax, LKOe, 0, label='t0', marker='X')
ax = ksp.add_orbit_point_to_plot(ax, LKOe, t1, label='t1', marker='o')



