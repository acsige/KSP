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
LKO = ksp.orbit(Kerbin, min_alt = 70000.1, e=0)
t,d = LKO.calc_min_distance_to(Mun, 0, LKO.T)
fig,ax = ksp.initialize_plot()
ax = ksp.add_orbit_to_plot(ax, LKO)
ax = ksp.add_orbit_to_plot(ax, Mun.orbit)
ax = ksp.add_orbit_point_to_plot(ax, Mun.orbit, t, label='Mun at minimum distance', marker='o')
ax = ksp.add_orbit_point_to_plot(ax, LKO, t, label='SV at minimum distance', marker='o')
# %%
LKO = ksp.orbit(Kerbin, min_alt = 70000.1, e=0)
t,d = LKO.calc_min_distance_to(Mun, 0, LKO.T)
r1,p1 = Mun.orbit.calc_polar(t)
r2,p2 = LKO.calc_polar(t)
print(abs(p2-p1)*180/pi)
print(abs(r2-r1)-d)
# %%
