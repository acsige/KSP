#%%
import matplotlib.pyplot as plt
plt.rc("figure", figsize=[12,8])
plt.rc("font", size=8)
from math import pi, sqrt
import KSP_module as ksp
from KSP_module import Kerbol, Kerbin, Duna, Eve, Mun, Minmus, Moho
#%%
# Testing during development
if __name__ == "__main__":
    # Plotting Hohmann transfer orbits
    # transfer1 = ksp.calc_window(Kerbin.orbit, Duna.orbit, 0)
    # fig,ax = ksp.plot_hohmann_orbit(Kerbin, Duna, transfer1)
    # plt.show()
    # print([ksp.pretty_date(transfer1.t_launch), transfer1.t_launch])
    
    # transfer2 = ksp.calc_window(Kerbin.orbit, Eve.orbit, 0)
    # fig,ax = ksp.plot_hohmann_orbit(Kerbin, Eve, transfer2)
    # plt.show()
    # print([ksp.pretty_date(transfer2.t_launch), transfer2.t_launch])

    LKO = ksp.orbit(Kerbin, min_alt = 70000.1, e=0)
    transfer3 = ksp.calc_window(LKO, Mun.orbit, 0)
    fig,ax = ksp.plot_hohmann_orbit(LKO, Mun, transfer3)
    ax = ksp.add_soi_to_plot(ax, Mun, transfer3.t_launch)
    plt.show()
    print([ksp.pretty_date(transfer3.t_launch), transfer3.t_launch])
# %%
LKOe = ksp.orbit(Kerbin, min_alt = 70000.1, e=0.2, nu0=2)
t1,d1 = LKOe.calc_min_distance_to(Kerbin, LKOe.T/2, LKOe.T)
print(t1,d1)
fig,ax = ksp.initialize_plot()
ax = ksp.add_orbit_to_plot(ax, LKOe)
ax = ksp.add_orbit_point_to_plot(ax, LKOe, LKOe.t0, label='t0', marker='X')
ax = ksp.add_orbit_point_to_plot(ax, LKOe, t1, label='t1', marker='X')
# ax = ksp.add_orbit_point_to_plot(ax, LKOe, LKOe.T/2, label='T/2', marker='X')
# ax = ksp.add_orbit_point_to_plot(ax, LKOe, LKOe.T, label='T', marker='X')
# ax = ksp.add_orbit_point_to_plot(ax, LKOe, t1, label='t1', marker='o')

# %%
LKO = ksp.orbit(Kerbin, min_alt = 70000.1, e=0)
transfer = ksp.calc_window(Kerbin.orbit, Moho.orbit, 0)
fig,ax = ksp.plot_hohmann_orbit(Kerbin, Moho, transfer)
ax = ksp.add_soi_to_plot(ax, Moho, transfer.t_launch)
plt.show()

min_dist = transfer.calc_distance_to(Moho, transfer.t_arrival)
# velocity relative to Kerbin when coming out of Kerbin SOI
v_soi = abs(transfer.calc_speed(transfer.t_launch) - Kerbin.orbit.calc_speed(transfer.t_launch))
# velocity needed at Kerbin periapsis to reach SOI with the speed necessary
v_lko_launch = sqrt(v_soi**2 + 2*Kerbin.GM*(1/LKO.ra - 1/Kerbin.soi))
dv = v_lko_launch - LKO.va
print(v_soi, v_lko_launch, dv)

# %%
def calc_hohmann_dv(hohmann_orbit, parking_orbit):
    """dv calculation: from parking orbit around a body to Hohmann orbit out of its SOI"""
    start_body = parking_orbit.primary
    # velocity change needed for the transfer orbit (relative to start body)
    t0 = hohmann_orbit.t_launch
    v_soi = abs(hohmann_orbit.calc_speed(t0) - start_body.orbit.calc_speed(t0))
    # velocity needed at parking orbit to reach SOI with the speed necessary
    v_parking_launch = sqrt(v_soi**2 + 2*start_body.GM*(1/parking_orbit.a - 1/start_body.soi))
    # velocity change needed to reach the launch velocity
    dv = v_parking_launch - parking_orbit.va
    return dv


#%%
# Example usage
LKO = ksp.orbit(Kerbin, min_alt = 70000.1, e=0)

t0 = 0
for i in range(3):
    transfer = ksp.calc_window(Kerbin.orbit, Moho.orbit, t0)
    dv = calc_hohmann_dv(transfer, LKO)
    min_dist = transfer.calc_distance_to(Moho, transfer.t_arrival)
    print(f'Time of launch: {ksp.pretty_date(transfer.t_launch)}')
    print(f'Time of arrival: {ksp.pretty_date(transfer.t_arrival)}')
    print(f"Minimum distance to Moho at arrival: {min_dist/1e6:.2f} M km")
    print(f"Delta-v needed for transfer from LKO to Moho: {dv:.2f} m/s")
    print('-------------------')
    t0 = transfer.t_launch + 1000
    print(t0)

# %%
transfer = ksp.calc_window(Kerbin.orbit, Moho.orbit, 6717200)
# %%
