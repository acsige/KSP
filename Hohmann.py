#%%
import numpy as np
from math import sqrt, pi
import matplotlib.pyplot as plt
plt.rc("figure", figsize=[12,8])
plt.rc("font", size=8)

from KSP_module import orbit, planetary_body, Kerbol, calc_orbit

#%%
def calc_hohmann(src_orbit, dst_orbit, t0):
    """
    Calculate semi-major axis, eccentricity, average angular velocity and time for Hohmann transfer.
    """
    assert(src_orbit.primary == dst_orbit.primary)
    primary = src_orbit.primary
    GM = primary.GM

    # position of src known at t0
    phi, r_src = src_orbit.calc_polar(t0)
    
    # initial guess for dest is using 0 time for Hohmann transfer
    hohmann_time = 0
    hohmann_time_prev = 10

    # Iterate a better Hohmann time by recalculating destination
    while abs(hohmann_time_prev-hohmann_time) > 1:
        phi, r_dst = dst_orbit.calc_polar(t0+hohmann_time)
        # semi-major axis
        a = (r_src + r_dst) / 2
        # Time of Hohmann transfer is half of the Hohmann orbit
        r_avg = 2*a/pi
        hohmann_period_sq = 4*pi**2*a**3/GM
        hohmann_time_prev = hohmann_time
        hohmann_time = sqrt(hohmann_period_sq) / 2    
        
    # eccentricity of Hohmann orbit
    if r_src > r_dst:
        e = 1 - 2/(r_src/r_dst + 1)
    else:
        e = 1 - 2/(r_dst/r_src + 1)

    # average angular velocity
    n = pi/hohmann_time
    
    return a,e,n,hohmann_time

def plot_orbits(orbit_list):
    """Plot planet orbits, 0: source, 1: destination, 2: transfer orbit """
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.clear()
    for p in orbit_list:
        ax.plot(p.phi, p.r)
    t_launch = orbit_list[2].t_launch
    t_arrival = t_launch + orbit_list[2].T/2

    ax.plot(orbit_list[0].calc_polar(t_launch)[0],orbit_list[0].calc_polar(t_launch)[1], 'o', label="Launch")
    ax.plot(orbit_list[1].calc_polar(t_launch)[0],orbit_list[1].calc_polar(t_launch)[1], 'o', label="Dest at launch")
    ax.plot(orbit_list[1].calc_polar(t_arrival)[0],orbit_list[1].calc_polar(t_arrival)[1], 'o', label="Arrival")
    ax.plot(orbit_list[2].calc_polar(t_launch)[0],orbit_list[2].calc_polar(t_launch)[1], 'x', label="Hohmann launch")
    ax.plot(orbit_list[2].calc_polar(t_arrival)[0],orbit_list[2].calc_polar(t_arrival)[1], 'x', label="Hohmann arrival")
    
    ax.set_rticks([0.5, 1, 1.5, 2])  # Less radial ticks
    ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
    ax.grid(True)
    ax.legend(loc=1)

def calc_window(src_orbit, dst_orbit, t0):
    # angle difference at zero time
    d_ang_0 = dst_orbit.calc_mean_anomaly(t0) - src_orbit.calc_mean_anomaly(t0)
    assert(src_orbit.primary == dst_orbit.primary)
    primary = src_orbit.primary

    # First iteration
    a,e,n,t_h = calc_hohmann(src_orbit, dst_orbit, t0)
    # Angle difference at ideal launch time, calculated with transfer time
    d_ang = pi - t_h*dst_orbit.n
    # If the target is already past the position, the next opportunity must be searched
    if d_ang_0 > d_ang:
        d_ang = d_ang + 2*pi
    # time until next position
    # delta_ang0 + (Eve.n - Kerbin.n) * t_launchdow = delta_ang
    t_launch = t0 + (d_ang - d_ang_0) / abs(dst_orbit.n - src_orbit.n)
    t_arrival = t_launch + t_h
    ang_launch,r = src_orbit.calc_polar(t_launch)
    # If the dest is on a lower orbit, then the SV starts from the apoapsis of the transfer orbit, 
    # so a half-orbit offset in its omega parameter and true anomaly is needed.
    if src_orbit.a > dst_orbit.a:
        ang_offset = pi
    else:
        ang_offset = 0

    Hohmann = orbit(primary, a,e, omega = (ang_launch+ang_offset)*180/pi, t0=ang_offset-t_launch*n)
    Hohmann.t_launch = t_launch
    
    # initial value for the while loop
    t_miss = 1000
    # Real iteration
    while abs(t_miss) > 100:
        # calculate timing error - by how much time we've missed the target when arriving
        t_arrival = t_launch + t_h
        ang_h = np.mod(Hohmann.calc_polar(t_arrival)[0], 2*np.pi)
        ang_dst = np.mod(dst_orbit.calc_polar(t_arrival)[0], 2*np.pi)
        t_miss = (ang_h - ang_dst)/dst_orbit.n
        # modify start time using calculated error
        # sign depends on relation of angular velocities of the src and dest orbits
        if src_orbit.a > dst_orbit.a:
            t_launch = t_launch + t_miss
        else:
            t_launch = t_launch - t_miss
        
        # recalculate transfer orbit
        a,e,n,t_h = calc_hohmann(src_orbit, dst_orbit, t_launch)
        ang_launch,r = src_orbit.calc_polar(t_launch)
        Hohmann = orbit(primary, a,e, omega = (ang_launch+ang_offset)*180/pi, t0=ang_offset-t_launch*n)
        Hohmann.t_launch = t_launch
    
    # Prepare output
    Hohmann.recalc_orbit_visu(t_launch,t_arrival)
    print("{0:.1f} s, {1:.2f} d".format(t_launch, t_launch/3600/6))
    return Hohmann
# %% 
# Testing during development
if __name__ == "__main__":
    #Moho = orbit(a=5263138304, e=0.2, omega=70)
    Eve = planetary_body(orbit(Kerbol, a=9832684544, e=0.01, omega=15, i=2.1),
                        GM=8.1717302e12, radius=7e5, atmo_height=9e4)
    Kerbin = planetary_body(orbit(Kerbol, a=13599840256, e=0),
                            GM=3.5316e12, radius=6e5, atmo_height=7e4)
    Duna = planetary_body(orbit(Kerbol, a=20726155264, e=0.051, omega=135.5, i=0.06),
                          GM=3.0136321e11, radius=3.2e5, atmo_height=5e4)
    Mun = planetary_body(orbit(Kerbin, a=1.2e6, e=0, t0=1.7),
                         GM=6.5138398e10, radius=2e5, atmo_height=0)
    Minmus = planetary_body(orbit(Kerbin, a=4.7e7, e=0, t0=0.9),
                         GM=1.7658e9, radius=6e4, atmo_height=0)
    
    plot_orbits([Kerbin.orbit, Eve.orbit, calc_window(Kerbin.orbit, Eve.orbit, 0)])
    plot_orbits([Kerbin.orbit, Duna.orbit, calc_window(Kerbin.orbit, Duna.orbit, 0)])

# %%
relay_orbit = calc_orbit(Minmus,T=3600*9, e=0.0)
relay_orbit.min_alt
temporary_orbit = calc_orbit(Minmus,T=3600*12, rp=relay_orbit.rp)
temporary_orbit.min_alt
temporary_orbit.calc_circularization_pe()
# %%
transfer_orbit = calc_window(Kerbin.orbit, Duna.orbit, 0)
# %%
transfer_orbit.calc_distance_to(Kerbin, 235.55*6*3600)
# %%
