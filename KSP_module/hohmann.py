import numpy as np 
from math import sqrt, pi
from KSP_module.main import orbit, planetary_body
from KSP_module.system import Kerbin

# tolerance for Hohmann transfer calculation
MISS_TOL = 5  # seconds

def calc_hohmann(src_orbit, dst_orbit, t0):
    """
    Calculate semi-major axis, eccentricity, average angular velocity and time for Hohmann transfer.
    """
    assert(src_orbit.primary == dst_orbit.primary)
    primary = src_orbit.primary
    GM = primary.GM

    # position of src known at t0
    r_src = src_orbit.calc_polar(t0)[0]
    
    # initial guess for dest is using 0 time for Hohmann transfer
    hohmann_time = 0
    hohmann_time_prev = 10

    # Iterate a better Hohmann time by recalculating destination
    while abs(hohmann_time_prev-hohmann_time) > 1:
        r_dst = dst_orbit.calc_polar(t0+hohmann_time)[0]
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

def calc_window(src_orbit, dst_orbit, t0):
    """Transfer window calculation for Hohmann transfer"""
    assert(src_orbit.primary == dst_orbit.primary)
    primary = src_orbit.primary

    # angle difference at zero time
    d_ang_0 = dst_orbit.calc_mean_anomaly(t0) - src_orbit.calc_mean_anomaly(t0)

    # calculate period of search: time to get to the same angle difference as at zero time
    t_window_period = 2*pi / abs(dst_orbit.n - src_orbit.n)
    # First iteration
    a,e,n,t_h = calc_hohmann(src_orbit, dst_orbit, t0)
    # Approximation of angle difference at ideal launch time
    # calculated with first iteration of transfer time and angular velocities assuming circular orbits
    d_ang = pi - t_h*dst_orbit.n
    # If the target is already past the position, the next opportunity must be searched
    if d_ang_0 > d_ang:
        d_ang = d_ang + 2*pi
    # time until next position
    t_launch = t0 + (d_ang - d_ang_0) / abs(dst_orbit.n - src_orbit.n)
    t_arrival = t_launch + t_h
    ang_launch = src_orbit.calc_polar(t_launch)[1]
    # If the dest is on a lower orbit, then the SV starts from the apoapsis of the transfer orbit, 
    # so a half-orbit offset in its omega parameter and true anomaly is needed.
    if src_orbit.a > dst_orbit.a:
        ang_offset = pi
    else:
        ang_offset = 0

    kwargs = {'a':a, 'e':e, 'omega':(ang_launch+ang_offset)*180/pi, 't0':t_launch, 'nu0':ang_offset}
    Hohmann = orbit(primary, **kwargs)
    Hohmann.t_launch = t_launch
    
    # initial value for the while loop
    t_miss = 1000
    # Real iteration
    while abs(t_miss) > MISS_TOL:
        # calculate timing error - by how much time we've missed the target when arriving
        t_arrival = t_launch + t_h
        ang_h = np.mod(Hohmann.calc_polar(t_arrival)[1], 2*np.pi)
        ang_dst = np.mod(dst_orbit.calc_polar(t_arrival)[1], 2*np.pi)
        ang_miss = ang_h - ang_dst
        
        # modify start time using calculated error
        # use the angular velocity of the faster orbit
        t_miss = ang_miss/max(src_orbit.n, dst_orbit.n)
        if src_orbit.n > dst_orbit.n:
            t_miss = -t_miss

        # iterate launch time with error
        t_launch = t_launch + t_miss
        
        # if the launch time is before t0, wrap around by adding the window period
        if t_launch < t0:
            t_launch = t_launch + t_window_period
        
        # recalculate transfer orbit
        a,e,n,t_h = calc_hohmann(src_orbit, dst_orbit, t_launch)
        ang_launch = src_orbit.calc_polar(t_launch)[1]
 
        kwargs = {'a':a, 'e':e, 'omega':(ang_launch+ang_offset)*180/pi, 't0':t_launch, 'nu0':ang_offset}
        Hohmann = orbit(primary, **kwargs)
        Hohmann.t_launch = t_launch
    
    # Prepare output
    Hohmann.recalc_orbit_visu(t_launch,t_arrival)
    Hohmann.t_launch = t_launch
    Hohmann.t_arrival = t_arrival
    Hohmann.is_hohmann = True
    Hohmann.dv = Hohmann.calc_speed(t_launch) - src_orbit.calc_speed(t_launch)  # delta-v needed for transfer
    return Hohmann
