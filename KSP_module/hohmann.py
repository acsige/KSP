import numpy as np 
from math import sqrt, pi
from KSP_module.main import orbit, planetary_body
from KSP_module.system import Kerbin

# tolerance for Hohmann transfer calculation
MISS_TOL = 5  # seconds

def calc_hohmann(src_orbit, dst_orbit, t0):
    """
    Calculate orbit for Hohmann transfer at t0. Destination position is found by iterating.
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
        T = sqrt(4*pi**2*a**3/GM)
        hohmann_time_prev = hohmann_time
        hohmann_time = T / 2    
        
    # initialize orbit
    # If the dest is on a lower orbit, then the SV starts from the apoapsis of the transfer orbit, 
    # so a half-orbit offset in its omega parameter and true anomaly is needed.

    if r_src > r_dst:
        kwargs = {'ra':r_src, 'rp':r_dst}
        ang_offset = pi
    else:   
        kwargs = {'ra':r_dst, 'rp':r_src}
        ang_offset = 0
    
    ang_launch = src_orbit.calc_polar(t0)[1]
    kwargs['omega'] = (ang_launch+ang_offset)*180/pi
    kwargs['t0'] = t0
    kwargs['nu0'] = ang_offset
    result_orbit = orbit(primary, **kwargs)
    
    result_orbit.t_launch = t0
    result_orbit.t_arrival = t0 + result_orbit.T/2
    result_orbit.is_hohmann = True

    return result_orbit

def calc_window(src_orbit, dst_orbit, t0):
    """Transfer window calculation for Hohmann transfer"""
    assert(src_orbit.primary == dst_orbit.primary)

    # calculate period of search: time to get to the same angle difference as at zero time
    t_window_period = 2*pi / abs(dst_orbit.n - src_orbit.n)    
    # First calculation - transfer orbit at start of search window
    Hohmann = calc_hohmann(src_orbit, dst_orbit, t0)
    # calculate time error - by how much time we've missed the target when arriving its orbit
    ang_h = np.mod(Hohmann.calc_polar(Hohmann.t_arrival)[1], 2*np.pi)
    ang_dst = np.mod(dst_orbit.calc_polar(Hohmann.t_arrival)[1], 2*np.pi)
    ang_miss = ang_h - ang_dst
    
    # modify start time using calculated error
    # use the angular velocity of the faster orbit
    t_miss = ang_miss/max(src_orbit.n, dst_orbit.n)
    if src_orbit.n > dst_orbit.n:
        t_miss = -t_miss

    # if the new launch time is before t0, wrap around by adding the window period
    if t_miss < 0:
        t_miss = t_miss + t_window_period

    # first launch time
    t_launch = t0 + t_miss
    # Real iteration
    while abs(t_miss) > MISS_TOL:

        # possible solution: wrap around by adding the window period
        #if t_launch < t0:
        #    t_launch = t_launch + t_window_period
        
        # recalculate Hohmann orbit
        Hohmann = calc_hohmann(src_orbit, dst_orbit, t_launch)
        # calculate timing error
        ang_h = np.mod(Hohmann.calc_polar(Hohmann.t_arrival)[1], 2*np.pi)
        ang_dst = np.mod(dst_orbit.calc_polar(Hohmann.t_arrival)[1], 2*np.pi)
        ang_miss = ang_h - ang_dst
        # convert miss angle between -pi and pi
        if ang_miss > pi:
            ang_miss = ang_miss - 2*pi
        elif ang_miss < -pi:
            ang_miss = ang_miss + 2*pi
            
        # modify start time using calculated error
        # use the angular velocity of the faster orbit
        t_miss = ang_miss/max(src_orbit.n, dst_orbit.n)
        if src_orbit.n > dst_orbit.n:
            t_miss = -t_miss
        # use miss time to recalculate launch time
        t_launch = t_launch + t_miss
        # if the launch time is before t0, wrap around
        t_launch = t_launch + t_window_period if t_launch < t0 else t_launch
        #  stop the iteration
        assert (t_launch >= t0), "Iterated launch time is before start time"
    
    # Prepare output
    Hohmann.recalc_orbit_visu(Hohmann.t_launch,Hohmann.t_arrival)
    Hohmann.dv = Hohmann.calc_speed(t_launch) - src_orbit.calc_speed(t_launch)  # delta-v needed for transfer
    return Hohmann
