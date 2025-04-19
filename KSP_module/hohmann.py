# List of functions in this file:
# 1. calc_hohmann(src_orbit, dst_orbit, t0): Calculates the Hohmann transfer orbit.
# 2. calc_window(src_orbit, dst_orbit, t0): Calculates the transfer window for a Hohmann transfer.
# 3. calc_hohmann_dv(src_orbit, dst_orbit): Calculates the delta-v for Hohmann escape or capture.

import numpy as np 
from math import sqrt, pi
from KSP_module.main import orbit, planetary_body

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
    kwargs['omega'] = np.mod((ang_launch+ang_offset)*180/pi, 360)
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
            
        # store previous miss time
        t_miss_prev = t_miss
        # calculate new miss time 
        # using the angular velocity of the faster orbit
        t_miss = ang_miss/max(src_orbit.n, dst_orbit.n)
        # sign depends on which is the faster: source or destination orbit
        if src_orbit.n > dst_orbit.n:
            t_miss = -t_miss

        # if the sign of the previous and current miss time is different,
        # then the launch time is between the two
        # under relaxation is used in this case to avoid unstable oscillations
        if t_miss_prev * t_miss < 0:
            t_miss = 0.8*t_miss

        # use miss time to recalculate launch time
        t_launch = t_launch + t_miss
        # if the launch time is before t0, wrap around
        t_launch = t_launch + t_window_period if t_launch < t0 else t_launch
    
    # Prepare output
    t_arrival = Hohmann.t_arrival
    Hohmann.recalc_orbit_visu(Hohmann.t_launch,Hohmann.t_arrival)
    # delta-v needed for leaving source body orbit
    Hohmann.leave_dv = abs(Hohmann.calc_speed(t_launch) - src_orbit.calc_speed(t_launch))
    # delta-v needed to match target body orbit
    Hohmann.enter_dv = abs(Hohmann.calc_speed(t_arrival) - dst_orbit.calc_speed(t_arrival))
    return Hohmann

def calc_hohmann_dv(src_orbit, dst_orbit):
    """Calculate delta-v for Hohmann escape or capture, from/to an orbit around the source or destination body"""
    
    # either source or destination must be a Hohmann orbit (exactly one)
    assert(src_orbit.is_hohmann ^ dst_orbit.is_hohmann), 'Exactly one orbit must be a Hohmann orbit'
    # chack which one is which
    if src_orbit.is_hohmann:
        hohmann_orbit = src_orbit
        parking_orbit = dst_orbit
        # t_calc is to check that the destination and the Hohmann orbit are at the same position
        # at the time of arrival, and calculate the dv
        t_calc = hohmann_orbit.t_arrival
    else:
        hohmann_orbit = dst_orbit
        parking_orbit = src_orbit
        t_calc = hohmann_orbit.t_launch

    # parking is around a body orbiting the same primary as the Hohmann orbit
    parking_body = parking_orbit.primary
    assert(parking_body.orbit.primary == hohmann_orbit.primary),\
        'Parking orbit must be around a body orbiting the same primary as the Hohmann orbit'
    
    #TODO check that the two orbits are plausible
    # Problem: inclination is not yet handled
    parking_p0 = parking_body.orbit.calc_xyz(t_calc)
    hohmann_p0 = hohmann_orbit.calc_xyz(t_calc)
    
    # calculate dv to/from Hohmann orbit, relative to the primary body of the parking orbit
    # this is the speed when leaving/entering the SOI
    v_soi = abs(hohmann_orbit.calc_speed(t_calc) - parking_body.orbit.calc_speed(t_calc))
    r_soi = parking_body.soi

    v_parking = parking_orbit.calc_speed(t_calc)
    r_parking = parking_orbit.a
    
    # calculate speed at parking orbit height using vis-viva equation
    GM = parking_body.GM
    v_hohmann_at_parking = sqrt(v_soi**2 + 2*GM/r_parking - 2*GM/r_soi)
    return v_hohmann_at_parking - v_parking