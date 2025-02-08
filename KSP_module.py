import numpy as np 
import matplotlib.pyplot as plt
from math import sqrt, pi, sin, cos, atan

class star:
    def __init__(self, name, GM=1.1723328e18, radius=2.616e8):
        self.__name__ = name
        self.GM = GM
        self.radius = radius
        self.secondaries = []
    
    def __str__(self):
        return self.__name__
    
    def __repr__(self):
        return self.__name__

    def calc_xy(self, time):
        """The star is at the origin"""
        return (0,0)

Kerbol = star('Kerbol')

class body:
    """Planets, moons, spacecrafts, etc."""
    def __init__(self, orbit):
        self.primary = orbit.primary
        self.orbit = orbit
        # self.name = 


class planetary_body(body):
    def __init__(self, name, orbit, GM, radius=0, atmo_height=0):
        self.__name__ = name
        super().__init__(orbit)
        self.GM = GM
        self.radius = radius
        self.year = orbit.T
        self.atmo_height = atmo_height
        self.soi = self.orbit.a*(self.GM/self.primary.GM)**(2/5)
        self.secondaries = []
        orbit.primary.secondaries.append(self)

    def __str__(self):
        return self.__name__
    
    def __repr__(self):
        return self.__name__

    def calc_orbital_velocity(self, altitude):
        if altitude > self.atmo_height:
            return np.sqrt(self.GM/(self.radius + altitude))
        else:
            raise ValueError("Altitude is inside atmosphere")
    
#todo: add RAAN
class orbit:
    """Default orbit is the circular orbit of Kerbin around Kerbol
    Inputs are a bit chaotic but conform to KSP's way of defining orbits:
    primary: body around which the orbit is
    a: semi-major axis in meters
    e: eccentricity
    i: inclination in degrees
    omega: argument of periapsis in degrees
    t0: time of orbit initialization
    nu0: true anomaly at initialization in radians"""
    def __init__(self, 
                 primary=Kerbol, 
                 a=13599840256, 
                 e=0,
                 i=0,
                 omega=0,
                 t0=0,
                 nu0=3.14):
        self.a = a # semi-major axis
        self.e = e # eccentricity
        self.i = i*np.pi/180 # inclination in radians
        self.omega = omega*np.pi/180 # argument of periapsis in radians
        self.t0 = t0 # time of orbit initialization
        self.primary = primary
        self.T = np.sqrt(4*np.pi**2*a**3/self.primary.GM) # orbital period, seconds
        self.n = 2*np.pi/self.T # angular velocity, rad/s
        if nu0 == 0:
            self.t_epoch = t0
        else:
            self.t_epoch = self.calc_epoch_time(nu0)
        self.nu0 = nu0 # true anomaly at epoch

        self.rp = a*(1-e) # periapsis
        self.ra = a*(1+e) # apoapsis
        self.vp = np.sqrt(primary.GM*(1+self.e)/(self.a*(1-self.e))) # velocity at periapsis
        self.va = np.sqrt(primary.GM*(1-self.e)/(self.a*(1+self.e))) # velocity at apoapsis
        self.min_alt = self.rp - primary.radius
        self.max_alt = self.ra - primary.radius
        
        self.is_elliptic = self.check_elliptic()
        self.is_hohmann = False # default setting

        # calculate orbit for visualization
        self.t = np.linspace(0, self.T, 100)
        self.phi, self.r = self.calc_polar(self.t)

    # If the mechanical energy is negative, the orbit is elliptic
    def check_elliptic(self):
        return self.vp**2/2-self.primary.GM/self.rp < 0

    def calc_mean_anomaly(self, time):
        return (time-self.t0)*self.n + self.nu0

#todo: scipy.optimize.newton
#todo: test function w/ calc_window
    def calc_eccentric_anomaly(self, time):
        """M: mean anomaly, e: eccentricity, E: eccentric anomaly"""
        def f(E,M,e): return E - e*np.sin(E) - M
        def dfdE(E,M,e): return 1 - e*np.cos(E)
        M = self.calc_mean_anomaly(time)
        E = M  # initial guess
        max_delta = 1 # initial error
        while max_delta > 1e-9:
            E_old = E
            E = E - f(E,M,self.e)/dfdE(E,M,self.e)
            max_delta = np.max(np.abs(E_old - E))
        return E

    def calc_true_anomaly(self, time):
        E = self.calc_eccentric_anomaly(time)
        beta = self.e/(1+np.sqrt(1-self.e*self.e))
        return E + 2*np.arctan( beta*np.sin(E) / (1-beta*np.cos(E)) )

    def calc_polar(self, time):
        """function to calculate both orbital distance and angle wrt reference angle"""
        nu = self.calc_true_anomaly(time)
        phi = self.omega + nu
        r = self.a*(1-self.e*self.e)/(1+self.e*np.cos(nu))
        return phi,r
    
    def calc_epoch_time(self, nu0):
        """Calculate epoch time from true anomaly"""
        E0 = np.arccos((self.e+np.cos(nu0)) / (1+self.e*np.cos(nu0)))
        M0 = E0 - self.e*np.sin(E0)
        return self.t0-M0/self.n

    def recalc_orbit_visu(self, start_time, end_time):
        """function to recalculate orbit for visualization"""
        self.t = np.linspace(start_time, end_time, 50)
        self.phi, self.r = self.calc_polar(self.t)

    def calc_speed(self, time):
        r = self.calc_polar(time)[1]
        return np.sqrt(self.primary.GM*(2/r - 1/self.a))
    
    def calc_zenith_angle(self, time):
        """Calculate the angle between the velocity vector and the radius vector"""
        r1 = self.calc_polar(time)[1]
        v1 = self.calc_speed(time)
        return np.arcsin((self.vp*self.rp) / (v1*r1))
        # nu = self.calc_true_anomaly(time)
        # return np.arccos(self.e + np.cos(nu))/(1 + self.e*np.cos(nu))
    
    def calc_circularization_ap(self):
        """Calculate the delta-v needed for circularization burn at apoapsis""" 
        v_target = np.sqrt(self.primary.GM/self.ra)
        return v_target - self.va

    def calc_circularization_pe(self):
        """Calculate the delta-v needed for circularization burn at periapsis""" 
        v_target = np.sqrt(self.primary.GM/self.rp)
        return self.vp - v_target
    
    def calc_xy(self, time):
        """Calculate Cartesian coordinates for a given time"""
        phi,r = self.calc_polar(time)
        return np.asarray(r*np.cos(phi), r*np.sin(phi))

    def calc_distance_to(self, other, time):
        """Calculate distance between two objects at a given time"""
        if isinstance(other, orbit):
            other_orbit = other
        elif isinstance(other, body) or isinstance(other, planetary_body):
            other_orbit = other.orbit
        else:
            raise(TypeError)
        return np.linalg.norm(self.calc_xy(time)-other_orbit.calc_xy(time))
    
    def do_maneuver(self, time, longitudinal_dv, lateral_dv=0, radial_dv=0):
        # radius at maneuver
        phi, r = self.calc_polar(time)
        # flight path angle at maneuver
        fphi = pi/2 - self.calc_zenith_angle(time)
        # original speed at maneuver
        v = self.calc_speed(time)
        # add the three components
        v = v + longitudinal_dv
        #TODO: lateral and radial to be implemented
        # recalculate orbital parameters
        e = sqrt((r*v**2/self.primary.GM - 1)**2 * cos(fphi)**2 + sin(fphi)**2)
        a = 1/(2/r - v**2/self.primary.GM)
        mu = atan( r*v**2/self.primary.GM*cos(fphi)*sin(fphi) / 
                  (r*v**2/self.primary.GM*cos(fphi)**2 - 1) )
        return a,e,mu

def calc_orbit(primary,a=False,e=False,T=False,rp=False,ra=False):
    """Calculate orbital parameters from given data"""

    if bool(T):
        a = (T/(2*np.pi))**(2/3)*primary.GM**(1/3)
    if bool(rp) and bool(ra):
        a = (rp + ra)/2
    if bool(rp) and bool(e):
        a = rp/(1-e)
    if bool(ra) and bool(e):
        a = ra/(1+e)
    
    if not(isinstance(e, bool)):
        rp = a*(1-e)
        ra = a*(1+e)

    if not(isinstance(ra, bool)):
        e = ra/a-1
        rp = a*(1-e)

    if not(isinstance(rp, bool)):
        e = 1-rp/a
        ra = a*(1+e)

    return orbit(primary, a, e)

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

def plot_hohmann_orbit(src_body, dst_body, transfer_orbit):
    """Plot Hohmann transfer, 0: source body, 1: destination body, 2: transfer orbit """
    fig,ax = initialize_plot()
    ax = add_planetary_body_to_plot(ax, src_body)
    ax = add_planetary_body_to_plot(ax, dst_body)
    ax = add_orbit_to_plot(ax, transfer_orbit)
    ax = add_point_to_plot(ax, transfer_orbit.calc_polar(transfer_orbit.t_launch))
    ax = add_point_to_plot(ax, transfer_orbit.calc_polar(transfer_orbit.t_arrival))
    ax = add_point_to_plot(ax, dst_body.orbit.calc_polar(transfer_orbit.t_launch), f"{dst_body} at launch")
    return fig,ax

def calc_window(src_orbit, dst_orbit, t0):
    """Transfer window calculation for Hohmann transfer"""
    assert(src_orbit.primary == dst_orbit.primary)
    primary = src_orbit.primary

    # angle difference at zero time
    d_ang_0 = dst_orbit.calc_mean_anomaly(t0) - src_orbit.calc_mean_anomaly(t0)

    # First iteration
    a,e,n,t_h = calc_hohmann(src_orbit, dst_orbit, t0)
    # Approximation of angle difference at ideal launch time
    # calculated with first iteration of transfer time and angular velocities assuming circular orbits
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

    Hohmann = orbit(primary, a,e, omega = (ang_launch+ang_offset)*180/pi, t0=t_launch, nu0=ang_offset)
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
        Hohmann = orbit(primary, a,e, omega = (ang_launch+ang_offset)*180/pi, t0=t_launch, nu0=ang_offset)
        Hohmann.t_launch = t_launch
    
    # Prepare output
    Hohmann.recalc_orbit_visu(t_launch,t_arrival)
    Hohmann.t_launch = t_launch
    Hohmann.t_arrival = t_arrival
    Hohmann.is_hohmann = True
    return Hohmann

def time(date):
    """Convert date to seconds for KSP. Input is a list of integers: [year, day, hour, minute, second]"""
    return (date[0]*Kerbin.year + date[1]*60*60*6 + date[2]*60*60 + date[3]*60 + date[4])

def date(time):
    """Convert seconds to date for KSP. Output is a list of integers: [year, day, hour, minute, second]"""
    #TODO: check if the year is correct
    year = time // Kerbin.year
    time = time % Kerbin.year
    day = time//(60*60*6)
    time = time%(60*60*6)
    hour = time//(60*60)
    time = time%(60*60)
    minute = time//60
    second = time%60
    return [year, day, hour, minute, second]

def pretty_date(time):
    [year, day, hour, minute, second] = date(time)
    print(f'{year:.0f} y {day:.0f} d {hour:.0f} h {minute:.0f} m {second:.0f} s')

def initialize_plot():
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.grid(True, alpha=0.2)
    ax.set_yticklabels([])  # Hide radial labels
    return fig, ax

def add_planetary_body_to_plot(ax, body):
    ax.plot(body.orbit.phi, body.orbit.r, label=body.__name__)
    ax.legend(loc=1)
    return ax

def add_orbit_to_plot(ax, orbit, label=None):
    ax.plot(orbit.phi, orbit.r, label=label)
    ax.legend(loc=1)
    return ax

def add_point_to_plot(ax, coordinates, label=None, marker='o'):
    (phi, r) = coordinates
    ax.plot(phi, r, marker, label=label)
    ax.legend(loc=1)
    return ax

# Planetary bodies
Eve = planetary_body('Eve', orbit(Kerbol, a=9832684544, e=0.01, omega=15, i=2.1),
                    GM=8.1717302e12, radius=7e5, atmo_height=9e4)
Kerbin = planetary_body('Kerbin', orbit(Kerbol, a=13599840256, e=0),
                        GM=3.5316e12, radius=6e5, atmo_height=7e4)
Duna = planetary_body('Duna', orbit(Kerbol, a=20726155264, e=0.051, omega=135.5, i=0.06),
                        GM=3.0136321e11, radius=3.2e5, atmo_height=5e4)
Mun = planetary_body('Mun', orbit(Kerbin, a=1.2e6, e=0, nu0=1.7),
                        GM=6.5138398e10, radius=2e5, atmo_height=0)
Minmus = planetary_body('Minmus', orbit(Kerbin, a=4.7e7, e=0, nu0=0.9),
                        GM=1.7658e9, radius=6e4, atmo_height=0)