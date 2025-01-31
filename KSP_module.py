import numpy as np 
import matplotlib.pyplot as plt
from math import sqrt, pi

class star:
    def __init__(self, GM=1.1723328e18, radius=2.616e8):
        self.GM = GM
        self.radius = radius
        self.secondaries = []

    def calc_xy(self, time):
        """The star is at the origin"""
        return (0,0)

Kerbol = star()

class body:
    """Planets, moons, spacecrafts, etc."""
    def __init__(self, orbit):
        self.primary = orbit.primary
        self.orbit = orbit
        # self.name = 


class planetary_body(body):
    def __init__(self, orbit, GM, radius=0, atmo_height=0):
        super().__init__(orbit)
        self.GM = GM
        self.radius = radius
        self.atmo_height = atmo_height
        self.soi = self.orbit.a*(self.GM/self.primary.GM)**(2/5)
        self.secondaries = []
        orbit.primary.secondaries.append(self)

    def calc_orbital_velocity(self, altitude):
        if altitude > self.atmo_height:
            return np.sqrt(self.GM/(self.radius + self.atmo_height))
        else:
            raise ValueError("Altitude is inside atmosphere")
    
class orbit:
    """Default orbit is the circular orbit of Kerbin around Kerbol"""
    def __init__(self, 
                 primary=Kerbol, 
                 a=13599840256, 
                 e=0,
                 i=0,
                 omega=0,
                 t0=3.14):
        self.a = a # semi-major axis
        self.e = e # eccentricity
        self.i = i*np.pi/180 # inclination in radians
        self.omega = omega*np.pi/180 # argument of periapsis in radians
        self.t0 = t0 # epoch
        self.primary = primary
        self.rp = a*(1-e) # periapsis
        self.ra = a*(1+e) # apoapsis
        self.vp = np.sqrt(primary.GM*(1+self.e)/(self.a*(1-self.e))) # velocity at periapsis
        self.va = np.sqrt(primary.GM*(1-self.e)/(self.a*(1+self.e))) # velocity at apoapsis

        self.min_alt = self.rp - primary.radius
        self.max_alt = self.ra - primary.radius

        self.T = np.sqrt(4*np.pi**2*a**3/self.primary.GM) # orbital period, seconds
        self.n = 2*np.pi/self.T # angular velocity, rad/s
        # calculate orbit for visualization
        self.t = np.linspace(0, self.T, 100)
        self.phi, self.r = self.calc_polar(self.t)

    def calc_mean_anomaly(self, time):
        return time*self.n + self.t0

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

    def recalc_orbit_visu(self, start_time, end_time):
        """function to recalculate orbit for visualization"""
        self.t = np.linspace(start_time, end_time, 50)
        self.phi, self.r = self.calc_polar(self.t)

    def calc_speed(self, time):
        r = self.calc_polar(time)[1]
        return np.sqrt(self.primary.GM*(2/r - 1/self.a))
    
    def calc_gamma(self, time):
        """Calculate the angle between the velocity vector and the radius vector"""
        r1 = self.calc_polar(time)[1]
        v1 = self.calc_speed(time)
        return np.arcsin((self.vp*self.rp)/v1*r1)
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