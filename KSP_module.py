import numpy as np 

class star:
    def __init__(self, GM=1.1723328e18, radius=2.616e8):
        self.GM = GM
        self.radius = radius

Kerbol = star()

class body:
    """Planets, moons, spacecrafts, etc."""
    def __init__(self, orbit):
        self.primary = orbit.primary
        self.orbit = orbit


class planetary_body(body):
    def __init__(self, orbit, GM, radius=0, atmo_height=0):
        super().__init__(orbit)
        self.GM = GM
        self.radius = radius
        self.atmo_height = atmo_height

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

    def calc_speed_at(self, time):
        r = self.calc_polar(time)[1]
        return np.sqrt(self.primary.GM*(2/r - 1/self.a))
    
    def calc_circularization_ap(self):
        """Calculate the delta-v needed for circularization burn at apoapsis""" 
        v_target = np.sqrt(self.primary.GM/self.ra)
        return v_target - self.va

    def calc_circularization_pe(self):
        """Calculate the delta-v needed for circularization burn at periapsis""" 
        v_target = np.sqrt(self.primary.GM/self.rp)
        return self.vp - v_target
    
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
