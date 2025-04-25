import numpy as np 
from math import sqrt, pi, sin, cos, tan, atan, atan2, acos

class star:
    def __init__(self, name, GM=1.1723328e18, radius=2.616e8):
        self.__name__ = name
        self.GM = GM
        self.radius = radius
        self.secondaries = []
        self.atmo_height = 0
    
    def __str__(self):
        return self.__name__
    
    def __repr__(self):
        return self.__name__

    def calc_xyz(self, time):
        """The star is at the origin"""
        return (0,0,0)

class body:
    """Planets, moons, spacecrafts, etc."""
    def __init__(self, orbit):
        self.primary = orbit.primary
        self.orbit = orbit

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

class orbit:
    """Class to hold orbital parameters and calculate orbital positions.
    Inputs are a bit chaotic but conform to KSP's way of defining orbits:
    primary: body around which the orbit is
    a: semi-major axis in meters
    e: eccentricity
    i: inclination in degrees
    omega: argument of periapsis in degrees
    OMEGA: RAAN in degrees
    t0: time of orbit initialization
    nu0: true anomaly at initialization in radians
    
    List of functions:
    1. calc_missing_parameters: Calculates missing parameters from given data.
    2. calc_orbit_from_burnout_parameters: Calculates orbital parameters from burnout parameters.
    3. check_elliptic: Checks if the orbit is elliptic.
    4. check_hyperbolic: Checks if the orbit is hyperbolic.
    5. calc_mean_anomaly: Calculates the mean anomaly.
    6. calc_eccentric_anomaly: Calculates the eccentric anomaly.
    7. calc_true_anomaly: Calculates the true anomaly.
    8. calc_polar: Calculates the polar coordinates.
    9. calc_xyz: Calculates the Cartesian coordinates.
    10. calc_epoch_time: Calculates the epoch time from true anomaly.
    11. recalc_orbit_visu: Recalculates the orbit for visualization.
    12. calc_speed: Calculates the speed at a given time.
    13. calc_zenith_angle: Calculates the zenith angle.
    14. calc_circularization_ap: Calculates the delta-v needed for circularization burn at apoapsis.
    15. calc_circularization_pe: Calculates the delta-v needed for circularization burn at periapsis.
    16. calc_distance_to: Calculates the distance to another object at a given time.
    17. calc_min_distance_to: Calculates the minimum distance to another object during a given time interval.
    18. calc_soi_leave_time: Calculates the time when the spacecraft leaves the SOI of the primary body.
    19. calc_soi_change: Calculates the time of SOI change and the new primary.
    20. calc_soi_enter: Calculates the time of entering a secondary's SOI.
    21. do_maneuver: Performs a maneuver and returns the new orbital parameters.
    22. calc_soi_change: Calculates the time of SOI change.
    """
    def __init__(self, primary, **kwargs):
        self.primary = primary
        # either initialize orbit from state vector...
        if 'r' in kwargs and 'v' in kwargs:
            orbit_kwargs_burnout = self.calc_orbit_from_burnout_parameters(**kwargs)
            self.orbit_kwargs = self.calc_missing_parameters(**orbit_kwargs_burnout)
        # or initialize from orbital parameters
        else:
            self.orbit_kwargs = self.calc_missing_parameters(**kwargs)

        self.a = self.orbit_kwargs['a'] # semi-major axis
        self.e = self.orbit_kwargs['e'] # eccentricity
        self.i = self.orbit_kwargs['i']*pi/180 # inclination in radians
        self.omega = self.orbit_kwargs['omega']*pi/180 # argument of periapsis in radians
        self.OMEGA = self.orbit_kwargs['OMEGA']*pi/180 # RAAN in radians
        self.OMEGA_projected = atan2(sin(self.OMEGA)*cos(self.i), cos(self.OMEGA))
        self.t0 = self.orbit_kwargs['t0'] # time of orbit initialization
        
        self.T = np.sqrt(4*np.pi**2*self.a**3/self.primary.GM) # orbital period, seconds
        self.n = 2*np.pi/self.T # angular velocity, rad/s
        if self.orbit_kwargs['nu0'] == 0:
            self.t_epoch = self.t0  # time of periapsis passage
        else:
            self.t_epoch = self.calc_epoch_time(self.orbit_kwargs['nu0'])
        self.nu0 = self.orbit_kwargs['nu0'] # true anomaly at epoch

        self.rp = self.a*(1-self.e) # periapsis
        self.ra = self.a*(1+self.e) # apoapsis
        self.vp = np.sqrt(primary.GM*(1+self.e)/(self.a*(1-self.e))) # velocity at periapsis
        self.va = np.sqrt(primary.GM*(1-self.e)/(self.a*(1+self.e))) # velocity at apoapsis
        
        self.min_alt = self.rp - primary.radius
        if self.min_alt < self.primary.atmo_height:
            raise ValueError("Minimum altitude is inside atmosphere")
        self.max_alt = self.ra - primary.radius
        
        self.is_elliptic = self.check_elliptic()
        self.is_hohmann = False # default setting

        self.recalc_orbit_visu(self.t0, self.t0+self.T)

    #TODO: this is just a stub for now, needs to be implemented
    def calc_orbit_from_burnout_parameters(self, v, r, t0=0, gamma=pi/2, beta=pi/2, delta=0, lambda2=0):
        """Calculate orbital parameters from burnout parameters. See http://www.braeunig.us/space/orbmech.htm
        Args:
            primary: planetary body
            v: velocity (magnitude only)
            r: radius (magnitude only)
            t0: time of orbit initialization
            gamma: zenith angle
            beta: azimuth heading measured in radians clockwise from north
            delta: geocentric latitude (or declination) of the burnout point in radians
            lambda2: longitude of the burnout point in radians (for simplification it's the celestial longitude)
        Returns:
            orbit_kwargs: dictionary with orbital parameters
        """
        #flight path angle
        phi = pi/2 - gamma
        # GM of the primary body
        GM = self.primary.GM

        # semi-major axis
        a = 1/(2/r - v**2/GM)
        # eccentricity
        e = sqrt((r*v**2/GM - 1)**2 * cos(phi)**2 + sin(phi)**2)
        # initial true anomaly
        nu0 = atan( r*v**2/GM*cos(phi)*sin(phi) /
                  (r*v**2/GM*cos(phi)**2 - 1) )
        # inclination
        i = acos(cos(delta)*sin(beta))
        # angular distance between the ascending node and the burnout point 
        # measured in the orbital plane
        l = atan2(tan(delta), cos(beta))
        # angular distance between the ascending node and the burnout point 
        # measured in the equatorial plane
        dl = atan(sin(delta)*tan(beta))
        # argument of periapsis in degrees
        omega = (l - nu0)*180/pi
        # longitude of the ascending node in radians
        lambda1 = lambda2 - dl
        # since longitudes are in celestial coordinates, lambda1 is RAAN
        # just convert here to degrees
        OMEGA = lambda1*180/pi
        
        return {
            'a': a,
            'e': e,
            't0': t0,
            'nu0': nu0,
            'i': i,
            'omega': omega,
            'OMEGA': OMEGA
        }

    def calc_missing_parameters(self, **kwargs):
        """Calculate parameters needed for orbit definition from given data.
        
        Parameters
        ----------
        primary : planetary_body
            The primary body around which the orbit is calculated.
        **kwargs : dict
            Keyword arguments for the orbit class.
        
        Returns
        -------
        orbit
            The calculated orbit.
        """
        
        # convert altitudes to apsis distances
        if 'min_alt' in kwargs:
            assert('rp' not in kwargs)
            kwargs['rp'] = self.primary.radius + kwargs['min_alt']

        if 'max_alt' in kwargs:
            assert('ra' not in kwargs)
            kwargs['ra'] = self.primary.radius + kwargs['max_alt']

        # calculate semi-major axis from different parameters
        if 'T' in kwargs:
            assert('a' not in kwargs)
            kwargs['a'] = (kwargs['T']/(2*pi))**(2/3)*self.primary.GM**(1/3)

        if 'rp' in kwargs and 'ra' in kwargs:
            assert('a' not in kwargs)
            assert('e' not in kwargs)
            kwargs['a'] = (kwargs['rp'] + kwargs['ra']) / 2
            kwargs['e'] = (kwargs['ra'] - kwargs['rp']) / (kwargs['ra'] + kwargs['rp'])

        if 'rp' in kwargs and 'e' in kwargs and 'a' not in kwargs:
            assert('ra' not in kwargs)
            kwargs['a'] = kwargs['rp'] / (1 - kwargs['e'])

        if 'ra' in kwargs and 'e' in kwargs and 'a' not in kwargs:
            assert('rp' not in kwargs)
            kwargs['a'] = kwargs['ra'] / (1 + kwargs['e'])
            
        # at this point we should have a and e
        # set defaults for missing values
        if 'i' not in kwargs:
            kwargs['i'] = 0

        if 'omega' not in kwargs:
            kwargs['omega'] = 0
        
        if 'OMEGA' not in kwargs:
            kwargs['OMEGA'] = 0

        if 't0' not in kwargs:
            kwargs['t0'] = 0

        if 'nu0' not in kwargs:
            kwargs['nu0'] = 0

        # define orbit with only the keys expected by the orbit class
        # no check for missing keys, the orbit class will raise an error
        orbit_input_keys = ['a', 'e', 'i', 'omega', 'OMEGA', 't0', 'nu0']
        orbit_kwargs = {k:kwargs[k] for k in orbit_input_keys}
        return orbit_kwargs

    # If the mechanical energy is negative, the orbit is elliptic
    def check_elliptic(self):
        return self.vp**2/2-self.primary.GM/self.rp < 0

    def check_hyperbolic(self):
        return not(self.check_elliptic())

    def calc_mean_anomaly(self, time):
        return (time-self.t_epoch)*self.n

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
        phi = self.OMEGA_projected + self.omega + nu
        r = self.a*(1-self.e*self.e)/(1+self.e*np.cos(nu))
        return r, phi
    
    def calc_xyz(self, time):
        """Calculate Cartesian coordinates for a given time"""
        r, phi = self.calc_polar(time)
        x = r*np.cos(phi)
        y = r*np.sin(phi)*cos(self.i)
        z = r*np.sin(phi)*sin(self.i)
        return np.asarray([x,y,z])

    def calc_epoch_time(self, nu0):
        """Calculate epoch time from true anomaly"""
        E0 = 2*np.arctan(np.tan(nu0/2)/np.sqrt((1+self.e)/(1-self.e)))
        M0 = E0 - self.e*np.sin(E0)
        return self.t0-M0/self.n

    def recalc_orbit_visu(self, start_time, end_time):
        """function to recalculate orbit for visualization"""
        self.t = np.linspace(start_time, end_time, 100)
        self.r, self.phi = self.calc_polar(self.t)

    def calc_speed(self, time):
        r = self.calc_polar(time)[0]
        return np.sqrt(self.primary.GM*(2/r - 1/self.a))
    
    def calc_zenith_angle(self, time):
        """Calculate the angle between the velocity vector and the radius vector"""
        r1 = self.calc_polar(time)[0]
        v1 = self.calc_speed(time)
        return np.arcsin((self.vp*self.rp) / (v1*r1))
    
    def calc_circularization_ap(self):
        """Calculate the delta-v needed for circularization burn at apoapsis""" 
        v_target = np.sqrt(self.primary.GM/self.ra)
        return v_target - self.va

    def calc_circularization_pe(self):
        """Calculate the delta-v needed for circularization burn at periapsis""" 
        v_target = np.sqrt(self.primary.GM/self.rp)
        return self.vp - v_target

    def calc_distance_to(self, other, time):
        """Calculate distance between two objects at a given time"""
        if isinstance(other, orbit):
            other_orbit = other
        elif isinstance(other, body) or isinstance(other, planetary_body):
            other_orbit = other.orbit
        else:
            raise(TypeError)
        
        # if both orbits are around the same primary
        if self.primary == other_orbit.primary:
            return np.linalg.norm(self.calc_xyz(time)-other_orbit.calc_xyz(time))
        # if one of the orbits is around the other
        elif self.primary == other:
            r,phi = self.calc_polar(time)
            return r
        else:
            raise ValueError("Usecase not implemented")
    
    def calc_min_distance_to(self, other, t_start, t_end):
        """Calculate minimum distance between two objects during a given time interval
        Args:
            other: other object (body or orbit)
            t_start: start time
            t_end: end time
        Returns:
            t: time of minimum distance
            d: minimum distance"""
        TOL = 1 # tolerance for time error

        t1, t2, delta_t = t_start, t_end, t_end - t_start
        d1 = self.calc_distance_to(other, t1)
        d2 = self.calc_distance_to(other, t2)

        while delta_t > TOL:
            t = (t1 + t2) / 2
            d = self.calc_distance_to(other, t)
            if d1 < d2:
                t2 = t
                d2 = d
            else:
                t1 = t
                d1 = d
            delta_t = t2 - t1
        return t, d
    
    def calc_soi_leave_time(self, t_start):
        """Calculate the time when the spacecraft leaves the SOI of the primary body
        Args:
            t_start: start time of calculation
        Returns:
            t_leave: time when spacecraft leaves the SOI"""

        # find initial t_end outside of SOI with timesteps doubling until outside
        dt = 1
        t_end = t_start + dt
        r_end = self.calc_polar(t_end)[0]
        r_soi = self.primary.soi
        while r_end < r_soi:
            dt *= 2
            t_end += dt
            r_end = self.calc_polar(t_end)[0]

        # find time when spacecraft leaves SOI with bisect method
        r = r_end
        while abs(r - r_soi) > 1:
            t = (t_start + t_end) / 2
            r = self.calc_polar(t)[0]
            if r < r_soi:
                t_start = t
            else:
                t_end = t

        return t
    
    def calc_soi_change(self, t_start):
        """Calculate the time of SOI change and the new primary
        Args:
            t0: start time of calculation
        Returns:
            t_change: time of SOI change
            new_primary: new primary body"""
        
        # if the orbit is  hyperbolic or highly eccentric, 
        # calculate time when leaving primary SOI
        if self.check_hyperbolic() or self.ra > self.primary.soi:
            t_end = t_start + self.calc_soi_leave_time()
            is_leaving = True
        # if elliptic and fully inside SOI, calc is for one orbit starting at t0
        else:
            t_end = t_start + self.T
            is_leaving = False
        
        # calculate minimum distance to all secondaries during transit time
        closest_secondary = None
        is_entering = False
        t_soi_change = t_end
        for secondary in self.primary.secondaries:
            t_min,d_min = self.calc_min_distance_to(secondary, t_start, t_end)
            if d_min < secondary.soi:
                is_entering = True
                t_soi = self.calc_soi_enter(secondary, t_start, t_min)
                # this check is for possible multiple SOI changes, to find the earliest one
                if t_soi < t_soi_change:
                    t_soi_change = t_soi
                    closest_secondary = secondary

        if is_entering:
            return t_soi_change, closest_secondary
        elif is_leaving:
            # if leaving, return time of leaving SOI and the primary of current primary
            return t_end, self.primary.orbit.primary
        else:
            # if not leaving or entering, return None
            return t_end, None

    def calc_soi_enter(self, secondary, t_start, t_min):
        """Calculate the time of entering a secondary's SOI
        Args:
            secondary: secondary body
            t_start: start time of calculation
            t_min: time of minimum distance to secondary
        Returns:
            t_soi: time of SOI change
        """
        
        t1 = t_start
        t2 = t_min
        while abs(t2-t1) > 1:
            t = (t1 + t2)/2
            if self.calc_distance_to(secondary, t) < secondary.soi:
                t2 = t
            else:
                t1 = t
        return t

    def do_maneuver(self, time, longitudinal_dv, lateral_dv=0, radial_dv=0):
        # radius at maneuver
        r = self.calc_polar(time)[0]
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
    