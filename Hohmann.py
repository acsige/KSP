import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, pi

class planet:
    def __init__(self, a , e, omega=0, t0 = 3.14):
        GM = 1.1723328e18 #Kerbol
        self.a = a # semi-major axis
        self.e = e # eccentricity 
        # longitude of ascending node, degrees converted to rad
        # argument of periapsis is 0, so periapsis is in the same position wrt reference angle
        self.omega = omega*np.pi/180 
        self.t0 = t0 # epoch
        self.T = sqrt(4*pi**2*a**3/GM) # orbital period, seconds
        self.n = 2*np.pi/self.T # angular velocity, rad/s
        # calculate orbit for visualization
        self.t = np.linspace(0, self.T, 100)
        self.phi, self.r = self.calc_polar(self.t)
        
    def calc_mean_anomaly(self, time):
        return time*self.n + self.t0 + self.omega

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
        beta = self.e/(1+sqrt(1-self.e*self.e))
        return E + 2*np.arctan( beta*np.sin(E) / (1-beta*np.cos(E)) )

    def calc_r(self, time):
        nu = self.calc_true_anomaly(time)
        return self.a*(1-self.e*self.e)/(1+self.e*np.cos(nu))

    # function to calculate both orbital distance and angle wrt reference angle
    def calc_polar(self, time):
        nu = self.calc_true_anomaly(time)
        phi = self.omega + nu
        r = self.a*(1-self.e*self.e)/(1+self.e*np.cos(nu))
        return phi,r

def calc_hohmann(src_planet, dst_planet, t0):
    """
    Calculate semi-major axis, eccentricity, average angular velocity and time for Hohmann transfer.
    """
    GM = 1.1723328e18 #Kerbol grav const
    
    # determine inner and outer planet
    planet_list = [src_planet, dst_planet]
    planet_list.sort(key = lambda x: x.a)
    inner_planet, outer_planet = planet_list[0], planet_list[1]

    # position of src known at t0
    phi, r_src = src_planet.calc_polar(t0)
    
    # initial guess for dest is using 0 time for Hohmann transfer
    hohmann_time = 0
    hohmann_time_prev = 10

    # Iterate a better Hohmann time
    while abs(hohmann_time_prev-hohmann_time) > 1:
        phi, r_dst = dst_planet.calc_polar(t0+hohmann_time)
        # semi-major axis
        a = (r_src + r_dst) / 2
        # Time of Hohmann transfer is half of the Hohmann orbit
        r_avg = 2*a/np.pi
        hohmann_period_sq = 4*np.pi**2*a**3/GM
        hohmann_time_prev = hohmann_time
        hohmann_time = sqrt(hohmann_period_sq) / 2    
        
    # eccentricity of Hohmann orbit
    if r_src > r_dst:
        e = 1 - 2/(r_src/r_dst + 1)
    else:
        e = 1 - 2/(r_dst/r_src + 1)

    # average angular velocity
    n = np.pi/hohmann_time
    
    return a,e,n,hohmann_time

# Testing during development
if __name__ == "__main__":
    #Moho = planet(a=5263138304, e=0.2, omega=70)
    Eve = planet(a=9832684544, e=0.01, omega=15)
    Kerbin = planet(a=13599840256, e=0)
    #Duna = planet(a=20726155264, e=0.051, omega=135.5)
    planets = [Eve, Kerbin]

    # Initial values for Kerbin -> Eve transfer calculated with this
    t_win = 0
    # angle difference at zero time
    d_ang_0 = Eve.calc_mean_anomaly(0) - Kerbin.calc_mean_anomaly(0)

    for i in range(1):
        a,e,n,t_h = calc_hohmann(Kerbin, Eve, t_win)
        # Angle difference at ideal launch time, calculated with transfer time
        d_ang = np.pi - t_h*Eve.n
        # If the target is already past the position, the next opportunity must be searched
        if d_ang_0 > d_ang:
            d_ang = d_ang + 2*np.pi
        # time until next position
        # delta_ang0 + (Eve.n - Kerbin.n) * t_window = delta_ang
        t_win = (d_ang - d_ang_0) / (Eve.n - Kerbin.n)

        #show iter results
        print(t_win/3600/6)

        ang_win,r = Kerbin.calc_polar(t_win)
        Hohmann = planet(a,e, omega = ang_win*180/np.pi+180, t0=-ang_win-t_win*n)
        planets = [Eve, Kerbin, Hohmann]

        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
        for p in planets:
            ax.plot(p.phi, p.r)
        
        t_arrival = t_win + t_h
        ax.plot(Kerbin.calc_polar(t_win)[0],Kerbin.calc_polar(t_win)[1], 'o', label="Kerbin launch")
        ax.plot(Hohmann.calc_polar(t_win)[0],Hohmann.calc_polar(t_win)[1], 'x', label="Hohmann launch")
        ax.plot(Hohmann.calc_polar(t_arrival)[0],Hohmann.calc_polar(t_arrival)[1], 'x', label="Hohmann arrival")
        ax.plot(Eve.calc_polar(t_arrival)[0],Eve.calc_polar(t_arrival)[1], 'o', label="Eve arrival")
        
        ax.set_rticks([0.5, 1, 1.5, 2])  # Less radial ticks
        ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
        ax.grid(True)
        ax.legend(loc=1)