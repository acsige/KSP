#%%
import numpy as np
from math import sqrt, pi
import matplotlib.pyplot as plt
plt.rc("figure", figsize=[12,8])
plt.rc("font", size=8)

class orbit:
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

def calc_hohmann(src_orbit, dst_orbit, t0):
    """
    Calculate semi-major axis, eccentricity, average angular velocity and time for Hohmann transfer.
    """
    GM = 1.1723328e18 #Kerbol grav const

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

def plot_orbits(orbit_list):
    """Plot planet orbits, 0: source, 1: destination, 2: transfer orbit """
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.clear()
    for p in orbit_list:
        ax.plot(p.phi, p.r)
    t_win = orbit_list[2].t_win
    t_arrival = t_win + orbit_list[2].T/2

    ax.plot(orbit_list[0].calc_polar(t_win)[0],orbit_list[0].calc_polar(t_win)[1], 'o', label="Launch")
    ax.plot(orbit_list[1].calc_polar(t_win)[0],orbit_list[1].calc_polar(t_win)[1], 'o', label="Dest at launch")
    ax.plot(orbit_list[1].calc_polar(t_arrival)[0],orbit_list[1].calc_polar(t_arrival)[1], 'o', label="Arrival")
    ax.plot(orbit_list[2].calc_polar(t_win)[0],orbit_list[2].calc_polar(t_win)[1], 'x', label="Hohmann launch")
    ax.plot(orbit_list[2].calc_polar(t_arrival)[0],orbit_list[2].calc_polar(t_arrival)[1], 'x', label="Hohmann arrival")
    
    ax.set_rticks([0.5, 1, 1.5, 2])  # Less radial ticks
    ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
    ax.grid(True)
    ax.legend(loc=1)

def calc_window(src_orbit, dest_orbit, t0):
    # angle difference at zero time
    d_ang_0 = dest_orbit.calc_mean_anomaly(t0) - src_orbit.calc_mean_anomaly(t0)
    
    # First iteration
    a,e,n,t_h = calc_hohmann(src_orbit, dest_orbit, t0)
    # Angle difference at ideal launch time, calculated with transfer time
    d_ang = np.pi - t_h*dest_orbit.n
    # TODO: It's only true if going in!
    # If the target is already past the position, the next opportunity must be searched
    if d_ang_0 > d_ang:
        d_ang = d_ang + 2*np.pi
    # time until next position
    # delta_ang0 + (Eve.n - Kerbin.n) * t_window = delta_ang
    t_win = t0 + (d_ang - d_ang_0) / (dest_orbit.n - src_orbit.n)
    t_arrival = t_win + t_h
    ang_win,r = src_orbit.calc_polar(t_win)
    # TODO: n is an average value, it's just for visuals so might be enough
    # If the dest is on a lower orbit, then the SV starts from the apoapsis of the transfer orbit, 
    # so an offset in its omega parameter is needed
    if src_orbit.a > dest_orbit.a:
        ang_offset = np.pi
    else:
        ang_offset = 0
    Hohmann = orbit(a,e, omega = (ang_win+ang_offset)*180/np.pi, t0=-ang_win-t_win*n)
    Hohmann.t_win = t_win

    # Real iteration
    for i in range(9):
        # calculate timing error - by how much time we've missed the target when arriving
        t_arrival = t_win + t_h
        t_miss = (Hohmann.calc_polar(t_arrival)[0] - dest_orbit.calc_polar(t_arrival)[0])/dest_orbit.n
        # modify start time using calculated error
        t_win = t_win + t_miss
        # recalculate transfer orbit
        a,e,n,t_h = calc_hohmann(src_orbit, dest_orbit, t_win)
        ang_win,r = src_orbit.calc_polar(t_win)
        Hohmann = orbit(a,e, omega = (ang_win+ang_offset)*180/np.pi, t0=-ang_win-t_win*n)
        Hohmann.t_win = t_win
        # Hohmann = orbit(a,e, omega = ang_win*180/np.pi+ang_offset, t0=-ang_win-t_win*n)
        print(t_miss)
    
    # Prepare output
    Hohmann.recalc_orbit_visu(t_win,t_arrival)
    return Hohmann
# %% 
# Testing during development
if __name__ == "__main__":
    #Moho = orbit(a=5263138304, e=0.2, omega=70)
    Eve = orbit(a=9832684544, e=0.01, omega=15)
    Kerbin = orbit(a=13599840256, e=0)
    Duna = orbit(a=20726155264, e=0.051, omega=135.5)
    
    # Hohmann = calc_window(Kerbin, Eve, 0)    
    plot_orbits([Kerbin, Eve, calc_window(Kerbin, Eve, 0)])
    # plot_orbits([Kerbin, Duna, calc_window(Kerbin, Duna, 0)])

# %%

# %%
