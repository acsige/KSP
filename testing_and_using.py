#%%
import matplotlib.pyplot as plt
plt.rc("figure", figsize=[12,8])
plt.rc("font", size=8)
from math import pi
import KSP_module as ksp
from KSP_module import orbit, planetary_body
from KSP_module import Kerbol, Kerbin, Duna, Eve, Mun, Minmus

# Testing during development
if __name__ == "__main__":
    # Plotting Hohmann transfer orbits
    transfer1 = ksp.calc_window(Kerbin.orbit, Duna.orbit, 0)
    fig,ax = ksp.plot_hohmann_orbit(Kerbin, Duna, transfer1)
    plt.show()
    print([ksp.pretty_date(transfer1.t_launch), transfer1.t_launch])
    
    transfer2 = ksp.calc_window(Kerbin.orbit, Eve.orbit, 0)
    fig,ax = ksp.plot_hohmann_orbit(Kerbin, Eve, transfer2)
    plt.show()
    print([ksp.pretty_date(transfer2.t_launch), transfer2.t_launch])
#%%
    LKO = orbit(Kerbin, min_alt = 70000, e=0)
    transfer3 = ksp.calc_window(LKO, Mun.orbit, 0)
    fig,ax = ksp.plot_hohmann_orbit(LKO, Mun, transfer3)
    plt.show()
    print([ksp.pretty_date(transfer3.t_launch), transfer3.t_launch])
# %%
LKO.do_maneuver(0,10)
# %%
def calc_orbit(primary, **kwargs):
    """Calculate an orbit around a primary body.
    
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
    usable_keys = ['a', 'e', 'i', 'omega', 't0', 'nu0', 'T', 'rp', 'ra', 'min_alt', 'max_alt']

    # convert altitudes to apsis distances
    if 'min_alt' in kwargs:
        assert('rp' not in kwargs)
        kwargs['rp'] = primary.radius + kwargs['min_alt']

    if 'max_alt' in kwargs:
        assert('ra' not in kwargs)
        kwargs['ra'] = primary.radius + kwargs['max_alt']

    # calculate semi-major axis from different parameters
    if 'T' in kwargs:
        assert('a' not in kwargs)
        kwargs['a'] = (kwargs['T']/(2*pi))**(2/3)*primary.GM**(1/3)

    if 'rp' in kwargs and 'ra' in kwargs:
        assert('a' not in kwargs)
        assert('e' not in kwargs)
        kwargs['a'] = (kwargs['rp'] + kwargs['ra']) / 2
        kwargs['e'] = (kwargs['ra'] - kwargs['rp']) / (kwargs['ra'] + kwargs['rp'])

    if 'rp' in kwargs and 'e' in kwargs:
        assert('a' not in kwargs)
        assert('ra' not in kwargs)
        kwargs['a'] = kwargs['rp'] / (1 - kwargs['e'])

    if 'ra' in kwargs and 'e' in kwargs:
        assert('a' not in kwargs)
        assert('rp' not in kwargs)
        kwargs['a'] = kwargs['ra'] / (1 + kwargs['e'])
        
    # at this point we should have a and e
    # set defaults for missing values
    if 'i' not in kwargs:
        kwargs['i'] = 0

    if 'omega' not in kwargs:
        kwargs['omega'] = 0

    if 't0' not in kwargs:
        kwargs['t0'] = 0

    if 'nu0' not in kwargs:
        kwargs['nu0'] = 0

    # define orbit with only the keys expected by the orbit class
    # no check for missing keys, the orbit class will raise an error
    orbit_input_keys = ['a', 'e', 'i', 'omega', 't0', 'nu0']
    orbit_kwargs = {k:kwargs[k] for k in orbit_input_keys}
    return orbit(primary, **orbit_kwargs)


lko = calc_orbit(Kerbin, min_alt = 70000, e=0)
    
# %%
