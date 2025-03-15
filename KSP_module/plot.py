import matplotlib.pyplot as plt
plt.rc("figure", figsize=[12,8])
plt.rc("grid", linestyle="--")
from KSP_module.main import orbit, planetary_body
from KSP_module.system import Kerbin

def plot_hohmann_orbit(src, dst, transfer_orbit):
    """Plot Hohmann transfer, 0: source body, 1: destination body, 2: transfer orbit """
    fig,ax = initialize_plot()
    if isinstance(src, planetary_body):
        ax = add_planetary_body_to_plot(ax, src)
    else:
        ax = add_orbit_to_plot(ax, src)

    if isinstance(dst, planetary_body):
        ax = add_planetary_body_to_plot(ax, dst)
    else:
        ax = add_orbit_to_plot(ax, dst)

    ax = add_orbit_to_plot(ax, transfer_orbit, "Transfer orbit")
    ax = add_point_to_plot(ax, transfer_orbit.calc_polar(transfer_orbit.t_launch))
    ax = add_point_to_plot(ax, transfer_orbit.calc_polar(transfer_orbit.t_arrival))
    ax = add_point_to_plot(ax, dst.orbit.calc_polar(transfer_orbit.t_launch), f"{dst} at launch")
    return fig,ax

def time(date):
    """Convert date to seconds for KSP. Input is a list of integers: [year, day, hour, minute, second]"""
    return ((date[0]-1)*Kerbin.year + (date[1]-1)*60*60*6 + date[2]*60*60 + date[3]*60 + date[4])

def date(time):
    """Convert seconds to date for KSP. Output is a list of integers: [year, day, hour, minute, second]
    Game starts at Year 1, Day 1, 0 h 0 m 0 s"""
    year = time // Kerbin.year + 1
    time = time % Kerbin.year
    day = time//(60*60*6) + 1
    time = time%(60*60*6)
    hour = time//(60*60)
    time = time%(60*60)
    minute = time//60
    second = time%60
    return [year, day, hour, minute, second]

def pretty_date(time):
    [year, day, hour, minute, second] = date(time)
    return f'Y{year:.0f} D{day:.0f} {hour:02.0f}:{minute:02.0f}:{second:02.0f}'

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
    if label is None:
        label = f'{orbit.primary.__name__} orbit'
    ax.plot(orbit.phi, orbit.r, label=label)
    ax.legend(loc=1)
    return ax

def add_point_to_plot(ax, coordinates, label=None, marker='o'):
    (r, phi) = coordinates
    ax.plot(phi, r, marker, label=label)
    ax.legend(loc=1)
    return ax
