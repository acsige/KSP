# List of functions in this file with short descriptions:
# - plot_hohmann_orbit: Plot a Hohmann transfer orbit, including source, destination, and transfer orbit.
# - time: Convert a date to seconds for KSP.
# - date: Convert seconds to a date for KSP.
# - pretty_date: Format a time in seconds into a readable KSP date string.
# - initialize_plot: Initialize a polar plot for orbital visualization.
# - add_planetary_orbit_to_plot: Add a planetary body's orbit to the plot.
# - add_orbit_to_plot: Add a generic orbit to the plot with an optional label.
# - add_orbit_point_to_plot: Add a specific point on an orbit to the plot.
# - calc_circle_for_polar: Calculate a circle's coordinates for a polar plot.
# - add_soi_to_plot: Add the sphere of influence (SOI) of a planetary body to the plot.
# - add_planetary_body_to_plot: Add a planetary body to the plot, including its radius.
# - add_point_to_plot: Add a point to the plot with optional label and marker.

__all__ = [
    "plot_hohmann_orbit",
    "time",
    "date",
    "pretty_date",
    "initialize_plot",
    "add_planetary_orbit_to_plot",
    "add_orbit_to_plot",
    "add_orbit_point_to_plot",
    "calc_circle_for_polar",
    "add_soi_to_plot",
    "add_planetary_body_to_plot",
    "add_point_to_plot",
]

import numpy as np
import matplotlib.pyplot as plt
plt.rc("figure", figsize=[12,8])
plt.rc("grid", linestyle="--")
from KSP_module.main import orbit, planetary_body
from KSP_module.system import Kerbin

def plot_hohmann_orbit(src, dst, transfer_orbit):
    """Plot Hohmann transfer, 0: source body, 1: destination body, 2: transfer orbit """
    fig,ax = initialize_plot()
    if isinstance(src, planetary_body):
        ax = add_planetary_orbit_to_plot(ax, src)
    else:
        ax = add_orbit_to_plot(ax, src)

    if isinstance(dst, planetary_body):
        ax = add_planetary_orbit_to_plot(ax, dst)
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

def add_planetary_orbit_to_plot(ax, body):
    ax.plot(body.orbit.phi, body.orbit.r, label=body.__name__)
    ax.legend(loc=1)
    return ax

def add_orbit_to_plot(ax, orbit, label=None):
    if label is None:
        label = f'{orbit.primary.__name__} orbit'
    ax.plot(orbit.phi, orbit.r, label=label)
    ax.legend(loc=1)
    return ax

def add_orbit_point_to_plot(ax, orbit, time, label=None, marker='o'):
    if label is None:
        label = f'{orbit.primary.__name__} orbit'
    (r, phi) = orbit.calc_polar(time)
    ax.plot(phi, r, marker, label=label)
    ax.legend(loc=1)
    return ax

def calc_circle_for_polar(r0, phi0, r):
    phi = np.linspace(0, 2*np.pi, 100)
    dx = r*np.cos(phi)
    x = dx + r0*np.cos(phi0)
    dy = r*np.sin(phi)
    y = dy + r0*np.sin(phi0)
    r_draw = np.sqrt(x**2 + y**2)
    phi_draw = np.arctan2(y, x)
    return r_draw, phi_draw

def add_soi_to_plot(ax, planetary_body, time, origin = False, label=None):
    if origin:
        r, phi = 0, 0
    else:
        (r, phi) = planetary_body.orbit.calc_polar(time)
    soi_r = planetary_body.soi
    soi_r_draw, soi_phi_draw = calc_circle_for_polar(r, phi, soi_r)
    if label: label=f'{planetary_body.__name__} SOI'
    ax.plot(soi_phi_draw, soi_r_draw, label=label)
    ax.legend(loc=1)
    return ax

def add_planetary_body_to_plot(ax, planetary_body, time, origin = False, label=None):
    if origin:
        r, phi = 0, 0
    else:
        (r, phi) = planetary_body.orbit.calc_polar(time)
    if label: label=f'{planetary_body.__name__}'
    r_body = planetary_body.radius
    r_draw, phi_draw = calc_circle_for_polar(r, phi, r_body)
    ax.plot(phi_draw, r_draw, label=label)
    ax.legend(loc=1)
    return ax

def add_point_to_plot(ax, coordinates, label=None, marker='o'):
    (r, phi) = coordinates
    ax.plot(phi, r, marker, label=label)
    ax.legend(loc=1)
    return ax
