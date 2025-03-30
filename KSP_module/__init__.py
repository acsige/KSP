# recommended usage: import KSP_module as ksp

from KSP_module.main import *
from KSP_module.system import *
from KSP_module.hohmann import *
from KSP_module.plot import *

print("KSP module loaded")
# Testing many things w/ Hohmann transfer orbits
max_error = 2*MISS_TOL
LKO = orbit(Kerbin, min_alt = 70000.1, e=0)
transfer1 = calc_window(Kerbin.orbit, Duna.orbit, 0)
transfer2 = calc_window(Kerbin.orbit, Eve.orbit, 0)
transfer3 = calc_window(LKO, Mun.orbit, 0)

transfer_dict = {
    # reference values calculated on 2025.03.30.
    'Duna': (transfer1, 5087190.73, 11464782.67),
    # reference values calculated on 2025.03.23.
    'Eve': (transfer2, 11824011.22, 15502404.90),
    # reference values calculated on 2025.03.23.
    'Mun': (transfer3, 1788.19, 28445.34)
}

for key, value in transfer_dict.items():
    t_launch = value[0].t_launch
    t_arrival = value[0].t_arrival
    t_launch_ref = value[1]
    t_arrival_ref = value[2]
    assert(abs(t_launch - t_launch_ref) < max_error), f'{key} t_launch: {t_launch:.2f}, ref: {t_launch_ref:.2f}'
    assert(abs(t_arrival - t_arrival_ref) < max_error), f'{key} t_arrival: {t_arrival:.2f}, ref: {t_arrival_ref:.2f}'

# Testing minimum distance calculation w/ LKO and Mun
t,d = LKO.calc_min_distance_to(Mun, 0, LKO.T)
r1,p1 = Mun.orbit.calc_polar(t)
r2,p2 = LKO.calc_polar(t)
#phase angle error less than 0.2 degrees
assert(abs(p2-p1)*180/pi < 0.2)
#distance error less than 2m
assert(abs(r1-r2-d) < 2)

#Testing that minimum distance is at periapsis
LKOe = orbit(Kerbin, min_alt = 70000.1, e=0.2, nu0=pi)
t,d = LKOe.calc_min_distance_to(Kerbin, 0, LKOe.T)
p = LKOe.calc_polar(t)[1]
ref_t = LKOe.T/2
ref_d = LKOe.min_alt + Kerbin.radius
ref_p = 0
#phase angle error less than 0.2 degrees
angle_error = (abs(p-ref_p)*180/pi)
angle_error = min(angle_error, 360-angle_error)
assert(angle_error < 0.2), f'phase angle: {p:.2f}, ref: {ref_p:.2f}'
#distance error less than 2m
assert(abs(d-ref_d) < 2), f'distance: {d:.2f}, ref: {ref_d:.2f}'
#time error less than 2s
assert(abs(t-ref_t) < 2), f'time: {t:.2f}, ref: {ref_t:.2f}'

print('All tests passed')