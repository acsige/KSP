# recommended usage: import KSP_module as ksp

from KSP_module.main import *
from KSP_module.system import *
from KSP_module.hohmann import *
from KSP_module.plot import *

print("KSP module loaded")

# Testing w/ Hohmann transfer orbits
# "good" values calculated on 2025.03.02.
max_error = 2*MISS_TOL
transfer1 = calc_window(Kerbin.orbit, Duna.orbit, 0)
# print(transfer1.t_launch, transfer1.t_arrival)
assert(abs(transfer1.t_launch - 5087925) < max_error)
assert(abs(transfer1.t_arrival - 11465643) < max_error)

transfer2 = calc_window(Kerbin.orbit, Eve.orbit, 0)
# print(transfer2.t_launch, transfer2.t_arrival)
assert(abs(transfer2.t_launch - 11823939) < max_error)
assert(abs(transfer2.t_arrival - 15502331) < max_error)

LKO = orbit(Kerbin, min_alt = 70000.1, e=0)
transfer3 = calc_window(LKO, Mun.orbit, 0)
# print(transfer3.t_launch, transfer3.t_arrival)
assert(abs(transfer3.t_launch - 1788) < max_error)
assert(abs(transfer3.t_arrival - 28445) < max_error)

# Testing minimum distance calculation w/ LKO and Mun
t,d = LKO.calc_min_distance_to(Mun, 0, LKO.T)
r1,p1 = Mun.orbit.calc_polar(t)
r2,p2 = LKO.calc_polar(t)
#phase angle error less than 0.2 degrees
assert(abs(p2-p1)*180/pi < 0.2)
#distance error less than 2m
assert(abs(r2-r1-d) < 2)

#TODO: Testing that minimum distance is at periapsis

print('All tests passed')