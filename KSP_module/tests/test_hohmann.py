import KSP_module as ksp
from KSP_module.system import *
from pytest import approx

def test_hohmann_windows():
    """
    Test Hohmann transfer windows for Duna, Eve, and Mun, comparison to reference values.
    """

    max_error = 2*ksp.MISS_TOL
    # LKO = ksp.orbit(Kerbin, min_alt = 70000.1, e=0)
    transfer1 = ksp.calc_window(Kerbin.orbit, Duna.orbit, 0)
    transfer2 = ksp.calc_window(Kerbin.orbit, Eve.orbit, 0)
    transfer3 = ksp.calc_window(LKO, Mun.orbit, 0)

    transfer_dict = {
        # reference values calculated on 2025.04.20.
        'Duna': (transfer1, 5087177, 11464768),
        # reference values calculated on 2025.04.20.
        'Eve': (transfer2, 11824001, 15502390),
        # reference values calculated on 2025.03.23.
        'Mun': (transfer3, 1788, 28445)
    }

    # compare calculated values to reference values
    for key, value in transfer_dict.items():
        transfer_orbit = value[0]
        t_launch_ref, t_arrival_ref = value[1], value[2]
        assert transfer_orbit.t_launch == approx(value[1], abs=max_error), "Launch time mismatch"
        assert transfer_orbit.t_arrival == approx(t_arrival_ref, abs=max_error), "Arrival time mismatch"