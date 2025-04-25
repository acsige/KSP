import numpy as np
import KSP_module as ksp
from KSP_module import Kerbol, Kerbin
import pytest

def test_star():
    """Test the star class"""
    assert isinstance(Kerbol, ksp.star), "Kerbol is not a star instance"
    assert Kerbol.__name__ == 'Kerbol', "Kerbol name is incorrect"
    assert np.linalg.vector_norm(Kerbol.calc_xyz(0)) == 0, "Kerbol position is not at the origin"
    print("Star test passed")

def compare_orbits(orbit1, orbit2):
    # Inner function to compare two orbits
    assert isinstance(orbit1, ksp.orbit), "orbit1 is not an orbit instance"
    assert isinstance(orbit2, ksp.orbit), "orbit2 is not an orbit instance"
    assert orbit1.primary == pytest.approx(orbit2.primary, 1e-9), "Orbits do not have the same primary body"
    assert orbit1.a == pytest.approx(orbit2.a, 1e-9), "Semi-major axes are not equal"
    assert orbit1.e == pytest.approx(orbit2.e, 1e-9), "Eccentricities are not equal"
    assert orbit1.omega == pytest.approx(orbit2.omega, 1e-9), "Argument of periapsis are not equal"
    assert orbit1.i == pytest.approx(orbit2.i, 1e-9), "Inclinations are not equal"
    assert orbit1.OMEGA == pytest.approx(orbit2.OMEGA, 1e-9), "Longitudes of ascending node are not equal"
    assert orbit1.nu0 == pytest.approx(orbit2.nu0, 1e-9), "True anomalies are not equal"
    assert orbit1.t0 == pytest.approx(orbit2.t0, 1e-9), "Epochs are not equal"

def test_orbit():
    """Test that defining the same orbit with different orbital parameters result in the same orbit"""
    # Test orbit is a 400 km circular orbit around Kerbin
    test_orbit = ksp.orbit(Kerbin, a=1e6, e=0)

    compare_orbits(ksp.orbit(Kerbin, min_alt=4e5, e=0), test_orbit)
    compare_orbits(ksp.orbit(Kerbin, min_alt=4e5, max_alt=4e5), test_orbit)
    compare_orbits(ksp.orbit(Kerbin, rp=1e6, ra=1e6), test_orbit)
    compare_orbits(ksp.orbit(Kerbin, min_alt=4e5, ra=1e6), test_orbit)
    compare_orbits(ksp.orbit(Kerbin, rp=1e6, max_alt=4e5), test_orbit)
    compare_orbits(ksp.orbit(Kerbin, T=test_orbit.T, e=0), test_orbit)