import numpy as np
from math import pi
import KSP_module as ksp
from KSP_module import Kerbol, Kerbin
import pytest

def test_star():
    """Test the star class"""
    assert isinstance(Kerbol, ksp.star), "Kerbol is not a star instance"
    assert Kerbol.__name__ == 'Kerbol', "Kerbol name is incorrect"
    assert np.linalg.vector_norm(Kerbol.calc_xyz(0)) == 0, "Kerbol position is not at the origin"
    # star should not have a primary body
    with pytest.raises(AttributeError):
        isinstance(Kerbol.primary, object)
    print("Star test passed")

def compare_orbits(orbit_1, orbit_2, REL_TOL=1e-9):
    """Compare two orbits to see if they are the same"""
    assert isinstance(orbit_1, ksp.orbit), "orbit1 is not an orbit instance"
    assert isinstance(orbit_2, ksp.orbit), "orbit2 is not an orbit instance"
    assert orbit_1.primary == pytest.approx(orbit_2.primary, rel=REL_TOL), "Orbits do not have the same primary body"
    assert orbit_1.a == pytest.approx(orbit_2.a, rel=REL_TOL), "Semi-major axes are not equal"
    assert orbit_1.e == pytest.approx(orbit_2.e, rel=REL_TOL), "Eccentricities are not equal"
    assert orbit_1.omega == pytest.approx(orbit_2.omega, rel=REL_TOL), "Argument of periapsis are not equal"
    assert orbit_1.i == pytest.approx(orbit_2.i, rel=REL_TOL), "Inclinations are not equal"
    assert orbit_1.OMEGA == pytest.approx(orbit_2.OMEGA, rel=REL_TOL), "Longitudes of ascending node are not equal"
    assert orbit_1.nu0 == pytest.approx(orbit_2.nu0, rel=REL_TOL), "True anomalies are not equal"
    assert orbit_1.t0 == pytest.approx(orbit_2.t0, rel=REL_TOL), "Epochs are not equal"

def test_orbit_missing_params():
    """Test that some parameters calculated during initialization are correct"""
    REL_TOL=1e-9
    #TODO: there should be a better way than copy-pasting the orbit definition
    # Test orbit is a 400 km equatorial orbit around Kerbin
    reference_orbit = ksp.orbit(Kerbin, a=1e6, e=0)
    # orbital speed calculated by the primary body
    h = 4e5
    v = reference_orbit.primary.calc_orbital_velocity(h)

    # Circular orbit periapsis and apoapsis are equal to the semi-major axis
    assert reference_orbit.ra == pytest.approx(reference_orbit.a, rel=REL_TOL), "ra is not equal to a"
    assert reference_orbit.rp == pytest.approx(reference_orbit.a, rel=REL_TOL), "rp is not equal to a"
    # Circular orbit max and min speed are the same, and equal to the orbital speed calculated by the primary body
    assert reference_orbit.va == pytest.approx(v, rel=REL_TOL), "va is not equal to v"
    assert reference_orbit.vp == pytest.approx(v, rel=REL_TOL), "vp is not equal to v"
    # Circular orbit max and min altitude are the same, and equal to the value defined earlier
    assert reference_orbit.max_alt == pytest.approx(h, rel=REL_TOL), "max_alt is not equal to h"
    assert reference_orbit.min_alt == pytest.approx(h, rel=REL_TOL), "min_alt is not equal to h"
    # orbit is elliptic
    assert reference_orbit.is_elliptic, "orbit is not elliptic"

def test_orbit_different_params():
    """Test that defining the same orbit with different orbital parameters result in the same orbit"""
    # Test orbit is a 400 km equatorial orbit around Kerbin
    reference_orbit = ksp.orbit(Kerbin, a=1e6, e=0)

    compare_orbits(ksp.orbit(Kerbin, min_alt=4e5, e=0), reference_orbit)
    compare_orbits(ksp.orbit(Kerbin, min_alt=4e5, max_alt=4e5), reference_orbit)
    compare_orbits(ksp.orbit(Kerbin, rp=1e6, ra=1e6), reference_orbit)
    compare_orbits(ksp.orbit(Kerbin, min_alt=4e5, ra=1e6), reference_orbit)
    compare_orbits(ksp.orbit(Kerbin, rp=1e6, max_alt=4e5), reference_orbit)
    compare_orbits(ksp.orbit(Kerbin, T=reference_orbit.T, e=0), reference_orbit)

def test_orbit_from_burnout():
    """Test that defining the same orbit with different orbital parameters result in the same orbit"""

    reference_orbits = []
    test_orbits = []
    # Test orbit is a 400 km orbit around Kerbin, with an inclination of 30 degrees
    reference_orbits.append(ksp.orbit(Kerbin, a=1e6, e=0, i=30, OMEGA=45))
    v_ref = reference_orbits[0].primary.calc_orbital_velocity(reference_orbits[0].min_alt)
    r_ref = 1e6
    
    # beta is the complement of the inclination
    # lambda2 is longitude of the burnout point in radians
    # Test 1: burnout point on the equator, at the ascending node
    test_orbits.append(ksp.orbit(Kerbin, v=v_ref, r=r_ref, beta=pi/3, lambda2=pi/4))

    # Test 2: burnout point at the northenmost point of the orbit
    reference_orbits.append(ksp.orbit(Kerbin, a=1e6, e=0, i=30, omega=90, OMEGA=45))
    test_orbits.append(ksp.orbit(Kerbin, v=v_ref, r=r_ref, beta=pi/2, delta = pi/6, lambda2=pi/4+pi/2))
    
    # Test 3: burnout point on the equator, at the descending node
    reference_orbits.append(ksp.orbit(Kerbin, a=1e6, e=0, i=30, omega=180, OMEGA=45))
    test_orbits.append(ksp.orbit(Kerbin, v=v_ref, r=r_ref, beta=2*pi/3, lambda2=pi/4+pi))
    
    # Test 4: burnout point at the southenmost point of the orbit
    reference_orbits.append(ksp.orbit(Kerbin, a=1e6, e=0, i=30, omega=270, OMEGA=45))
    test_orbits.append(ksp.orbit(Kerbin, v=v_ref, r=r_ref, beta=pi/2, delta = -pi/6, lambda2=3*pi/4+pi))

    # do the comparisons
    for r,t in zip(reference_orbits, test_orbits):
        compare_orbits(r, t)
