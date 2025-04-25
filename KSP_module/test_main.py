import numpy as np
import KSP_module as ksp
from KSP_module import Kerbol

def test_star():
    # Test the star class
    assert isinstance(Kerbol, ksp.star), "Kerbol is not a star instance"
    assert Kerbol.__name__ == 'Kerbol', "Kerbol name is incorrect"
    assert np.linalg.vector_norm(Kerbol.calc_xyz(0)) == 0, "Kerbol position is not at the origin"
    print("Star test passed")