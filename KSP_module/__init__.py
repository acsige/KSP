# recommended usage: import KSP_module as ksp

from KSP_module.main import *
from KSP_module.system import *
from KSP_module.hohmann import *
from KSP_module.plot import *

print("KSP module loaded.")

import sys
if 'pytest' in sys.modules:
    print("Module loaded in pytest")
else:
    print("Module loaded in normal mode, running tests...")
    import pytest
    pytest.main(["KSP_module/tests", "-q"])
    print("...done.")