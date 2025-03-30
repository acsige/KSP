import KSP_module as ksp
from KSP_module import Kerbol, Kerbin, Duna, Eve, Mun, Minmus, Moho

transfer = ksp.calc_window(Kerbin.orbit, Moho.orbit, 6717200.0)

print(ksp.pretty_date(transfer.t_launch))