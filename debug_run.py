import KSP_module as ksp
from KSP_module import Kerbol, Kerbin, Duna, Eve, Mun, Minmus, Moho, Jool, Dres

print(ksp.time([2, 1, 0, 0, 0]))  # Y2 D1 00:00:00
# transfer = ksp.calc_window(Kerbin.orbit, Jool.orbit, 4139600)
transfer = ksp.calc_window(Kerbin.orbit, Dres.orbit, ksp.time([2, 1, 0, 0, 0]))

print(f'{ksp.pretty_date(transfer.t_launch)}, {transfer.t_launch}')