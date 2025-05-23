from KSP_module.main import star, planetary_body, orbit

Kerbol = star('Kerbol')

# Planetary bodies
Moho = planetary_body('Moho', orbit(Kerbol, a=5263138304, e=0.2, omega=15, i=7, nu0=3.14),
                    GM=8.1717302e12, radius=7e5, atmo_height=9e4)

Eve = planetary_body('Eve', orbit(Kerbol, a=9832684544, e=0.01, omega=0, OMEGA=15, i=2.1, nu0=3.14),
                    GM=8.1717302e12, radius=7e5, atmo_height=9e4)
Gilly = planetary_body('Gilly', orbit(Eve, a=3.15e7, e=0.55, omega=10, OMEGA=80, i=12, nu0=0.9),
                    GM=8289449.8, radius=1.3e4, atmo_height=0)

Kerbin = planetary_body('Kerbin', orbit(Kerbol, a=13599840256, e=0, nu0=3.14),
                        GM=3.5316e12, radius=6e5, atmo_height=7e4)
Mun = planetary_body('Mun', orbit(Kerbin, a=1.2e7, e=0, nu0=1.7),
                        GM=6.5138398e10, radius=2e5, atmo_height=0)
Minmus = planetary_body('Minmus', orbit(Kerbin, a=4.7e7, e=0, omega=38, OMEGA=78, nu0=0.9),
                        GM=1.7658e9, radius=6e4, atmo_height=0)

Duna = planetary_body('Duna', orbit(Kerbol, a=20726155264, e=0.051, omega=0, OMEGA=135.5, i=0.06, nu0=3.14),
                        GM=3.0136321e11, radius=3.2e5, atmo_height=5e4)
Ike = planetary_body('Ike', orbit(Duna, a=3.2e6, e=0.03, i=0.2, nu0=1.7),
                        GM=1.8568369e10, radius=1.3e5, atmo_height=0)

Dres = planetary_body('Dres', orbit(Kerbol, a=40839348203, e=0.145, omega=90, OMEGA=280, i=5.0, nu0=3.14),
                        GM=2.1484489e10, radius=1.38e5, atmo_height=0)

Jool = planetary_body('Jool', orbit(Kerbol, a=68773560320, e=0.05, omega=0, OMEGA=52, i=1.304, nu0=0.1),
                        GM=2.82528e14, radius=6e6, atmo_height=2e5)

# Important orbits
# LKO = Low Kerbin Orbit
LKO = orbit(Kerbin, min_alt=70000.1, e=0)
# LMO = Low Mun Orbit
LMO = orbit(Mun, min_alt=15000, e=0)

# LDO = Low Duna Orbit
LDO = orbit(Duna, min_alt=50000.1, e=0)
# LIO = Low Ike Orbit
LIO = orbit(Ike, min_alt=15000, e=0)

# LEO = Low Eve Orbit
LEO = orbit(Eve, min_alt=90000.1, e=0)