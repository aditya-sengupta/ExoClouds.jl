
# TODO check the units of all these parameters
# Buck, 1981 (J. Atmos. Sci., 20, p. 1527)
buck = (
    BAL = 6.1121e3, 
    BBL = 18.729, 
    BCL = 257.87, 
    BDL = 227.3,
    BAI = 6.1115e3,
    BBI = 23.036, 
    BCI = 279.82,
    BDI = 333.7 
)

# (see Flatau et al., 1992, J. Appl. Meteor. p. 1507)
wexler = (
    GG0 =-0.29912729e4,
    GG1 =-0.60170128e4,
    GG2 = 0.1887643854e2,
    GG3 =-0.28354721e-1,
    GG4 = 0.17838301e-4,
    GG5 =-0.84150417e-9,
    GG6 = 0.44412543e-12,
    GG7 = 0.28584870e1,
    HH0 = -0.58653696e4,
    HH1 =  0.2224103300e2,
    HH2 =  0.13749042e-1,
    HH3 = -0.34031775e-4,
    HH4 =  0.26967687e-7,
    HH5 =  0.6918651
)

methane = (
    AMR = 16.043 / 8.3143,
	TCRIT = 90.68,
	PCRIT = .11719,
	AS = 2.213 - 2.650,
	AL = 2.213 - 3.370,
	ALS = 611.10,
	ALV = 552.36
)