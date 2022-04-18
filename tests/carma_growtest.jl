## This code is to test condensational growth.
## Based on carma_growtest.F90
## @author  Chuck Bardeen 
## @version May-2009

function carma_growtest()
    println("Growth Test")
  
    ## Just have one grid box. In that grid box, put an initial concentration
    ## of drops at the smallest size, then allow that to grow using a gas. The
    ## total mass of drops + gas should be conserved.

    using ExoClouds
    using Unitful
    using Unitful: Length
    using Unitful: cm

    const NZ = 1
    const NELEM = 1
    const NBIN = 24
    const NPARTICLE = 1 # NGROUP in the fortran
    const NSOLUTE = 0
    const NGAS = 1
    const NWAVE = 0
    const Δx = 100.0m
    const Δy = 100.0m
    const Δz = 100.0m
    const zmin = 3000.0m

    atm = Atmosphere(...)

    w = water()
    i = ice()
    ice_crystal = Particle([i], 1e-4cm, 2, is_ice=true)

    const processes = [
      Growth(i, w)
    ]

    # Horizonal centers
    dx = Δx * ones(NZ)
    dy = Δy * ones(NZ)
    xc, yc = dx ./ 2, dy ./ 2
    zc = zmin .+ Δz * (0.5:NZ-0.5)
    zl = zmin .+ Δz * (0:NZ-1)
    
    # Set up an arbitray mass mixing ratio of water vapor, so there is someting to
    # grow the particles.
    mmr_gas = 1e-2 .* ones(NZ, NGAS)
    mmr_cond = zeros(NZ, NELEM, NBIN)
          
    # Start with some initial water drops in the smallest bin, which can then grow
    # to larger sizes in the presence of the water vapor.
    mmr_cond[:,:,1] = 1e-6

    p[1] = 90.0hPA
    zc[1] = 17.0km
    t[1] = 190.0K
    zl[1] = zc[1] - Δz
    zl[2] = zc[1] + Δz
    rho[1] = p[1] / (gas_constant(atm) * t[1])
    pl[1] = p[1] - (zl[1] - zc[1]) * rho[1] * atm.surface_gravity
    pl[2] = p[1] - (zl[2] - zc[1]) * rho[1] * atm.surface_gravity
    mmr_gas[:,:] .= 3.5e-6
    mmr_cond[:,1,1] = ice_crystal.
    
    t_orig = t[1]
      
    ...
end
  
  