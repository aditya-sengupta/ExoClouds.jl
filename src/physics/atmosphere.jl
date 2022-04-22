using Unitful: Acceleration, Length, Mass, Temperature, Pressure # dimensions
using Unitful: cm, K, u # units
using Unitful: k, R, G # constants
using Interpolations: LinearInterpolation, Extrapolation

"""
The structural model aspects of a CARMA/virga run.
"""
struct Atmosphere
    planet_radius::FloatLeng
    surface_gravity::FloatAccl
    # these are kiiiind of treated as dynamic variables by CARMA
    # I'm going to assume they're provided as discrete functions
    # and the Extrapolation object will let us linearly interpolate on that
    # this simulation does not discretize on z just yet, but there's no other way to get P, T etc profiles
    # so we interface with these as if they're continuous, write PDEs on them, then discretize those back.

    # I don't love this because it's an abstract type and for performance I should really specify it statically
    zp::Vector{FloatLeng}
    mw::Extrapolation # Molar weight(z) 
    P::Extrapolation # Pressure(z)
    rho::Extrapolation # Atmosphere density(z) (not calling it ρ because it l.ooks too similar to p)
    # rlheat::Extrapolation # Latent heat(z)
    mmw::FloatMass
    cₚ::Float64
    mh::Float64
    ϵₖ::FloatTemp
    d_molecule::FloatLeng
    zref::FloatLeng

    function Atmosphere(
            planet_radius::FloatLeng, 
            surface_gravity::FloatAccl, 
            zp::Vector{FloatLeng}, Pp::Vector{FloatPres}, mwp::Vector{FloatMass};
            # rlheatp::Vector{TemperatureFlux{Float64}},
            cₚ::Float64=3.5, mh::Float64=1.0, ϵₖ::FloatTemp = 59.7 * K, d_molecule::FloatLeng = 2.827e-8cm, kwargs...
        )
        zref = zp[length(zp)÷2]
        mw = LinearInterpolation((zp,), mwp, extrapolation_bc=Flat())
        P =  LinearInterpolation((zp,), Pp, extrapolation_bc=Flat())
        rho_p = Pp .* (mwp ./ mol) ./ (R * Tp) # ideal gas law
        rho =  LinearInterpolation((zp,), rho_p, extrapolation_bc=Flat())
        # rlheat =  LinearInterpolation((zp,), rlheatp, extrapolation_bc=Flat())
        new(planet_radius, surface_gravity, zp, mw, P, rho, 
        # rlheat, 
        mw(zref), cₚ, mh, ϵₖ, d_molecule, zref)
    end
end

gravity(mass::Mass, radius::Length) = G * mass / (radius ^ 2)
gas_constant(atm::Atmosphere, z::Length) = R / (atm.mw(z) / mol)

function mean_free_path(atm::Atmosphere, T::Temperature, p::Pressure)
    n_atmos = p / (k * T)
    1 / (sqrt(2) * (n_atmos * π * atm.molecule_diameter^2))
end

scale_height(atm::Atmosphere, T::Temperature) = gas_constant(atm, atm.zref) * T / atm.surface_gravity

"""
Mass mixing ratio of saturated vapor of element e.
"""

function qvs(atm::Atmosphere, e::Element, T::Temperature, p::Pressure)
    denom = p * gas_constant(atm, atm.zref) / T
    (atm.supsat + 1) * vaporpressure(e, T, p) / ((R / (molecular_weight(e) / mol)) * T) / denom
end

mixing_ratio(e::Element, atm::Atmosphere) = mixing_ratio(e, atm.mmw, atm.mh)
    # atmospheric viscosity
    # EQN B2 in A & M 2001, originally from Rosner+2000
    # Rosner, D. E. 2000, Transport Processes in Chemically Reacting Flow Systems (Dover: Mineola)
viscosity(atm::Atmosphere, T::Temperature) = (5/16) * sqrt(π * k * T * (atm.mw / Na)) / (π * atm.d_molecule^2) / (1.22 * (T / atm.ϵₖ))