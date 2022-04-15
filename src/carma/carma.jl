module CARMA

    using Parameters
    using Unitful: Length, Acceleration, Mass, Pressure, Temperature
    using Interpolations: LinearInterpolation, Extrapolation

    @derived_dimension TemperatureFlux 𝚯*𝐓^-1 true

    include("coordinates.jl")

    """
    The structural model aspects of a CARMA run.
    """
    struct Atmosphere{Nxy,Nz} # unify this with Virga's at some point
        planet_radius::Length
        surface_gravity::Acceleration
        xy::Horizontal{Nxy}
        z::Vertical{Nz}
        # these are kiiiind of treated as dynamic variables by CARMA
        # I'm going to assume they're provided as discrete functions
        # and the Extrapolation object will let us linearly interpolate on that
        # this simulation does not discretize on z just yet, but there's no other way to get P, T etc profiles
        # so we interface with these as if they're continuous, write PDEs on them, then discretize those back.
        # p::SVector{Nz,Pressure{Float64}}
        # T::SVector{Nz,Temperature{Float64}} # T(z, t)
        P::Extrapolation # Pressure(z)
        rhoₐ::Extrapolation # Atmosphere density(z)
        rlheat::Extrapolation # Latent heat(z)
        mw_atm::Mass{Float64}
        # metallicity::Float64
        cₚ::Float64

        function Atmosphere(Nxy::Int64, Nz::Int64, 
                planet_radius::Length{Float64}, 
                surface_gravity::Acceleration{Float64}, 
                xycoords::Horizontal, zcoords::Vertical, 
                zp::Vector{Length{Float64}}, pp::Vector{Pressure{Float64}}, rho_p::Vector{Density{Float64}}, rlheatp::Vector{TemperatureFlux{Float64}},
                mw_atm::Mass{Float64}=2.2u, cₚ::Float64=3.5
            )
            P = LinearInterpolation(zp, pp)
            rhoₐ = LinearInterpolation(zp, rho_p)
            rlheat = LinearInterpolation(zp, rlheatp)
            new{Nxy,Nz}(planet_radius, surface_gravity, xycoords, zcoords, P, rhoₐ, rlheat, mw_atm, cₚ)
        end
    end

    # the CARMA state is an implicit function of t and z
    # when we discretize along z, we'll get States for each z
    mutable struct State{NBIN,NELEM,NGAS}
        particle_concentration::MArray{(NBIN,NELEM),Float64}
        gas_concentration::MArray{(NGAS),Float64}
        temperature::Float64
    end

    ## I think these are just constants for everything I'm doing
    # cloud_fraction::Float64
    # relative_humidity::Float64

end