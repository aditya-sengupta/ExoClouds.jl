module CARMA
   include("utils.jl")

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

        mw::Extrapolation # Molar weight(z) 
        P::Extrapolation # Pressure(z)
        rho::Extrapolation # Atmosphere density(z) (not calling it ρ because it looks too similar to p)
        rlheat::Extrapolation # Latent heat(z)
        # metallicity::Float64
        cₚ::Float64

        function Atmosphere(Nxy::Int64, Nz::Int64, 
                planet_radius::Length{Float64}, 
                surface_gravity::Acceleration{Float64}, 
                xycoords::Horizontal, zcoords::Vertical, 
                zp::Vector{Length{Float64}}, Pp::Vector{Pressure{Float64}}, mwp::Vector{Mass{Float64}}, rlheatp::Vector{TemperatureFlux{Float64}},
                cₚ::Float64=3.5
            )
            mw = LinearInterpolation(zp, mwp)
            P = LinearInterpolation(zp, Pp)
            rho_p = Pp .* (mwp ./ mol) ./ (R * Tp) # ideal gas law
            rho = LinearInterpolation(zp, rho_p)
            rlheat = LinearInterpolation(zp, rlheatp)
            new{Nxy,Nz}(planet_radius, surface_gravity, xycoords, zcoords, mw, P, rho, rlheat, cₚ)
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

    abstract type Microphysics end

    function dynamics(::Microphysics, state::State)
        println("This function should provide a contribution to the CARMA differential equation source term due to the relevant piece of physics, in a form that I've yet to work out.")
    end

    """
    Describes an ice or liquid element that grows in the presence of a gas condensing onto it.
    """
    struct Nucleation <: Microphysics
        particle::Particle        # the particle being nucleated onto
        condensate::Element     # the condensing gas
    end

    function dynamics(n::Nucleation, state::State)
        
    end

    """
    Describes phase changes due to diffusion of condensate molecules to and away from the cloud particle.
    """
    struct Diffusion <: Microphysics
        particle::Particle
    end


    struct Heating <: Microphysics
        particle::Particle
    end

    struct Coagulation <: Microphysics
        p1::Particle
        p2::Particle
    end

    struct VerticalTransport <: Microphysics
        particle::Particle
    end
end