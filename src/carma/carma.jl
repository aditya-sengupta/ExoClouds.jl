module CARMA
    include("../units.jl")
    include("utils.jl")

    using StaticArrays: MArray

    using ..Elements
    using ..Physics

    # the CARMA state is an implicit function of t and z
    # when we discretize along z, we'll get States for each z
    mutable struct State{S <: Tuple, G <: Tuple}
        particle_concentration::MArray{S,Float64}
        gas_concentration::MArray{G,Float64}
        temperature::Temperature{Float64}

        function State(particle_concentration::AbstractArray, gas_concentration::AbstractArray, temperature::Temperature{Float64})
            new{Tuple{size(particle_concentration)...},Tuple{size(gas_concentration)...}}(
                particle_concentration, gas_concentration, temperature
            )
        end
    end

    include("particle.jl")

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