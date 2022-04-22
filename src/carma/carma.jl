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
        temperature::FloatTemp

        function State(particle_concentration::AbstractArray, gas_concentration::AbstractArray, temperature::FloatTemp)
            @assert length(size(particle_concentration)) == 2 "particle concentration should be (nbin, nelem)"
            @assert length(size(gas_concentration)) == 1 "gas concentration should be (ngas)"
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

    function source(::Microphysics, pc::AbstractArray)
        println("This function should provide a contribution to the CARMA differential equation source term due to the relevant piece of physics, in a form that I've yet to work out.")
    end

    function loss(::Microphysics, pc::AbstractArray)
        println("This function should provide a contribution to the CARMA differential equation loss term due to the relevant piece of physics, in a form that I've yet to work out.")
    end

    """
    Describes an ice or liquid element that grows in the presence of a gas condensing onto it.
    """
    struct Nucleation <: Microphysics
        particle::Particle        # the particle being nucleated onto
        condensate::Element       # the condensing gas
        cond_index::Int64         # the state index to look up for the conc. of condensate
        are_equal::Bool           # is condensate = the nucleus element being nucleated onto?
    end

    function source!(n::Nucleation, state::State, dpc_dt::AbstractArray)
        for (i, r) in enumerate(particle.radii)
            if n.are_equal
                dpc_dt[i] += hom_nucleation(n.condensate, state.gas_concentration[n.cond_index], state.temperature)
            else
                dpc_dt[i] += het_nucleation(n.condensate, r, state.gas_concentration[n.cond_index], state.temperature)
            end
        end
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
        source1::Particle
        source2::Particle
        destination::Particle
    end

    struct VerticalTransport <: Microphysics
        particle::Particle
    end

    function gas_source() end
end