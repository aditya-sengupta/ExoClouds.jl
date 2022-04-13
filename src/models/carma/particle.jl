"""
A Particle is a CARMA Group, but with a better name.

CARMA has what I'll call a "group-first" structure, where groups are the smallest units you operate on, and elements are defined by their relations to groups. This is convenient when you're writing the main loop, but inconvenient if you want to modify the physics because it means you have to nest calls awkwardly even when you're asking for invariant facts about chemical elements. Put simply, calls that should have nothing to do with groups end up having to do with groups, and that makes the code awkward. It's also not very amenable to writing other models within the same physics. So what I'm doing is writing code for molecules first, writing that in its own module, and having particles know what molecules they've got within them.
"""
@with_kw struct Particle 
    elements::Vector{Molecule}
    is_ice::Bool=false
end

function vaporpressure(particle::Particle, m::Molecule, T::Temperature, p::Pressure=1*bar, mh::Real=1.0)
    if particle.is_ice
        return vaporpressure_ice(m, T, p, mh)
    else
        return vaporpressure(m, T, p, mh)
    end
end

function latent_heat(p::Particle, m::Molecule, T::Temperature)
    if p.is_ice
        rlh = latent_heat_evap(m, T) + latent_heat_melt(m, T)
    else
        rlh = latent_heat_evap(m, T)
    end
end

function akelvin(p::Particle, m::Molecule, T::Temperature)
    if p.is_ice
        return akelvin_ice(m, T)
    else
        return akelvin(m, T)
    end
end

function supersaturation(p::Particle, m::Molecule, concentration::Density, T::Temperature; relative_humidity=1, cloud_frac=1)
    ssl, ssi = supersaturations()
end