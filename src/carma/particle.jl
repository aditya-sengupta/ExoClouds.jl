"""
A Particle is a CARMA Group, but with a better name.

CARMA has what I'll call a "group-first" structure, where groups are the smallest units you operate on, and elements are defined by their relations to groups. This is convenient when you're writing the main loop, but inconvenient if you want to modify the physics because it means you have to nest calls awkwardly even when you're asking for invariant facts about chemical elements. Put simply, calls that should have nothing to do with groups end up having to do with groups, and that makes the code awkward. It's also not very amenable to writing other models within the same physics. So what I'm doing is writing code for elements first, writing that in its own module, and having particles know what elements they've got within them.

I've realized the difference here is I'm thinking about the problem in a Lagrangian (particle-first) way, and CARMA is an Eulerian (field-first) solver. I think Lagrangian is probably the better way for reasons I wrote in the essay, but I'll try not to rock the boat too much just yet.
"""
struct Particle 
    # elements::Vector{Element} 
    # you always have to pass this in anyway to specify which one you're talking about
    radius_min::FloatLeng
    mass_ratio::Float64
    radii::Vector{FloatLeng}
    aspect_ratio::Float64 # eshape
    is_ice::Bool
    ref_idx::Float64
    is_cloud::Bool

    function Particle(
        # density is tmp_rhop in Fortran
        radius_min::FloatLeng, mass_ratio::Float64, density::Density, nbin::Int64, is_ice::Bool=false; aspect_ratio::Float64=1.0, ref_idx::Float64=1.0, is_cloud::Bool=true
    )
        rmassmin = (4/3 * π) * density * radius_min^3
        masses = rmassmin .* (mass_ratio) .^ (0:nbin-1)
        radii = (masses ./ (density / (4/3 * π))) .^ (1/3)
        new(radius_min, mass_ratio, radii, aspect_ratio, is_ice, ref_idx, is_cloud)
    end
end

function radius_to_bin(particle::Particle, radius::Length)
    i = findfirst(radius .< particle.radii)
    if isnothing(i)
        return length(particle.radii)
    else
        return i - 1
    end
end
