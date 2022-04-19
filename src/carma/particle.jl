"""
A Particle is a CARMA Group, but with a better name.

CARMA has what I'll call a "group-first" structure, where groups are the smallest units you operate on, and elements are defined by their relations to groups. This is convenient when you're writing the main loop, but inconvenient if you want to modify the physics because it means you have to nest calls awkwardly even when you're asking for invariant facts about chemical elements. Put simply, calls that should have nothing to do with groups end up having to do with groups, and that makes the code awkward. It's also not very amenable to writing other models within the same physics. So what I'm doing is writing code for elements first, writing that in its own module, and having particles know what elements they've got within them.

I've realized the difference here is I'm thinking about the problem in a Lagrangian (particle-first) way, and CARMA is an Eulerian (field-first) solver. I think Lagrangian is probably the better way for reasons I wrote in the essay, but I'll try not to rock the boat too much just yet.
"""
@with_kw struct Particle 
    # elements::Vector{Element} 
    # you always have to pass this in anyway to specify which one you're talking about
    radius_min::Length
    radius_ratio::Float64
    aspect_ratio::Float64=1 # eshape
    is_ice::Bool=false
    ref_idx::Float64=1.0
    is_cloud::Bool=true
end

function radius_to_bin(particle::Particle, radius::Length)
    if radius < particle.radius_min
        return 1
    else
        return 2 + Int(floor(log(radius / particle.radius_min) / log(particle.radius_ratio)))
    end
end
