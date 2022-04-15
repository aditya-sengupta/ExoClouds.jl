abstract type Microphysics end

function dynamics(::Microphysics, state::State)
    println("This function should provide a contribution to the CARMA differential equation source term due to the relevant piece of physics, in a form that I've yet to work out.")
end

"""
Describes an ice or liquid element that grows in the presence of a gas condensing onto it.
"""
struct Growth <: Microphysics
    element::Element        # the growing element
    condensate::Element     # the condensing gas
end

function dynamics(g::Growth, state::State)

end