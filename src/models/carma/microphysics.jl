abstract type Microphysics end

function dynamics(::Microphysics, state::State)
    println("This function should provide a contribution to the CARMA differential equation source term due to the relevant piece of physics, in a form that I've yet to work out.")
end

struct Growth <: Microphysics
    element::Element
    condensate::Element
end

function dynamics(g::Growth, state::State)
    
end