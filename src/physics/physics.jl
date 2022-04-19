module Physics
    using ..Elements

    include("../units.jl")
    include("utils.jl")
    include("atmosphere.jl")
    export Atmosphere

    include("diffusion.jl")
    include("coagulation.jl")
    # include("heating.jl") - later
    include("mie.jl")
    export mie_efficiencies

    include("nucleation.jl")
    include("sulf_nucleation.jl")
    include("vertical.jl")

end