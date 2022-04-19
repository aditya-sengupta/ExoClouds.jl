module Physics
    using ..Elements

    include("../units.jl")
    include("atmosphere.jl")
    export Atmosphere
    
    include("coagulation.jl")
    include("diffusion.jl")
    # include("heating.jl") - later
    include("mie.jl")
    export mie_efficiencies

    include("nucleation.jl")
    include("sulf_nucleation.jl")
    include("vertical.jl")

end