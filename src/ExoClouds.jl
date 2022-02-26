module ExoClouds
    using Reexport

    const PROJECT_ROOT = pkgdir(ExoClouds)
    # "/"*relpath((@__FILE__)*"/../..","/")

    include("molecules/Molecules.jl")
    include("scattering/mie.jl")
    include("models/virga/virga.jl")

    export mie_efficiencies
    @reexport using .Molecules
    @reexport using .Virga

end # module
