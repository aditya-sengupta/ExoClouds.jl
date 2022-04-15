module ExoClouds
    using Reexport

    const PROJECT_ROOT = pkgdir(ExoClouds)
    # "/"*relpath((@__FILE__)*"/../..","/")

    include("elements/Elements.jl")
    include("scattering/mie.jl")
    include("virga/virga.jl")
    include("carma/carma.jl")

    export mie_efficiencies
    @reexport using .Elements
    @reexport using .Virga
    @reexport using .CARMA

end # module
