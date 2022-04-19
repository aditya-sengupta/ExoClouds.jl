module ExoClouds
    using Reexport

    const PROJECT_ROOT = pkgdir(ExoClouds)
    # "/"*relpath((@__FILE__)*"/../..","/")

    include("elements/Elements.jl")
    include("physics/physics.jl")
    include("virga/virga.jl")
    include("carma/carma.jl")

    @reexport using .Elements
    @reexport using .Physics
    @reexport using .Virga
    @reexport using .CARMA

end # module
