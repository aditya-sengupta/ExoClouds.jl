module ExoClouds

const PROJECT_ROOT = pkgdir(ExoClouds)
# "/"*relpath((@__FILE__)*"/../..","/")

include("molecules/Molecules.jl")
include("scattering/mie.jl")

export mie_efficiencies

end # module
