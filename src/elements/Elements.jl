"""
Defines various elements that may exist in atmospheres. A specific element is a struct that subtypes Element.
"""
module Elements
    include("../utils.jl")
    include("../utils.jl")
    include("constants.jl")
    include("element.jl")
    include("properties.jl")
    include("vaporpressure.jl")
    export Element, available_elements
    export TiO₂, Cr, ZnS, NH₃, Na₂S, MnS, MgSiO₃, Mg₂SiO₄, KCl, H₂O, Fe, CH₄, Al₂O₃, H₂SO₄
    export vaporpressure
    export molecular_weight, mixing_ratio, density
end