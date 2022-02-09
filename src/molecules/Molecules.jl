"""
Defines various molecules that may exist in atmospheres. A specific molecule is a struct that subtypes Molecule.
Note that molecules have no attributes of their own, and are instead only to be used for function dispatch.
"""
module Molecules
    export TiO₂, Cr, ZnS, NH₃, Na₂S, MnS, MgSiO₃, Mg₂SiO₄, KCl, H₂O, Fe, CH₄, Al₂O₃
    export vaporpressure
    export mw, mmr, ρ
end