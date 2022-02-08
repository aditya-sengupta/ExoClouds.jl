"""
Defines various molecules that may exist in atmospheres. A specific molecule is a struct that subtypes Molecule.
Note that molecules have no attributes of their own, and are instead only to be used for function dispatch.
"""

abstract type Molecule end

struct TiO₂ <: Molecule end
struct Cr <: Molecule end
struct ZnS <: Molecule end
struct NH₃ <: Molecule end
struct Na₂S <: Molecule end
struct MnS <: Molecule end
struct MgSiO₃ <: Molecule end
struct Mg₂SiO₄ <: Molecule end
struct KCl <: Molecule end
struct H₂O <: Molecule end
struct Fe <: Molecule end
struct CH₄ <: Molecule end
struct Al₂O₃ <: Molecule end
