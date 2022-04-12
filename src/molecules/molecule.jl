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
struct H₂SO₄ <: Molecule end

available_molecules()::Vector{Molecule} = map(x -> x(), [TiO₂, Cr, ZnS, NH₃, Na₂S, MnS, MgSiO₃, Mg₂SiO₄, KCl, H₂O, Fe, CH₄, Al₂O₃])