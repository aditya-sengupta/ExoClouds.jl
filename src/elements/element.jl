abstract type Element end

struct TiO₂ <: Element end
struct Cr <: Element end
struct ZnS <: Element end
struct NH₃ <: Element end
struct Na₂S <: Element end
struct MnS <: Element end
struct MgSiO₃ <: Element end
struct Mg₂SiO₄ <: Element end
struct KCl <: Element end
abstract type H₂O <: Element end
struct Water <: H₂O end
struct Ice <: H₂O end
struct Fe <: Element end
struct CH₄ <: Element end
struct Al₂O₃ <: Element end
struct H₂SO₄ <: Element end

available_elements()::Vector{Element} = map(x -> x(), [TiO₂, Cr, ZnS, NH₃, Na₂S, MnS, MgSiO₃, Mg₂SiO₄, KCl, H₂O, Fe, CH₄, Al₂O₃])