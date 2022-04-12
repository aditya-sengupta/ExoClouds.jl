using Unitful: g, cm, u
using Unitful: Mass

molecular_weight(m::Molecule) = dimless_weight(m) * u
molar_weight(m::Molecule) = dimless_weight(m) * g

function mixing_ratio(m::Molecule, mw_atmos::Mass, mh::Real, gas_mmr=nothing)
    if isnothing(gas_mmr)
        ratio = mmr_prop(m)
    else
        ratio = gas_mmr
    end
    if mh == 1
        return ratio * molecular_weight(m) / mw_atmos
    else
        throw("Alert: no M/H dependence in $(typeof(m)) routine. Consult your local theorist to determine next steps.")
    end
end

dimless_weight(::TiO₂) = 80.0
mmr_prop(::TiO₂) = 1.69e-7
density(::TiO₂) = 4.25g/cm^3
collision_diameter(::TiO₂) = 3.92e-8cm

dimless_weight(::CH₄) = 16.0
mmr_prop(::CH₄) = 4.9e-4
density(::CH₄) = 0.49g/cm^3

dimless_weight(::NH₃) = 17.0
mmr_prop(::NH₃) = 1.34e-4
density(::NH₃) = 0.84g/cm^3

dimless_weight(::H₂O) = 18.016
mmr_prop(::H₂O) = 7.54e-4
density(::H₂O) = 0.93g/cm^3 # this is water ice; liquid water is 1
collision_diameter(::H₂O) = 3.11e-8cm

dimless_weight(::Fe) = 55.845
mmr_prop(::Fe) = 5.78e-5
density(::Fe) = 7.875g/cm^3
collision_diameter(::Fe) = 4.54e-8cm

dimless_weight(::KCl) = 74.5
function mixing_ratio(m::KCl, mw_atmos::Mass, mh::Real)
    if mh == 1
        ratio = 2.2627e-7
    elseif mh == 10
        ratio = 2.1829e-6
    elseif mh == 50
        ratio = 8.1164e-6
    else
        throw("KCl gas properties can only be computed for 1, 10, and 50x solar metallicity")
    end
    ratio * molecular_weight(m) / mw_atmos
end
density(::KCl) = 1.988g/cm^3
collision_diameter(::KCl) = 3.31e-8cm

dimless_weight(::MgSiO₃) = 100.4
mixing_ratio(::MgSiO₃, mw_atmos::Mass, mh::Real) = 2.75e-3 * mh # this should be manually override-able?
density(::MgSiO₃) = 3.192g/cm^3

dimless_weight(::Mg₂SiO₄) = 140.69
mmr_prop(::Mg₂SiO₄) = 59.36e-6
density(::Mg₂SiO₄) = 3.214g/cm^3
# collision_diameter(::Mg₂SiO₄) = 6.63e-8cm

dimless_weight(::MnS) = 87.003
mmr_prop(::MnS) = 6.32e-7
density(::MnS) = 4.0g/cm^3
# collision_diameter(::MnS) = 5.22e-8cm

dimless_weight(::ZnS) = 97.474
mmr_prop(::ZnS) = 8.40e-8
density(::ZnS) = 4.04g/cm^3
# collision_diameter(::ZnS) = 2.0604e-8cm

dimless_weight(::Cr) = 51.996
function mixing_ratio(m::Cr, mw_atmos::Mass, mh::Real)
    if mh == 1
        ratio = 8.87e-7
    elseif mh == 10
        ratio = 8.6803e-6
    elseif mh == 50
        ratio = 4.1308e-5
    else
        throw("Cr gas properties can only be computed for 1, 10, and 50x solar metallicity")
    end
    ratio * molecular_weight(m) / mw_atmos
end
density(::Cr) = 7.15g/cm^3
collision_diameter(::Cr) = 4.46e-8cm

dimless_weight(::Al₂O₃) = 101.961
mmr_prop(::Al₂O₃) = 2.51e-6
density(::Al₂O₃) = 3.987g/cm^3
collision_diameter(::Al₂O₃) = 4.46e-8cm

dimless_weight(::Na₂S) = 78.0452
mmr_prop(::Na₂S) = 3.97e-6
density(::Na₂S) = 1.856g/cm^3
collision_diameter(::Na₂S) = 6.538e-8cm
