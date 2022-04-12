# For now, this will consist of just the function signatures for every subroutine that goes into CARMA's microfast.

"""
Calculates the total amount of condensate (ice and liquid) in a particular particle type associated with a particular gas.
This is one loop iteration of CARMA/totalcondensate.F90, to be looped/mapped in a main function.

Parameters
----------
state : ...

condensate : Molecule

gas : Molecule
    The gas undergoing growth and inducing condensation.

Returns
-------
total_ice : Mass
total_liquid : Mass
"""
function condensate_quantities(state::State, condensate::Molecule, gas::Molecule)::Tuple{Mass,Mass} end

"""
Calculates supersaturation for a gas over liquid water/ice.
"""
function supersaturation(gas::Molecule, concentration::Density, T::Temperature, pvap::Pressure)
    return concentration * (R / molar_weight(gas)) * T
end