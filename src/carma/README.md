# CARMA.jl

A Julia port of the Community Aerosol and Radiation Model for Atmospheres (CARMA) bin microphysical model.

For now, this is sort of a living design document.

This implementation contains physics modules that ideally have no link to the structure of CARMA.jl. These are then incorporated into a CARMA.jl simulation via partial aliasing (arguments passed in from CARMA structs and types) in the functions `dynamics(::Microphysics, ...)`. This is done in order to write down physical laws in the way they would be written by hand/in a paper, without any reference to the overall structure and ensuring that the same physics can be employed by other simulation methods. This creates an abstraction barrier in CARMA.jl: the modeling side (`carma/physics`) ideally has no "knowledge" of the simulation side (top-level `carma`), reducing the likelihood of errors in the physics due to adapting it for specific models.

`carma/purgatory` is temporary, for files I don't think I need but might be useful later. I can also look through the commit history for these, so purgatory may or may not exist.