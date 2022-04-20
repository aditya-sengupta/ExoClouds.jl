# virga/direct_mmr_solver.py

"""
Given an atmosphere and condensates, calculate size and concentration
of condensates in balance between eddy diffusion and sedimentation.

Returns
-------
qc : ndarray 
condenstate mixing ratio (g/g)
qt : ndarray 
gas + condensate mixing ratio (g/g)
rg : ndarray
geometric mean radius of condensate  (cm) 
reff : ndarray
droplet effective radius (second moment of size distrib, cm)
ndz : ndarray 
number column density of condensate (cm^-3)
qc_path : ndarray 
vertical path of condensate 
pres : ndarray
Pressure at each layer (dyn/cm^2)
temp : ndarray
Temperature at each layer (K)
z : ndarray
altitude of each layer (cm)
"""
function direct_solver(atm::Atmosphere, cache::VirgaCache; kwargs...)

    ngas = length(cache.gases)

    pres_out = pressure
    temp_out = temperature
    z_out = LinearInterpolation(pres, z, extrapolation_bc=Flat()).(pres_out)
    nz = length(pressure)

    mixl_out = zeros(nz, ngas)
    qc_out = zeros(nz, ngas)
    qt_out = zeros(nz, ngas)
    rg_out = zeros(nz, ngas)
    reff_out = zeros(nz, ngas)
    ndz_out = zeros(nz, ngas)
    qc_path = zeros(ngas)

    # find mmr and particle distribution for every condensible
    # perform calculation on refined TP profile but output values corresponding to initial profile
    for (i, gas) in enumerate(cache.gases)
        qc, qt, rg, reff, ndz, dz, mixl, qc_path = optical_for_layer(atm, cache, gas)

        qc_out[:,i] = LinearInterpolation((pres,), qc).(pres_out)
        qt_out[:,i] = LinearInterpolation((pres,), qt).(pres_out)
        rg_out[:,i] = LinearInterpolation((pres,), rg).(pres_out)
        reff_out[:,i] = LinearInterpolation((pres,), reff).(pres_out)
        mixl_out[:,i] = LinearInterpolation((pres,), mixl).(pres_out)
        ndz_out[:,i] = LinearInterpolation(pres, ndz ./ dz).(pres_out) .* dz
    end

    return (qc_out, qt_out, rg_out, reff_out, ndz_out, qc_path, pres_out, temp_out, z_out,mixl_out)
end
