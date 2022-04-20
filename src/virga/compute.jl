# virga/justdoit.py and some of virga/direct_mmr_solver.py

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
function direct_solver(atm::Atmosphere, cache::VirgaCache, rmin, nrad; tol = 1e-15)

    ngas = length(cache.gases)

    pres_out = pressure
    temp_out = temperature
    z_out = LinearInterpolation(pres, z, extrapolation_bc=Flat()).(pres_out)

    mixl_out = zeros((len(pres_out), ngas))
    qc_out = zeros((len(pres_out), ngas))
    qt_out = zeros((len(pres_out), ngas))
    rg_out = zeros((len(pres_out), ngas))
    reff_out = zeros((len(pres_out), ngas))
    ndz_out = zeros((len(pres_out), ngas))
    qc_path = zeros(ngas)

    # find mmr and particle distribution for every condensible
    # perform calculation on refined TP profile but output values corresponding to initial profile
    for (i, gas) in enumerate(cache.gases)
        optical_for_layer(atm, cache, gas)
    end
        qc, qt, rg, reff, ndz, dz, qc_path[i], mixl = calc_qc(z, P_z, T_z, T_P, kz,
    gravity, gas_name, gas_mw[i], gas_mmr[i], rho_p[i], mw_atmos, mh, fsed, sig, rmin, nrad, 
    d_molecule,eps_k,c_p_factor,
    tol,og_vfall, analytical_rg)

    # generate qc values for original pressure data
    qc_out[:,i] = interp1d(pres, qc)(pres_out)
    qt_out[:,i] = interp1d(pres, qt)(pres_out)
    rg_out[:,i] = interp1d(pres, rg)(pres_out)
    reff_out[:,i] = interp1d(pres, reff)(pres_out)
    mixl_out[:,i] = interp1d(pres, mixl)(pres_out)

    ndz_temp = ndz/dz
    dz_new = insert(-(diff(z_out)), length(z_out)-1, 1e-8)
    ndz_out[:,i] = interp1d(pres, ndz_temp)(pres_out) * dz_new

    return (qc_out, qt_out, rg_out, reff_out, ndz_out, qc_path, pres_out, temp_out, z_out,mixl_out)
end
