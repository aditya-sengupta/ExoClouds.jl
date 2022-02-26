using Roots

"""
Given an atmosphere and condensates, calculate size and concentration
of condensates in balance between eddy diffusion and sedimentation.
Parameters
----------
t_top : ndarray
    Temperature at each layer (K)
p_top : ndarray
    Pressure at each layer (dyn/cm^2)
t_mid : ndarray
    Temperature at each midpoint (K)
p_mid : ndarray 
    Pressure at each midpoint (dyn/cm^2)
condensibles : ndarray or list of str
    List or array of condensible gas names
b : float
    Denominator of exponential in sedimentation efficiency  (if param is 'exp')
z_top : float
    Altitude at each layer
z_alpha : float
    Altitude at which fsed=alpha for variable fsed calculation
param : str
    fsed parameterisation
    'const' (constant), 'exp' (exponential density derivation)
sig : float 
    Width of the log normal particle distribution
og_vfall : bool , optional
    optional, default = True. True does the original fall velocity calculation. 
    False does the updated one which runs a tad slower but is more consistent.
    The main effect of turning on False is particle sizes in the upper atmosphere 
    that are slightly bigger.
do_virtual : bool,optional 
    optional, Default = True which adds a virtual layer if the 
    species condenses below the model domain.
supsat : float, optional
    Default = 0 , Saturation factor (after condensation)
Returns
-------
qc : ndarray 
    condenstate mixing ratio (g/g)
qt : ndarray 
    gas + condensate mixing ratio (g/g)
rg : ndarray
    geometric mean radius of condensate  cm 
reff : ndarray
    droplet effective radius (second moment of size distrib, cm)
ndz : ndarray 
    number column density of condensate (cm^-3)
qc_path : ndarray 
    vertical path of condensate 
"""
function eddysed(
    atm::Atmosphere,
    condensibles::Vector{Molecule}
    t_top, p_top, t_mid, p_mid, 
    b, z_top, z_alpha, z_min, param,
    sig, rmin, nrad;
    og_vfall=true, do_virtual=true, supsat=0, verbose=true
)
    
    #default for everything is false, will fill in as True as we go

    did_gas_condense = falses(size(condensibles))
    t_bot = t_top[end]
    p_bot = p_top[end]
    z_bot = z_top[end]
    ngas =  length(condensibles)
    nz = length(t_mid)
    qc = zeros(nz,ngas)
    qt  = zeros(nz, ngas)
    rg = zeros(nz, ngas)
    reff = zeros(nz, ngas)
    ndz = zeros(nz, ngas)
    fsed_layer = zeros(nz,ngas)
    qc_path = zeros(ngas)
    z_cld_out = zeros(ngas)

    for (i, igas) in enumerate(condensibles)
        q_below = mixing_ratio(igas, mw_atmos, mh)

        #include decrease in condensate mixing ratio below model domain
        if do_virtual
            z_cloud = nothing
            qvs_factor = (supsat+1) * molecular_weight(igas) /mw_atmos
            pvap = vaporpressure(igas, t_bot, p_bot, mh=mh)
            qvs = qvs_factor*pvap/p_bot   
            if qvs <= q_below   
                # find the pressure at cloud base 
                # parameters for finding root 
                p_lo = p_bot
                p_hi = p_bot * 1e3

                #temperature gradient 
                dtdlnp = (t_top[end-1] - t_bot) / log(p_bot/p_top[end-1])

                #   load parameters into qvs_below common block
                qv_dtdlnp = dtdlnp
                qv_p = p_bot
                qv_t = t_bot
                qv_gas_name = igas
                qv_factor = qvs_factor
                p_base = nothing
                try
                    model = p -> qvs_below_model(p, qv_dtdlnp, qv_p, qv_t, qv_factor, qv_gas_name, mh, q_below)
                    p_base = find_zero(
                        model,
                        (p_lo, p_hi),
                        Roots.Brent()
                    )
                    if verbose 
                        println("Virtual Cloud Found: $(qv_gas_name)")
                    end
                    root_was_found = true
                catch e
                    println(e)
                    @warn "this might not be the error you want to catch"
                    root_was_found = false
                end 
                if root_was_found
                    #Yes, the gas did condense (below the grid)
                    did_gas_condense[i] = true
                    t_base = t_bot + log(p_bot/p_base)*dtdlnp
                    z_base = z_bot + scale_h[-1] * log(p_bot/p_base) 
                    
                    #   Calculate temperature and pressure below bottom layer
                    #   by adding a virtual layer 

                    p_layer_virtual = 0.5 * (p_bot + p_base)
                    t_layer_virtual = t_bot + log10(p_bot/p_layer_virtual) * dtdlnp

                    #we just need to overwrite 
                    #q_below from this output for the next routine
                    qc_v, qt_v, rg_v, reff_v,ndz_v,q_below, z_cloud, fsed_layer_v = layer(igas, rho_p[i], 
                        #t,p layers, then t.p levels below and above
                        t_layer_virtual, p_layer_virtual, t_bot,t_base, p_bot, p_base,
                        kz[-1], mixl[-1], gravity, mw_atmos, gas_mw[i], q_below,
                        supsat, fsed, b, eps, z_bot, z_base, z_alpha, z_min, param,
                        sig,mh, rmin, nrad, d_molecule,eps_k,c_p_factor, #all scalaers
                        og_vfall, z_cld
                    )
                end # if root_was_found
            end # if
        end # do_virtual
        z_cld=None
        for iz in (nz-1):-1:1 #goes from BOA to TOA
            qc[iz,i], qt[iz,i], rg[iz,i], reff[iz,i], ndz[iz,i], q_below, z_cld, fsed_layer[iz,i]  = layer( igas, rho_p[i], 
                #t,p layers, then t.p levels below and above
                t_mid[iz], p_mid[iz], t_top[iz], t_top[iz+1], p_top[iz], p_top[iz+1],
                kz[iz], mixl[iz], gravity, mw_atmos, gas_mw[i], q_below,  
                supsat, fsed, b, eps, z_top[iz], z_top[iz+1], z_alpha, z_min, param,
                sig,mh, rmin, nrad, d_molecule,eps_k,c_p_factor, #all scalars
                og_vfall, z_cld
            )

            qc_path[i] = (qc_path[i] + qc[iz,i] *
                            (p_top[iz+1] - p_top[iz]) / gravity)
        end # for
        z_cld_out[i] = z_cld

    return qc, qt, rg, reff, ndz, qc_path,mixl, z_cld_out
