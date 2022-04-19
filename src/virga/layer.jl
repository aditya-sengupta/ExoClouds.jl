"""
Calculate layer condensate properties by iterating on optical depth
in one model layer (convering on optical depth over sublayers)
gas_name : str 
    Name of condenstante 
rho_p : float 
    density of condensed vapor (g/cm^3)
t_layer : float 
    Temperature of layer mid-pt (K)
p_layer : float 
    Pressure of layer mid-pt (dyne/cm^2)
t_top : float 
    Temperature at top of layer (K)
t_bot : float 
    Temperature at botton of layer (K)
p_top : float 
    Pressure at top of layer (dyne/cm2)
p_bot : float 
    Pressure at botton of layer 
kz : float 
    eddy diffusion coefficient (cm^2/s)
mixl : float 
    Mixing length (cm)
gas_mw : float 
    Gas molecular weight 
q_below : float 
    total mixing ratio (vapor+condensate) below layer (g/g)
supsat : float 
    Super saturation factor
z_top : float
    Altitude at top of layer
z_bot : float
    Altitude at bottom of layer
sig : float 
    Width of the log normal particle distribution 
mh : float 
    Metallicity NON log soar (1=1xSolar)
rmin : float 
    Minium radius on grid (cm)
nrad : int 
    Number of radii on Mie grid
d_molecule : float 
    diameter of atmospheric element (cm) (Rosner, 2000)
    (3.711e-8 for air, 3.798e-8 for N2, 2.827e-8 for H2)
    Set in Atmosphere constants 
Returns
-------
qc_layer : ndarray 
    condenstate mixing ratio (g/g)
qt_layer : ndarray 
    gas + condensate mixing ratio (g/g)
rg_layer : ndarray
    geometric mean radius of condensate  cm 
reff_layer : ndarray
    droplet effective radius (second moment of size distrib, cm)
ndz_layer : ndarray 
    number column density of condensate (cm^-3)
q_below : ndarray 
    total mixing ratio (vapor+condensate) below layer (g/g)
nsub : Int
    number of levels of grid refinement used
"""
function layer(atm::Atmosphere, cache::VirgaCache, gas::Element, ρₚ::Density, zs::Tuple, ps::Tuple, Ts::Tuple,
    sig; f::fsedParam=fsedConst(),
    vfallsolver::VfallSolution=OGVfall(),
    nsub_max = 128)

    z_bot, _, z_top = zs
    p_bot, p_layer, p_top, ps
    t_bot, t_layer, t_top = Ts

    r_atmos = gas_constant(atm)

    #specific gas constant for cloud
    r_cloud = R / molecular_weight(gas)

    #   pressure thickness of layer
    dp_layer = p_bot - p_top
    dlnp = log(p_bot / p_top)

    #   temperature gradient 
    dtdlnp = (t_top - t_bot) / dlnp
    #   atmospheric scale height
    scale_h = scale_height(atm, t_layer)

    #   convective velocity scale from mixing length theory
    w_convect = atm.kz / mixing_length(atm) 
    visc = viscosity(atm, t_layer)

    #   --------------------------------------------------------------------
    #   Top of convergence loop  
    nsub = 1  
    converge = false
    while !converge
        #   Zero cumulative values
        qc_layer = 0.
        qt_layer = 0.
        ndz_layer = 0.
        opd_layer = 0.        

        # total mixing ratio and pressure at bottom of sub-layer

        qt_bot_sub = q_below
        p_bot_sub = p_bot
        z_bot_sub = z_bot

        #SUBLAYER 
        dp_sub = dp_layer / nsub

        for _ in 1:nsub
            qt_below = qt_bot_sub
            p_top_sub = p_bot_sub - dp_sub
            dz_sub = scale_h * log(p_bot_sub / p_top_sub) # width of layer
            p_sub = 0.5 * (p_bot_sub + p_top_sub)
            #################### CHECK #####################
            z_top_sub = z_bot_sub + dz_sub
            z_sub = z_bot_sub + scale_h * log(p_bot_sub / p_sub) # midpoint of layer 
            ################################################
            qt_top, qc_sub, qt_sub, reff_sub, ndz_sub = optical_for_layer(atm, cache, e, qt_below, atm.fsed, mixlength,
            z_sub, ρₚ; f=f, vfallsolver=vfallsolver)

            #   vertical sums
            qc_layer = qc_layer + qc_sub*dp_sub/gravity
            qt_layer = qt_layer + qt_sub*dp_sub/gravity
            ndz_layer = ndz_layer + ndz_sub

            if reff_sub > 0.
                opd_layer = (opd_layer + 1.5*qc_sub*dp_sub/atm.surface_gravity/(ρₚ*reff_sub))
            end
    
            # Increment values at bottom of sub-layer

            qt_bot_sub = qt_top
            p_bot_sub = p_top_sub
            z_bot_sub = z_top_sub
        end

        #    Check convergence on optical depth
        if nsub_max == 1
            converge = true
        elseif nsub == 1
            opd_test = opd_layer
        elseif (opd_layer == 0.) || (nsub >= nsub_max)
            converge = true
        elseif (abs(1. - opd_test/opd_layer ) <= 1e-2)
            converge = true
        else
            opd_test = opd_layer
        end
        nsub = nsub * 2
    end # while
    #   Update properties at bottom of next layer

    q_below = qt_top

    #Get layer averages

    if opd_layer > 0.0
        reff_layer = 1.5*qc_layer / (rho_p*opd_layer)
        lnsig2 = 0.5*log( sig )^2
        rg_layer = reff_layer*exp( -5*lnsig2 )
    else 
        reff_layer = 0.
        rg_layer = 0.
    end

    qc_layer = qc_layer*gravity / dp_layer
    qt_layer = qt_layer*gravity / dp_layer

    return qc_layer, qt_layer, rg_layer, reff_layer, ndz_layer,q_below, z_cld, fsed_layer
end