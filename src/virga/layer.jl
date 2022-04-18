
using Unitful: Na, k, R
using Unitful: erg, K, g

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
gravity : float 
    Gravity of planet cgs 
mw_atmos : float 
    Molecular weight of the atmosphere 
gas_mw : float 
    Gas molecular weight 
q_below : float 
    total mixing ratio (vapor+condensate) below layer (g/g)
supsat : float 
    Super saturation factor
fsed : float
    Sedimentation efficiency coefficient (unitless) 
b : float
    Denominator of exponential in sedimentation efficiency  (if param is 'exp')
eps: float
    Minimum value of fsed function (if param=exp)
z_top : float
    Altitude at top of layer
z_bot : float
    Altitude at bottom of layer
z_alpha : float
    Altitude at which fsed=alpha for variable fsed calculation
param : str
    fsed parameterisation
    'const' (constant), 'exp' (exponential density derivation)
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
eps_k : float 
    Depth of the Lennard-Jones potential well for the atmosphere 
    Used in the viscocity calculation (units are K) (Rosner, 2000)
c_p_factor : float 
    specific heat of atmosphere (erg/K/g) . Usually 7/2 for ideal gas
    diatomic molecules (e.g. H2, N2). Technically does slowly rise with 
    increasing temperature
og_vfall : bool 
    Use original or new vfall calculation
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
function layer(gas, rho_p, t_layer, p_layer, t_top, t_bot, p_top, p_bot,
    q_below, supsat, b, z_top, z_bot, z_alpha, z_min, param,
    sig, mh, rmin, nrad, d_molecule, og_vfall, z_cld;
    nsub_max = 128, nsub = 1)

    r_atmos = gas_constant(atm)

    #specific gas constant for cloud (erg/K/g)
    r_cloud = R / molecular_weight(gas)

    #   specific heat of atmosphere (erg/K/g)
    c_p = atm.cₚ * r_atmos

    #   pressure thickness of layer
    dp_layer = p_bot - p_top
    dlnp = log(p_bot / p_top)

    #   temperature gradient 
    dtdlnp = (t_top - t_bot) / dlnp
    lapse_ratio = (t_bot - t_top) / dlnp / (t_layer / atm.cₚ)

    #   atmospheric density (g/cm^3)
    rho_atmos = p_layer / (r_atmos * t_layer)

    #   atmospheric scale height
    scale_h = r_atmos * t_layer / atm.gravity    

    #   convective velocity scale from mixing length theory
    w_convect = atm.kz / mixing_length(atm) 

    mfp = mean_free_path(atm, t_layer, p_layer)

    # atmospheric viscosity
    # EQN B2 in A & M 2001, originally from Rosner+2000
    # Rosner, D. E. 2000, Transport Processes in Chemically Reacting Flow Systems (Dover: Mineola)
    visc = (5/16 * sqrt(pi*k*t_layer*(mw_atmos/Na)) / (pi * d_molecule^2) / (1.22 * (t_layer / atm.ϵₖ)^(-0.16)))

    #   --------------------------------------------------------------------
    #   Top of convergence loop    
    converge = false
    while !converge
        #   Zero cumulative values
        qc_layer = 0.
        qt_layer = 0.
        ndz_layer = 0.
        opd_layer = 0.        

        #   total mixing ratio and pressure at bottom of sub-layer

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
            t_sub = t_bot + log(p_bot/p_sub) * dtdlnp
            qt_top, qc_sub, qt_sub, rg_sub, reff_sub,ndz_sub, z_cld, fsed_layer = calc_qc(
                    gas_name, supsat, t_sub, p_sub,r_atmos, r_cloud,
                        qt_below, dz_sub,visc,
                        rho_p,w_convect, fsed, b, eps, param, z_bot_sub, z_sub, z_alpha, z_min,
                        sig,mh, rmin, nrad, og_vfall,z_cld)


            #   vertical sums
            qc_layer = qc_layer + qc_sub*dp_sub/gravity
            qt_layer = qt_layer + qt_sub*dp_sub/gravity
            ndz_layer = ndz_layer + ndz_sub

            if reff_sub > 0.
                opd_layer = (opd_layer + 1.5*qc_sub*dp_sub/gravity/(rho_p*reff_sub))
            end
    
            #   Increment values at bottom of sub-layer

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
        elseif (abs( 1. - opd_test/opd_layer ) <= 1e-2)
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