stwtp = [0, 23.8141, 38.0279, 40.6856, 45.335, 52.9305, 56.2735,
    59.8557, 66.2364, 73.103, 79.432, 85.9195, 91.7444, 97.6687, 100]
    
stc0 = [117.564, 103.303, 101.796, 100.42, 98.4993, 91.8866, 
    88.3033, 86.5546, 84.471, 81.2939, 79.3556, 75.608, 70.0777,
    63.7412, 61.4591]

stc1 = [-0.153641, -0.0982007, -0.0872379, -0.0818509,           
    -0.0746702, -0.0522399, -0.0407773, -0.0357946, -0.0317062,   
    -0.025825, -0.0267212, -0.0269204, -0.0276187, -0.0302094,    
    -0.0303081]

dnwtp = [0., 1., 5., 10., 20., 25., 30., 35., 40.,
41., 45., 50., 53., 55., 56., 60., 65., 66., 70.,
72., 73., 74., 75., 76., 78., 79., 80., 81., 82.,
83., 84., 85., 86., 87., 88., 89., 90., 91., 92.,
93., 94., 95., 96., 97., 98., 100.]

dnc0 = [1, 1.13185, 1.17171, 1.22164, 1.3219, 1.37209,
1.42185, 1.4705, 1.51767, 1.52731, 1.56584, 1.61834, 1.65191,
1.6752, 1.68708, 1.7356, 1.7997, 1.81271, 1.86696, 1.89491,
1.9092, 1.92395, 1.93904, 1.95438, 1.98574, 2.00151, 2.01703,
2.03234, 2.04716, 2.06082, 2.07363, 2.08461, 2.09386, 2.10143,
2.10764, 2.11283, 2.11671, 2.11938, 2.12125, 2.1219, 2.12723, 
2.12654, 2.12621, 2.12561, 2.12494, 2.12093]

dnc1 = [0,  -0.000435022, -0.000479481, -0.000531558, -0.000622448,
-0.000660866, -0.000693492, -0.000718251, -0.000732869, -0.000735755, 
-0.000744294, -0.000761493, -0.000774238, -0.00078392, -0.000788939,  
-0.00080946, -0.000839848, -0.000845825, -0.000874337, -0.000890074,  
-0.00089873, -0.000908778, -0.000920012, -0.000932184, -0.000959514,  
-0.000974043, -0.000988264, -0.00100258, -0.00101634, -0.00102762,    
-0.00103757, -0.00104337, -0.00104563, -0.00104458, -0.00104144,      
-0.00103719, -0.00103089, -0.00102262, -0.00101355, -0.00100249,      
-0.00100934, -0.000998299, -0.000990961, -0.000985845, -0.000984529,  
-0.000989315]

water, acid = H₂O(), H₂SO₄()

function sulfate_density(wtp::Real, T::Temperature)
    if (wtp < 0.0 || wtp > 100.0) then
        throw("Illegal value for wtp: $(wtp). Occurred at temp = $(T)")
    end
    i = findfirst(x -> wtp > x, dnwtp)
    den2 = dnc0[i] + dnc1[i] * (T / K)

    if (i == 1 || wtp == dnwtp[i]) then
        return den2
    end

    sig1 = dnc0[i-1] + dnc1[i-1] * temp
    frac = (dnwtp[i] - wtp)/(dnwtp[i] - dnwtp[i-1])
    return sig1 * frac + sig2 * (1.0-frac) * g/cm^3 # TODO check this 
end
    
function sulfate_surf_tens(wtp::Real, T::Temperature)
    if (wtp < 0.0 || wtp > 100.0) then
        throw("Illegal value for wtp: $(wtp). Occurred at temp = $(T)")
    end
      
    i = findfirst(x -> wtp > x, stwtp)
    sig2 = stc0[i] + stc1[i] * (T / K)
  
    if (i == 1 || wtp == stwtp[i]) then
      return sig2
    end
  
    sig1 = stc0[i-1] + stc1[i-1] * temp
    frac = (stwtp[i] - wtp)/(stwtp[i] - stwtp[i-1])
    return sig1 * frac + sig2 * (1.0-frac) * erg/cm^2
end

function wtpct_tabaz(T::Temperature, concentration::Density)
    # Get number density of water (/cm3) from mass concentration (g/cm3)
    h2o_num = concentration / molecular_weight(water)

    # Get partial pressure of water (dynes/cm2) from concentration (/cm3)
    # Ideal gas law: P=nkT
    p_h2o = h2o_num * k * T
    vp_h2o = vaporpressure(water, T)

    #  Prevent a NaN calculation  
    #  In the upper thermosphere p_h2o can be very low and vp_h2o can be very high
    if (p_h2o < 1.0e-10*mbar && vp_h2o > 0*Pa) 
        p_h2o=1.e-10*mbar
    end
   
    #  Activity = water pp in mb / water eq. vp over pure water in mb
    activ = p_h2o/vp_h2o
!	write(*,*) activ,p_h2o,vp_h2o
 
    if (activ < 0.05)
      activ = max(activ, 1.e-6)    # restrict minimum activity
      atab1 	= 12.37208932	
      btab1 	= -0.16125516114
      ctab1 	= -30.490657554
      dtab1 	= -2.1133114241
      atab2 	= 13.455394705	
      btab2 	= -0.1921312255
      ctab2 	= -34.285174607
      dtab2 	= -1.7620073078
    elseif (activ > 0.05 && activ < 0.85)
      atab1 	= 11.820654354
      btab1 	= -0.20786404244
      ctab1 	= -4.807306373
      dtab1 	= -5.1727540348
      atab2 	= 12.891938068	
      btab2 	= -0.23233847708
      ctab2 	= -6.4261237757
      dtab2 	= -4.9005471319
    elseif (activ > 0.85)
      activ = min(activ, 1.)      # restrict maximum activity
      atab1 	= -180.06541028
      btab1 	= -0.38601102592
      ctab1 	= -93.317846778
      dtab1 	= 273.88132245
      atab2 	= -176.95814097
      btab2 	= -0.36257048154
      ctab2 	= -90.469744201
      dtab2 	= 267.45509988
    end

    contl = atab1*(activ^btab1)+ctab1*activ+dtab1
    conth = atab2*(activ^btab2)+ctab2*activ+dtab2
    
    contt = contl + (conth-contl) * ((temp -190.)/70.)
    conwtp = (contt*98.) + 1000.

    wtpct_tabaz = (100. * contt * 98.)/conwtp
    wtpct_tabaz = min(max(wtpct_tabaz,1.),100.) # restrict between 1 and 100 %
      
    #  Note: restricting activity to 1.e-6 minimum allows for a maximum of
    #  98.5 wtpct at T=650K, 95.8 wtpct at T=300K, and 90.9 wtpct at 180K.
  
    return wtpct_tabaz
end



function sulfuric_nucleation_rate(conc_water::Density, conc_acid::Density, T::Temperature, max_radius::Length)::Tuple{Length,SpecificParticleRate}
    dnpot = zeros(46)
    dnwf = dnwtp ./ 100
    dnpot = 4.184 .* (23624.8 .- 1.14208e8 ./ ((dnwtp .- 105.318) .^ 2 .+ 4798.69))

    # Molecular weight ratio of H2SO4 / H2O:
    wtmolr = molar_weight(acid) / molar_weight(water)

    # Compute H2O and H2SO4 concentrations in molecules/cm3 
    # (factor of Avogadro's constant removed because Unitful should handle it, as we're using molecular weight)
    conc_water_mol   = conc_water / molecular_weight(water)
    conc_acid_mol    = conc_acid  / molecular_weight(acid)

    # Compute ln of H2O and H2SO4 ambient vapor pressures [dynes/cm2]
    h2oln   = log(h2o_cgs   * (R / molar_weight(water)) * T)
    h2so4ln = log(h2so4_cgs * (R / molar_weight(acid)) * T)

    # loop through wt pcts and calculate vp/composition for each
    dens = dnc0 .+ dnc1 .* T

    # Calc. water eq.vp over solution using (Lin & Tabazadeh eqn 5, JGR, 2001)  
    cw = 22.7490 .+ 0.0424817 .* dnwtp .- 0.0567432 * dnwtp .^ 0.5 .- 0.000621533 .* dnwtp .^ 2
    dw = -5850.24 .+ 21.9744 .* dnwtp .- 44.5210 .* dnwtp .^ 0.5 .- 0.384362 .* dnwtp .^ 2
    # pH20 | eq[mb]
    wvp = exp(cw .+ dw ./ T)
    # Ln(pH2O | eq [dynes/cm2])
    wvpln = log(wvp * 1000)
    # Save the water eq.vp over solution at each wt pct into this array:
    # Ln(pH2O/pH2O|eq) with both terms in dynes/cm2
    pb = h2oln .- wvpln
    # Calc. sulfuric acid eq.vp over solution using (Ayers et. al., GRL, V.7, No.6, June 1980)
    #
    # T0 set in the low end of the Ayers measurement range (338-445K)
    t0_kulm = 340K
    seqln   = -10156 / t0_kulm + 16.259
    # Now calc. Kulmala correction (J. CHEM. PHYS. V.93, No.1, 1 July 1990)
    #
    # Critical temperature = 1.5 * Boiling point
    t_crit_kulm = 905K
    factor_kulm = -1 / T + 1 / t0_kulm + 0.38 / (t_crit_kulm - t0_kulm) * (1.0 + log(t0_kulm / T) - t0_kulm / T)

    # For pure sulfuric acid
    seqln = seqln + 10156 * factor_kulm

    # Now adjust vp based on weight % composition using parameterization of Giauque 1960
    #
    # Adjust for WTPCT composition
    seqln = seqln .- dnpot / (R * T)

    # Convert atmospheres => dynes/cm2
    seqln = seqln + log(1013250)

    # Save the sulfuric acid eq.vp over solution at each wt pct into this array:
    #
    # Ln(pH2SO4/pH2SO4|eq) with both terms in dynes/cm2
    pa = h2so4ln .- seqln

    # Create 2-component solutions of varying composition c1 and c2
    c1 = pa. - pb .* wtmolr
    c2 = pa .* dnwf .+ pb * (1 .- dnwf) .* wtmolr

    # Now loop through until we find the c1+c2 combination with minimum Gibbs free energy
    fct = zeros(46)
    dw2     = dnwtp[46] - dnwtp[45]
    dens1   = (dens[46] - dens[45]) / dw2
    fct[46] = c1[46] + c2[46] * 100 * dens1 / dens[46]
    dens12 = dens1
        
    i = 45 # in Julia loops are LLVMd away anyway
    while (i > 1)
        dw1    = dw2
        dens11 = dens12
        dw2    = dnwtp[i] - dnwtp[i-1]
        dens12 = (dens[i] - dens[i-1]) / dw2
        dens1  = (dens11 * dw2 + dens12 * dw1) / (dw1 + dw2)

        fct[i] = c1[i] + c2[i] * 100 * dens1 / dens[i]

        # Find saddle where fct[i]<0<fct[i+1]
        if (fct[i] * fct[i+1] <= 0) 
            break
        end
    end
    
    if (i == 1)
        dens1  = (dens[2] - dens[1]) / (dnwtp[2] - dnwtp[1])
        fct[1] = c1[1] + c2[1] * 100 * dens1 / dens[1]
    end
        
    # Possibility 1: loop finds no saddle, so no nucleation occurs:
    if (fct[i] * fct[i+1] > 0) then
        return 0.0cm, 0.0/cm^3/s
    # Possibility 2: loop crossed the saddle; interpolate to find exact value:
    elseif (fct[i] * fct[i+1] < 0)
        xfrac = fct[i+1] / (fct[i+1] - fct[i])
        wstar = dnwtp[i+1] * (1.0 - xfrac) + dnwtp[i] * xfrac
        dstar = dens[i+1]  * (1.0 - xfrac) + dens[i]  * xfrac
        rhln  = pb[i+1] * (1.0 - xfrac) + pb[i] * xfrac
        raln  = pa[i+1] * (1.0 - xfrac) + pa[i] * xfrac
    # throwing out Possibility 3: loop found the saddle point exactly, because that should just be an interpolation with (1, 0)
    end

    # Critical weight fraction
    wfstar = wstar / 100

    if ((wfstar < 0) || (wfstar > 1)) then
        throw("sulfnuc: wstar out of bounds!")
    end
        
    # Critical surface tension  [erg/cm2]
    sigma = sulfate_surf_tens(wstar, T)  
        
    # Critical Y (eqn 13 in Zhao & Turco 1993) [erg/cm3]
    ystar = dstar * R * T * (wfstar / molar_weight(acid) * raln + (1 - wfstar) / molar_weight(water) * rhln)
    if (ystar < 1.e-20 * erg/cm^3) then
        return 0.0cm, 0.0/cm^3/s
    end

    # Critical cluster radius [cm]        
    rstar = 2 * sigma / ystar 
    rstar = max(rstar, 0.0cm)
    r2    = rstar * rstar

    # Critical Gibbs free energy [erg]
    gstar = (4π / 3) * r2 * sigma 
    
    #   kT/(2*Pi*M) = [erg/mol/K]*[K]/[g/mol] = [erg/g] = [cm2/s2]
    #   RB[erg/mol] = RGAS[erg/mol/K] * T[K] / (2Pi)
    rb = R * T / (2π)
        
    # Beta[cm/s] = sqrt(RB[erg/mol] / WTMOL[g/mol])
    beta1 = sqrt(rb / molar_weight(acid)) 
    beta2 = sqrt(rb / molar_weight(water))

    # RPR[molecules/s] = 4Pi * R2[cm2] * H2O[molecules/cm3] * Beta[cm/s]
    rpr = 4π * r2 * conc_water_mol * beta1

    # RPRE[/cm3/s] = RPR[/s] * H2SO4[/cm3]; first part of Zhao & Turco eqn 16
    rpre = rpr * conc_acid_mol

    # Zeldovitch non-equilibrium correction factor [unitless]
    # Jaecker-Voirol & Mirabel, 1988 (not considered in Zhao & Turco) 
    fracmol = 1 /(1 + wtmolr * (1 - wfstar) / wfstar)
    zphi    = atan(fracmol)
    zeld    = 0.25 / (sin(zphi))^2

    # Empirical correction factor:
    cfac = 0.0

    # Gstar exponential term in Zhao & Turco eqn 16 [unitless]
    ahom = (-gstar / BK / T) + cfac
    if (ahom < -500) then
        exhom = 0.0
    else
        # exhom = exp(min(ahom, 28.0))
        exhom = exp(ahom)
    end

    #   Calculate mass of critical nucleus
    rho_H2SO4_wet = sulfate_density(wtpct, T)
    rmstar = (4π / 3) * rho_H2SO4_wet * r2 * rstar
        
    # Calculate dry mass of critical nucleus
    nuc_radius = rmstar * wfstar

    # If none of the bins are large enough for the critical radius, then
    # no nucleation will occur.
    if (nuc_radius > max_radius)
        return 0.0cm, 0.0/cm^3/s
    else
        # Calculate the nucleation rate [#/cm3/s], Zhao & Turco eqn 16.
        nuc_rate = rpre * zeld * exhom #* redugrow(igash2so4)
        # literally every reference I can find to redugrow either calls it or the one exception sets it to 1
        # so I am leaving it out as I don't know what it does
        # and it seemingly doesn't matter
    end
    
  return nuc_radius, nuc_rate
end