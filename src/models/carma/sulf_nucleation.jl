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