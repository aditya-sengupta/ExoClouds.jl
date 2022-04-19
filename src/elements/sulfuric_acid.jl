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

function density(::H₂SO₄; wtp::Float64=0.0, T::Temperature=0K)
    @warn "Make sure you passed in the correct kwargs to density for H₂SO₄! The default call will probably give something nonsensical."
    if (wtp < 0.0 || wtp > 100.0) then
        throw("Illegal value for wtp: $(wtp). Occurred at temp = $(T)")
    end
    i = findfirst(x -> wtp > x, dnwtp)
    den2 = dnc0[i] + dnc1[i] * (T / K)

    if (i == 1 || wtp == dnwtp[i]) then
        return den2
    end

    sig1 = dnc0[i-1] + dnc1[i-1] * T
    frac = (dnwtp[i] - wtp)/(dnwtp[i] - dnwtp[i-1])
    return sig1 * frac + sig2 * (1.0-frac) * g/cm^3 # TODO check this 
end
    
function surface_tension(::H₂SO₄, T::Temperature=0K; wtp::Float64=0.0)
    @warn "Make sure you passed in the correct kwargs to surface_tension for H₂SO₄! The default call will probably give something nonsensical."
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
    w = water()
    # Get number density of water (/cm3) from mass concentration (g/cm3)
    h2o_num = concentration / molecular_weight(w)

    # Get partial pressure of water (dynes/cm2) from concentration (/cm3)
    # Ideal gas law: P=nkT
    p_h2o = h2o_num * k * T
    vp_h2o = vaporpressure(w, T)

    #  Prevent a NaN calculation  
    #  In the upper thermosphere p_h2o can be very low and vp_h2o can be very high
    if (p_h2o < 1.0e-10*mbar && vp_h2o > 0*Pa) 
        p_h2o=1.e-10*mbar
    end
   
    #  Activity = water pp in mb / water eq. vp over pure water in mb
    activ = p_h2o/vp_h2o
 
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
