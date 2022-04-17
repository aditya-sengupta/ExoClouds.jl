## This code is to test condensational growth.
##
## Upon execution, a text file (carma_growtest.txt) is generated.
## The text file can be read with the IDL procedure read_growtest.pro.
##
## @author  Chuck Bardeen
## @version May-2009

function carma_growtest()
    println("Growth Test")
  
    ## Just have one grid box. In that grid box, put an initial concentration
    ## of drops at the smallest size, then allow that to grow using a gas. The
    ## total mass of drops + gas should be conserved.

    using ExoClouds
    using Unitful
    using Unitful: Length
    using Unitful: cm

    const NZ = 1
    const NZP1 = NZ + 1
    const NELEM = 1
    const NBIN = 24
    const NPARTICLE = 1 # NGROUP in the fortran
    const NSOLUTE = 0
    const NGAS = 1
    const NWAVE = 0
    const Δx = 100.0m
    const Δy = 100.0m
    const Δz = 100.0m
    const zmin = 3000.0m

    w = water()
    i = ice()
    ice_crystal = Particle([i], 1e-4cm, 2, is_ice=true)

    const processes = [
      Growth(i, w)
    ]

    # Horizonal centers
    dx = Δx * ones(NZ)
    dy = Δy * ones(NZ)
    xc, yc = dx ./ 2, dy ./ 2
    zc = zmin .+ Δz * (0.5:NZ-0.5)
    zl = zmin .+ Δz * (0:NZ-1)
    
    # Set up an arbitray mass mixing ratio of water vapor, so there is someting to
    # grow the particles.
    mmr_gas(:,:) = 1e-2_f
          
    # Start with some initial water drops in the smallest bin, which can then grow
    # to larger sizes in the presence of the water vapor.
    mmr(:,:,1) = 1e-6_f

    # Try TTL Conditions ...
    #
    # p, T, z, mmrgas, rmin, particle concentration
    # 90 hPa, 190 K, 17 km, 3.5e-6 g/g, 5 um, 0.1/cm^3.
    p[1]         = 90hPA
    zc[1]        = 17km
    t[1]         = 190K
    zl(1)        = zc(1) - Δz
    zl(2)        = zc(1) + Δz
    rho(1)       = (p(1) * 10._f) / (R_AIR * t(1)) * (1e-3_f * 1e6_f)
    pl(1)        = p(1) - (zl(1) - zc(1)) * rho(1) * (GRAV / 100._f)
    pl(2)        = p(1) - (zl(2) - zc(1)) * rho(1) * (GRAV / 100._f)
    mmr_gas(:,:) = 3.5e-6_f
    mmr(:,1,1)   = (0.1_f * rmass(1) * (1e-3_f * 1e6_f)) / rho(:)
    
    t_orig = t(1)
      
    
    write(lun,*) 0

    write(lun,*) 0._f, 0._f

    do ielem = 1, NELEM
      do ibin = 1, NBIN
        write(lun,'(2i4,e10.3)') ielem, ibin, real(mmr(1,ielem,ibin))
      end do
    end do

    do igas = 1, NGAS
      write(lun,'(i4,3e10.3)') igas, real(mmr_gas(1,igas)), 0., 0.
    end do

      
    # Iterate the model over a few time steps.
    do istep = 1, nstep
    
      # Calculate the model time.
      time = (istep - 1) * dtime

        # Create a CARMASTATE for this column.
        call CARMASTATE_Create(cstate, carma_ptr, time, dtime, NZ, &
                            I_CART, I_CART, lat, lon, &
                            xc(:), dx(:), &
                            yc(:), dy(:), &
                            zc(:), zl(:), &
                            p(:),  pl(:), &
                            t(:), rc)
      if (rc /=0) stop "    *** CARMASTATE_Create FAILED ***"
      
        # Send the bin mmrs to CARMA
        do ielem = 1, NELEM
          do ibin = 1, NBIN
          call CARMASTATE_SetBin(cstate, ielem, ibin, mmr(:,ielem,ibin), rc)
      if (rc /=0) stop "    *** CARMASTATE_SetBin FAILED ***"
          end do
        end do
        
      # Send the gas mmrs to CARMA
      do igas = 1, NGAS
        call CARMASTATE_SetGas(cstate, igas, mmr_gas(:,igas), rc)
        if (rc /=0) stop "    *** CARMASTATE_SetGas FAILED ***"
      end do

      # Execute the step
      call CARMASTATE_Step(cstate, rc)
      if (rc /=0) stop "    *** CARMASTATE_Step FAILED ***"
        
      # Get the updated bin mmr.
      do ielem = 1, NELEM
        do ibin = 1, NBIN
          call CARMASTATE_GetBin(cstate, ielem, ibin, mmr(:,ielem,ibin), rc)
          if (rc /=0) stop "    *** CARMASTATE_GetBin FAILED ***"
        end do
      end do

      # Get the updated gas mmr.
      do igas = 1, NGAS
        call CARMASTATE_GetGas(cstate, igas, &
                              mmr_gas(:,igas), rc, &
                              satliq=satliq(:,igas), &
                              satice=satice(:,igas))
        if (rc /=0) stop "    *** CARMASTATE_GetGas FAILED ***"
      end do

      # Get the updated temperature.
      call CARMASTATE_GetState(cstate, rc, t=t(:), rlheat=rlheat(:))
      if (rc /=0) stop "    *** CARMASTATE_Get FAILED ***"

      # Write output for the falltest
      write(lun,'(f12.0)') istep*dtime
      
      write(lun,'(2g16.5)') t(1) - t_orig, rlheat(1)

      do ielem = 1, NELEM
        do ibin = 1, NBIN
          write(lun,'(2i4,e12.3)') ielem, ibin, real(mmr(1,ielem,ibin))
        end do
      end do
    
      do igas = 1, NGAS
        write(lun,'(i4,3e12.3)') igas, real(mmr_gas(1,igas)), satliq(1,igas), satice(1,igas)
      end do
    end do   # time loop
    
    println("Done")
end
  
  