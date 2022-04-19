# virga/optics.py

using DelimitedFiles

function calc_optics(gases::Vector{Element}, nwave, qc, rg, ndz, radii, dr, qext, qscat, cos_qscat, sig, rmin)
    nz = size(qc, 1)
    ngas = size(qc, 2)
    nrad = length(radii)

    opd_layer = zeros((nz, ngas))
    scat_gas = zeros((nz,nwave,ngas))
    ext_gas = zeros((nz,nwave,ngas))
    cqs_gas = zeros((nz,nwave,ngas))
    opd = zeros((nz,nwave))
    opd_gas = zeros((nz,ngas))
    w0 = zeros((nz,nwave))
    g0 = zeros((nz,nwave))
    for iz in 1:nz
        for (igas, gas) in enumerate(gases)
            # Optical depth for conservative geometric scatterers 
            if ndz[iz,igas] > 0
                if log10(rg[iz,igas] / rmin) < 0.75*sig
                     @warn "Take caution in analyzing results. There have been a calculated particle radii off the Mie grid, which has a min radius of $(rmin) and distribution of $(sig). Error at$(rg[iz,igas]) for $(gas) at the $(iz)th grid point."
                end

                r2 = rg[iz,igas]^2 * exp(2*log(sig)^2 )
                opd_layer[iz,igas] = 2π * r2 * ndz[iz,igas]

                #  Calculate normalization factor (forces lognormal sum = 1.0)
                rsig = sig
                norm = 0.
                for irad in 1:nrad
                    rr = radii[irad]
                    arg1 = dr[irad] / ( sqrt(2π)*rr*log(rsig))
                    arg2 = -log(rr/rg[iz,igas] )^2 / ( 2*log(rsig)^2)
                    norm = norm + arg1*exp( arg2 )
                    #print (rr, rg[iz,igas],rsig,arg1,arg2)
                end

                # normalization
                norm = ndz[iz,igas] / norm

                for irad in 1:nrad
                    rr = radius[irad]
                    arg1 = dr[irad] / ( sqrt(2π)*log(rsig) )
                    arg2 = -log( rr/rg[iz,igas] )^2 / ( 2*log(rsig)^2 )
                    pir2ndz = norm * π * rr * arg1 * exp(arg2)         
                    for iwave in 1:nwave
                        scat_gas[iz,iwave,igas] = scat_gas[iz,iwave,igas]+qscat[iwave,irad,igas]*pir2ndz
                        ext_gas[iz,iwave,igas] = ext_gas[iz,iwave,igas]+qext[iwave,irad,igas]*pir2ndz
                        cqs_gas[iz,iwave,igas] = cqs_gas[iz,iwave,igas]+cos_qscat[iwave,irad,igas]*pir2ndz
                    end
                end
            end
        end
    end

    #TO DO ADD IN CLOUD SUBLAYER KLUGE LATER 

    #Sum over gases and compute spectral optical depth profile etc
    for iz in 1:nz
        for iwave in 1:nwave
            opd_scat = 0.
            opd_ext = 0.
            cos_qs = 0.
            for igas in 1:ngas
                opd_scat = opd_scat + scat_gas[iz,iwave,igas]
                opd_ext = opd_ext + ext_gas[iz,iwave,igas]
                cos_qs = cos_qs + cqs_gas[iz,iwave,igas]

                if opd_scat > 0
                    opd[iz,iwave] = opd_ext
                    w0[iz,iwave] = opd_scat / opd_ext
                    g0[iz,iwave] = cos_qs / opd_scat
                end
            end
        end
    end
                    
    #cumulative optical depths for conservative geometric scatterers
    for igas in 1:ngas
        opd_gas[:,igas] = cumsum(opd_layer[:,igas])
    end
    return opd, w0, g0, opd_gas
end

function refr_ind(e::Element; path=joinpath(PROJECT_ROOT, data, refrinds)) 
    DataFrame(
        readdlm(
            joinpath(path, unsubscribe(typeof(e)) * ".refrind")
        )[:, 2:4],
        ["λ", "nn", "kk"]
    )
end