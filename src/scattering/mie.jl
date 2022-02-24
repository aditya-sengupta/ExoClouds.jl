# so there's already a Julia package that kinda does this, but it doesn't quite meet my use case
# therefore, this'll be a port of the Python PyMieScatt package
# with some other stuff maybe idk
# this may be split off into its own package at a later date, tbd

using SpecialFunctions
using Unitful: Length
using LinearAlgebra: ⋅

function mie_efficiencies(
    refractive_index::Complex,
    wavelength::Length,
    diameter::Length;
    n_medium::Real = 1.0
)
#  http://pymiescatt.readthedocs.io/en/latest/forward.html#MieQ
    refractive_index /= n_medium
    wavelength /= n_medium
    x = pi * diameter / wavelength
    if x == 0
        return 0, 0, 0, 1.5, 0, 0, 0
    elseif x <= 0.05
        # rayleigh scattering
        LL = (refractive_index^2-1)/(refractive_index^2+2) # Lorentz-Lorenz term
        qsca = 8 * abs2(LL) * (x^4)/3 # B&H eq 5.8
        qabs = 4 * x * imag(LL) # B&H eq. 5.11
        qext = qsca + qabs
        qback = 1.5 * qsca # B&H eq. 5.9
        qratio = 1.5
        g = 0
        qpr = qext
        return qext, qsca, qabs, g, qpr, qback, qratio
    else
        # the Mie expansion
        nmax = round(2 + x + 4 * (x^(1/3))) |> Int64
        n = 1:nmax
        n1 = 2 .* n .+ 1
        n2 = n .* (n .+ 2) ./ (n .+ 1)
        n3 = n1 ./ (n .* (n .+ 1))
        x2 = x^2

        an, bn = mie_external_coefficients(refractive_index, x) #Mie_ab in the Python

        q_extinction = (2/x2) * (n1 ⋅ (real.(an) .+ real.(bn)))
        q_scattering = (2/x2) * (n1 ⋅ (abs2.(an) .+ abs2.(bn)))
        qabs = q_extinction - q_scattering

        coeff_components = [real(an), imag(an), real(bn), imag(bn)]
        g0 = map(x -> view(x, 1:nmax-1), coeff_components)
        g1 = map(x -> view(x, 2:nmax), coeff_components)

        g = (4/(q_scattering * x2)) * (
            n2[1:end-1] ⋅ sum(g0[i] .* g1[i] for i in 1:4) + n3 ⋅ ((real(an) .* real(bn) .+ imag(an) .* imag(bn)))
        )
        qpr = q_extinction - q_scattering * g
        qback = (1/x2) * abs2((n1 ⋅ (((-1) .^ n) .* (an .- bn))))
        
        qratio = qback / q_scattering
        
        return q_extinction, q_scattering, qabs, g, qpr, qback, qratio
    end # if
end # function

function mie_external_coefficients(
    m::Complex, # refractive index
    x::Real # size parameter
)
#  http://pymiescatt.readthedocs.io/en/latest/forward.html#Mie_ab
    mx = m * x
    nmax = 2 + x + 4*(x^(1/3)) |> round |> Int
    nmx = max(nmax, abs(mx))+16 |> round |> Int
    n = 1:nmax
    ν = n .+ 0.5

    sx = sqrt(π * x / 2)

    px = sx * besselj.(ν, x)
    p1x = cat(sin(x), px[1:nmax-1], dims = 1)

    chx = -sx * bessely.(ν, x)
    ch1x = cat(cos(x), chx[1:nmax-1], dims = 1)

    gsx = px - 1im * chx
    gs1x = p1x - 1im * ch1x

    # B&H Equation 4.89
    Dn = zeros(ComplexF64, nmx)
    for i in nmx:-1:2
        Dn[i-1] = (i/mx)-(1/(Dn[i]+i/mx))
    end

    D = Dn[1:nmax] # Dn(mx), drop terms beyond nMax
    da = D ./ m .+ n ./ x
    db = m .* D .+ n ./ x

    an = (da .* px .- p1x) ./ (da .* gsx .- gs1x)
    bn = (db .* px .- p1x) ./ (db .* gsx .- gs1x)

    return an, bn
end