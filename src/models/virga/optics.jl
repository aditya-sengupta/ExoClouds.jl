# virga/optics.py

using DelimitedFiles

function optics(
    condensibles::Vector{<:Molecule};
     nrad=40, rmin=1e-10, read_mie=true
)
    # skipping, as this currently has no return value
end

function refr_ind(m::Molecule; path=joinpath(PROJECT_ROOT, data, refrinds)) 
    DataFrame(
        readdlm(
            joinpath(path, unsubscribe(typeof(m)) * ".refrind")
        )[:, 2:4],
        ["λ", "nn", "kk"]
    )
end