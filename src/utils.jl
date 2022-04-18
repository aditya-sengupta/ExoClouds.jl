# you may wonder why i bother having the subscripts if I'm just going to raise them, and to that I simply answer: The Aesthetic
Interp = Gridded(Linear())
Extrap = Flat()

import Interpolations.extrapolate

"""
Raises the subscripts in element names.
"""
function unsubscribe(element::String)::String
    raiser = Dict("₂" => "2", "₃" => "3", "₄" => "4")
    map(x -> get(raiser, x, x), split(element, "")) |> join
end

"""
Standardises on a package-wide extrapolation scheme
"""
function Interpolations.extrapolate(x::AbstractArray, y::AbstractArray)
    itp = interpolate((x,), y, Interp)
    return extrapolate(itp, Extrap)
end

# indexing into these objects is by the z value not by index, 
# so this is for if you want to loop over the actual values
each(ex::Extrapolation) = ex.itp.coefs 