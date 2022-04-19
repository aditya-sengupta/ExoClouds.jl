using Interpolations: Extrapolation

# indexing into these objects is by the z value not by index, 
# so this is for if you want to loop over the actual values
each(ex::Extrapolation) = ex.itp.coefs 

findnearest(arr, val) = argmin(abs.(arr .- val))
