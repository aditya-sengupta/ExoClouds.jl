using Distributions: UnivariateDistribution
using QuadGK: quadgk
using Interpolations: LinearInterpolation, extrapolate

# surprised Julia doesn't have a canonical implementation for this
moment(d::UnivariateDistribution, n::Float64, lbd::Float64=0, ubd::Float64=Inf) = quadgk(x -> x^n * pdf(d, x), lbd, ubd)[1]

"""
Linear interpolation that matches the behaviour of np.interp,
i.e. instead of throwing a bounds error, just return f(lower bound) or f(upper bound)
currently this remakes the interpolation N times, not ideal
"""
function fixed_interp(x, xp, fp) 
    if x < minimum(xp)
        return fp[argmin(xp)]
    elseif x > maximum(xp)
        return fp[argmax(xp)]
    else
        return LinearInterpolation(xp, fp)(x)
    end
end
