using Distributions: UnivariateDistribution
using QuadGK: quadgk

# surprised Julia doesn't have a canonical implementation for this
moment(d::UnivariateDistribution, n::Real, lbd::Real=0, ubd::Real=Inf) = quadgk(x -> x^n * pdf(d, x), lbd, ubd)[1]