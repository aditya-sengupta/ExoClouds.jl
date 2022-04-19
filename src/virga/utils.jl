using QuadGK: quadgk

# surprised Julia doesn't have a canonical implementation for this
moment(d::UnivariateDistribution, n::Float64, lbd::Float64=0.0, ubd::Float64=Inf) = quadgk(x -> x^n * pdf(d, x), lbd, ubd)[1]
