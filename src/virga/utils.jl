using QuadGK: quadgk

# surprised Julia doesn't have a canonical implementation for this
moment(f, n::Float64, lbd::Float64=0.0, ubd::Float64=Inf) = quadgk(x -> x^n * f(x), lbd, ubd)[1]

"""
Warning
-------
Original code from A&M code. 
Discontinued function. See 'get_r_grid'.

Get spacing of radii to run Mie code

r_min : float 
    Minimum radius to compute (cm)

n_radii : int
    Number of radii to compute 
    """
function get_r_grid(r_min=1e-8cm, n_radii=60; vrat=2.2, pw=1/3)
    f1 = ( 2.0*vrat / ( 1.0 + vrat) )^pw
    f2 = (( 2.0 / ( 1.0 + vrat ) )^pw) * (vrat^(pw-1.0))

    radius = r_min * vrat .^ ((0:(1-1/n_radii):n_radii-1) ./ 3)
    rup = f1*radius
    dr = f2*radius

    return radius, rup, dr
end