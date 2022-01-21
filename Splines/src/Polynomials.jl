"""
    Polynomials functions for spline interpolation: basis functions for b-splines
    and bernstein polynomials for bezier curves.

    Most polynomials have two implementations: one in which the polynomial is 
    computed over the entier parameter range, the other for evaluating the 
    polynomial at a single parameter value.
"""
module Polynomials

    # ---------------------------------------------------------------------------- #
    #                           B-Spline BASIS FUNCTIONS                           #
    # ---------------------------------------------------------------------------- #
    """
        N_0(τ, k; i=0)

    Zero-th order basis function for b-splines for the i-th knot.
        for t∈[k[i], k[i+1]] = 1
        else: 0

    The additional "+1" index in the code is because in the bspline maths 
    the knots are indexed starting from 0.
    """
    function N_0(τ; k::AbstractArray, i::Int=0)
        if i+2 < length(k)
            return k[i+1] .<= τ .< k[i+2]
        else
            return k[i+1] .<= τ .<= k[i+2]
        end
    end

    """
        ω(τ, j)

    Function used within N_D to compute the factors for the linear interpolation
    of the two lover level basis functions
    """
    function ω(τ, j, d, k)
        Δ = k[j+d]-k[j]
        @. Δ!=0 ? (τ - k[j])/Δ : zero(τ)
    end

    """
        N_D(τ, k; i=0, d=0)

    D-th order basis function for b-splines for the i-th knot. 
    Built recursively on lower order basis functions.
    See: https://www.geometrictools.com/Documentation/BSplineCurveLeastSquaresFit.pdf
    and: https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-basis.html

    The "j=i+1" index in the code is because in the bspline maths 
    the knots are indexed starting from 0.
    """
    function N_D(τ; k::AbstractArray, i::Int=0, d::Int=0)
        ω(τ, i+1, d, k) .* N_basis(τ; k=k, i=i, d=d-1) .+ (1 .- ω(τ, i+1+1, d, k)) .* N_basis(τ; k=k, i=i+1, d=d-1)
    end

    """
        N(τ, k; i=0, d=0)

    B-spline basis function for the i-th knot and d-th order.
    Calls either N_0 or N_D depending on the value of d.
    """
    N_basis(τ; k::AbstractArray, i::Int=0, d::Int=0) = (d == 0) ? N_0(τ; k, i=i) : N_D(τ; k, i=i, d=d)


    # ---------------------------------------------------------------------------- #
    #                             BERNSTEIN POLYNOMIALS                            #
    # ---------------------------------------------------------------------------- #
    """
        bernstein(t::Float64; i::Int, n::Int)

    Evaluate the bernstein polynomial at parameter value `t` given the index `i` and the number
    of polynomials `n`.
    """
    bernstein(τ; i::Int, n::Int) = @. binomial(n, i) * τ^i * (1 - τ)^(n-i) 


end