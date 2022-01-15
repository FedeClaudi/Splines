"""
    Polynomials functions for spline interpolation: basis functions for b-splines
    and bernstein polynomials for Bezier curves.

    Most polynomials have two implementations: one in which the polynomial is 
    computed over the entier parameter range, the other for evaluating the 
    polynomial at a single parameter value.
"""
module Polynomials
    include("Types.jl")

    import .Types: Point, Points, Knots

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
    N_0(τ::Vector{Float64}; k::Knots, i::Int=0)::Vector{Float64} = (k[i+1] .<= τ .< k[i+2])
    N_0(t::Number; k::Knots, i::Int=0)::Number = (k[i+1] <= t < k[i+2]) ? 1 : 0

    """
        N_D(τ, k; i=0, d=0)

    D-th order basis function for b-splines for the i-th knot. 
    Built recursively on lower order basis functions.
    See: https://www.geometrictools.com/Documentation/BSplineCurveLeastSquaresFit.pdf
    and: https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-basis.html

    The "j=i+1" index in the code is because in the bspline maths 
    the knots are indexed starting from 0.
    """
    function N_D(τ::Vector{Float64}; k::Knots, i::Int=0, d::Int=0)::Vector{Float64}
        j = i + 1

        kⱼ = k[j]
        k̂ = k[j + d + 1]

        α = (τ .- kⱼ)/(k[j+d] - kⱼ + eps()) 
        β = (k̂ .- τ)/(k̂ - k[j+1] + eps())
        
        α .* N(τ; k=k, i=i, d=d-1) .+ β .* N(τ; k=k, i=i+1, d=d-1)
    end

    function N_D(t::Number; k::Knots, i::Int=0, d::Int=0)::Number
        j = i + 1

        kⱼ = k[j]
        k̂ = k[j + d + 1]

        α = (t - kⱼ)/(k[j+d] - kⱼ + eps()) 
        β = (k̂ - t)/(k̂ - k[j+1] + eps())
        
        α * N(t; k=k, i=i, d=d-1) + β .* N(t; k=k, i=i+1, d=d-1)
    end

    """
        N(τ, k; i=0, d=0)

    B-spline basis function for the i-th knot and d-th order.
    Calls either N_0 or N_D depending on the value of d.
    """
    N(τ::Vector{Float64}; k::Knots, i::Int=0, d::Int=0)::Vector{Float64} = (d == 0) ? N_0(τ; k, i=i) : N_D(τ; k, i=i, d=d)
    N(t::Number; k::Knots, i::Int=0, d::Int=0)::Number = (d == 0) ? N_0(t; k, i=i) : N_D(t; k, i=i, d=d)



    # ---------------------------------------------------------------------------- #
    #                             BERNSTEIN POLYNOMIALS                            #
    # ---------------------------------------------------------------------------- #
    """
        bernstein(t::Float64; i::Int, n::Int)

    Evaluate the bernstein polynomial at parameter value `t` given the index `i` and the number
    of polynomials `n`.
    """
    bernstein(t::Number; i::Int, n::Int)::Float64 = binomial(n, i) * t^i * (1 - t)^(n-i) 

    """
        bernstein(τ::Vector{Float64}; i::Int, n::Int)::Vector{Float64}

    Evaluate the polynomial over an entire parametr range: `τ = 0:δt:1`
    """
    bernstein(τ::AbstractArray; i::Int, n::Int)::Vector{Float64} = @. binomial(n, i) * τ^i * (1 - τ)^(n-i) 


    """
        bernstein(τ::AbstractArray; i::AbstractArray, n::Int)
    
    Evaluate the polynomial over an entire parameter range and for all degrees `i ∈ [0, n]``
    """
    function Bernstein2(τ::AbstractArray; i::AbstractArray, n::Int)::AbstractArray
        binom = [binomial(n, j) for j in i]'
        @. binom * τ^(i') * (1 - τ)^(n-i)'
    end

end