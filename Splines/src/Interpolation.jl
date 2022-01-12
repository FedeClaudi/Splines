module Interpolation
    import StaticArrays: @MMatrix


    include("Geometry.jl")
    import .Geometry: Point, Points

    export lerp, PiecewiseLinear, PiecewiseLinear, BSpline, Bezier, Bezier!

    const Knots = Vector{Float64}


    # ---------------------------------------------------------------------------- #
    #                             LINEAR interpolation                             #
    # ---------------------------------------------------------------------------- #

    """
        lerp(x₀, x₁, p)

    Linear interpolation between two values, given an interpolation factor
    `p ∈ [0, 1]`. `x₀` and `x₁` can be number of `Point`, but they need to have the same
    Type and shape (if `Point`).
    """

    lerp(x₀::Point, x₁::Point, p::Float64)::Point = x₀ * (1 - p) + x₁ * p
    lerp(x₀::Number, x₁::Number, p::Float64)::Number = x₀ * (1 - p) + x₁ * p
    lerp(x₀::Point, x₁::Point, p::AbstractArray)::AbstractArray = @. x₀' * (1 - p) + x₁' * p


    """
        PiecewiseLinear(nodes[; η=100])

    Piecewise linear interpolation between consecutive pairs of nodes 
    (called nodes, d x N array, `Points`) ∈ ℝⁿ to return
    a curve connecting them. The curve is specified by (η*N)-many points (d x η*N Points array), the 
    parameter η specifies the number of points in each linear segment.
    """
    function PiecewiseLinear(nodes::Points; η::Int64=100, closed::Bool=false)::Points
        # compute quantities
        N = size(nodes, 2)
        n_nodes = closed ? N : N - 1
        ηₜ = η * n_nodes
        curve = @MMatrix zeros(size(nodes, 1), ηₜ)
        return PiecewiseLinear!(curve, nodes; η=η, closed=closed)
    end

    """
        PiecewiseLinear!(curve, nodes[; η=100])

    In-place piecewise linear interpolation between consecutive pairs of nodes called nodes, given a
    pre-allocated curve array (d x η*N Points array).
    """
    function PiecewiseLinear!(curve::Points, nodes::Points; η::Int64=100, closed::Bool=false)::Points
        # compute quantities
        N = size(nodes, 2)
        n_nodes = closed ? N : N - 1
        P = range(0, 1, length=η)  # for interpolation
        
        # create curve by linear interpolation of each segment
        for n in range(1, length=n_nodes)
            if n == N
                k₀, k₁ = nodes[:, end], nodes[:, 1]
            else
                k₀, k₁ = nodes[:, n], nodes[:, n+1]
            end
            curve[:, η * (n - 1) + 1: η * n] = lerp(k₀, k₁, P)'
        end
    
        return curve
    end


    # ---------------------------------------------------------------------------- #
    #                            B-SPLINE interpolation                            #
    # ---------------------------------------------------------------------------- #

    # --------------------------- knots initialization --------------------------- #
    """
        uniform(n, d)

    Defines d+n+1 uniform knots ∈[0,1] for b-spline interpolation. 
    For i ∈ 0:n+d+1:
        if i ≤ d: k=0
        if d+1 ≤ i ≤ n: k=(i-d)/(n+1-d)
        else: k = 1

    See: https://www.geometrictools.com/Documentation/BSplineCurveLeastSquaresFit.pdf
    """
    function uniform(n::Int, d::Int)::Knots
        knots = zeros(n+d+1+1)  # second 1 is because start from 0
        for i in 0:n+d+1
            if i ≤ d
                knots[i+1] = 0
            elseif (d+1) ≤ i ≤ n
                knots[i+1] = (i-d)/(n+1-d)
            else
                knots[i+1] = 1
            end
        end
        return knots
    end


    """
        periodic(n, d)

    Defines d+n+1 periodic knots ∈[0,1] for b-spline interpolation. 

    See: https://www.geometrictools.com/Documentation/BSplineCurveLeastSquaresFit.pdf
    """
    periodic(n::Int, d::Int)::Knots = map((i)->(i-d)/(n+1-d), 0:n+d+1)

    # -------------------------- b-spline basis function ------------------------- #

    """
        N_0(τ, k; i=0)

    Zero-th order basis function for b-splines for the i-th knot.
        for t∈[k[i], k[i+1]] = 1
        else: 0

    The additional "+1" index in the code is because in the bspline maths 
    the knots are indexed starting from 0.
    """
    N_0(τ::Vector{Float64}; k::Knots, i::Int=0)::Vector{Float64} = (k[i+1] .<= τ .< k[i+2])

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


    """
        N(τ, k; i=0, d=0)

    B-spline basis function for the i-th knot and d-th order.
    Calls either N_0 or N_D depending on the value of d.
    """
    N(τ::Vector{Float64}; k::Knots, i::Int=0, d::Int=0)::Vector{Float64} = (d == 0) ? N_0(τ; k, i=i) : N_D(τ; k, i=i, d=d)


    # ----------------------------- B-SPLINE FUNCTION ---------------------------- #
    """
        BSpline(X; d, [δt=0.01, knots_type=:uniform])

    Computes the b-spline of degree 'd' given a set of control nodes (d x N, of type ::Points).
    The paramtert 'δt' specifies how densly to sample the paramter interval τ=[0,1] (i.e. how many
    points in the bspline curve).

    Value of a d-dimensional b-spline defined by a set of nkots `k` at `t ∈ τ = [0,1]` given
    control points `X` (dxN array of points).

        `S(t) = ∑ᵢⁿ Nᵢ(t)Xᵢ`
    """
    function BSpline(nodes::Points; d::Int, δt::Float64=.01, knots_type::Symbol=:uniform, closed::Bool=false)::Points
        ndim = size(nodes, 1)
        τ = Array(0:δt:1-δt)  # domain

        curve = @MMatrix zeros(ndim, length(τ)-1)

        return BSpline!(curve, nodes; d=d, δt=δt, knots_type=knots_type, closed=closed)
    end

    """
        In place implementation of BSpline function. See `BSpline`
    """
    function BSpline!(curve::Points, nodes::Points; d::Int, δt::Float64=.01, knots_type::Symbol=:uniform, closed::Bool=false)::Points
        if d == 1
            @warn "For b-splines with d=1, `PiecewiseLinear` offers a more efficient implementation"
        end

        if closed
            nodes = [nodes nodes[:, 1]] # repeat first control point to make it loop
        end

        ndim = size(nodes, 1)
        n = size(nodes, 2) - 1 # number of control points
        k = eval(:($knots_type($n, $d)))
        τ = Array(0:δt:1-δt)  # domain

        B(i) = N(τ; k=k, i=i, d=d)' .* nodes[:, i+1]
        curve = sum(B, 0:n)
    end

    # ---------------------------------------------------------------------------- #
    #                                    Bezier                                    #
    # ---------------------------------------------------------------------------- #

    """
        Bernstein(t::Float64; i::Int, n::Int)

    Evaluate the Bernstein polynomial at paramter value `t` given the index `i` and the number
    of polynomials `n`.
    """
    Bernstein(τ::Vector{Float64}; i::Int, n::Int)::Vector{Float64} = @. binomial(n, i) * τ^i * (1 - τ)^(n-i) 


    """
        Bezier(nodes::Points;   δt::Float64=.01, knots_type::Symbol=:uniform, closed::Bool=false)

    Compute the Bezier curve given a set of `nodes` (d x N array of points)
    """
    function Bezier(nodes::Points; δt::Float64=.01, closed::Bool=false)::Points
        ndim = size(nodes, 1)
        τ = Array(0:δt:1-δt)  # paramter interval 
        curve = @MMatrix zeros(ndim, length(τ)-1)

        return Bezier!(curve, nodes; δt=δt, closed=closed)
    end

    function Bezier!(curve::Points, nodes::Points; δt::Float64=.01, closed::Bool=false)::Points
        if closed
            nodes = [nodes nodes[:, 1]] # repeat first control point to make it loop
        end

        n = size(nodes, 2) - 1 # number of control points
        τ = Array(0:δt:1-δt)

        B(i) = Bernstein(τ; i=i, n=n)' .* nodes[:, i+1]
        curve = sum(B, 0:n)
    end


    # ---------------------------------------------------------------------------- #
    #                                Rational Bezier                               #
    # ---------------------------------------------------------------------------- #

end

