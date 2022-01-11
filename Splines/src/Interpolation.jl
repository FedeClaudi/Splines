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
        N_0(t, k; i=0)

    Zero-th order basis function for b-splines for the i-th knot.
        for t∈[k[i], k[i+1]] = 1
        else: 0

    The additional "+1" index in the code is because in the bspline maths 
    the knots are indexed starting from 0.
    """
    N_0(t::Number; k::Knots, i::Int=0) = (k[i+1] <= t < k[i+2]) ? 1 : 0

    """
        N_D(t, k; i=0, d=0)
    
    D-th order basis function for b-splines for the i-th knot. 
    Built recursively on lower order basis functions.
    See: https://www.geometrictools.com/Documentation/BSplineCurveLeastSquaresFit.pdf
    and: https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-basis.html

    The "j=i+1" index in the code is because in the bspline maths 
    the knots are indexed starting from 0.
    """
    function N_D(t::Float64; k::Knots, i::Int=0, d::Int=0)::Number
        j = i + 1

        kⱼ = k[j]
        k̂ = k[j + d + 1]
    
        α = (t - kⱼ)/(k[j+d] - kⱼ + eps()) 
        β = (k̂ - t)/(k̂ - k[j+1] + eps())
        
        α * N(t; k=k, i=i, d=d-1) + β * N(t; k=k, i=i+1, d=d-1)
    end


    """
        N(t, k; i=0, d=0)

    B-spline basis function for the i-th knot and d-th order.
    Calls either N_0 or N_D depending on the value of d.
    """
    N(t::Number; k::Knots, i::Int=0, d::Int=0)::Number = (d == 0) ? N_0(t; k, i=i) : N_D(t; k, i=i, d=d)


    # ----------------------------- B-SPLINE FUNCTION ---------------------------- #

    """
        bspline(t::Number, k::Knots, X::Points; d::Int=1) 
    
    Value of a d-dimensional b-spline defined by a set of nkots `k` at `t ∈ [0,1]` given
    control points `X` (dxN array of points).

        `S(t) = ∑ᵢⁿ Nᵢ(t)Xᵢ`
    """
    function eval_bspline(t::Number; k::Knots, nodes::Points, d::Int=1)
        # Sum the product of the n-many basis functions with the control points.
        out = zeros(size(nodes, 1))
        for i in 0:size(nodes, 2)-1  # n-many basis functions
            out += N(t; k, i=i, d=d) .* nodes[:, i+1]
        end
        return out
    end

    """
        eval_bspline(τ::Vector{Float64}, k::Knots, X::Points; d::Int=1)::Points

    Method for evalutaing the bspline function over a vector of values
    """
    function eval_bspline(τ::Vector{Float64}; k::Knots, nodes::Points, d::Int=1)
        hcat(eval_bspline.(τ; k, nodes, d=d)...)
    end


    """
        BSpline(X; d, [δt=0.01, knots_type=:uniform])

    Computes the b-spline of degree 'd' given a set of control nodes (d x N, of type ::Points).
    The paramtert 'δt' specifies how densly to sample the paramter interval τ=[0,1] (i.e. how many
    points in the bspline curve).
    """
    function BSpline(nodes::Points; d::Int, δt::Float64=.01, knots_type::Symbol=:uniform, closed::Bool=false)::Points
        if d == 1
            @warn "For b-splines with d=1, `PiecewiseLinear` offers a more efficient implementation"
        end

        if closed
            nodes = [nodes nodes[:, 1]] # repeat first control point to make it loop
        end

        n = size(nodes, 2) - 1 # number of control points
        k = eval(:($knots_type($n, $d)))
        τ = Array(0:δt:1)  # domain

        return eval_bspline(τ; k, nodes, d=d)[:, 1:end-1]
    end

    """
        In place implementation of BSpline function.
    """
    function BSpline!(curve::Points, nodes::Points; d::Int, δt::Float64=.01, knots_type::Symbol=:uniform, closed::Bool=false)::Points
        curve = BSpline(nodes; d=d, δt=δt, knots_type=knots_type, closed=closed)
    end

    # ---------------------------------------------------------------------------- #
    #                                    Bezier                                    #
    # ---------------------------------------------------------------------------- #

    """
        Bernstein(t::Float64; i::Int, n::Int)

    Evaluate the Bernstein polynomial at paramter value `t` given the index `i` and the number
    of polynomials `n`.
    """
    Bernstein(t::Float64; i::Int, n::Int)::Float64 = binomial(n, i) * t^i * (1 - t)^(n-i) 
    Bernstein(τ::Vector{Float64}; i::Int, n::Int)::Vector{Float64} = @. binomial(n, i) * τ^i * (1 - τ)^(n-i) 

    """
        eval_bezier(t::Float64; i::Int, n::Int)
    
    Compute the value of a Bezier curve through `n` control `nodes`
    """
    eval_bezier(t::Float64; n::Int, nodes::Points) = sum((i)->Bernstein(t; i=i, n=n)* nodes[:, i+1], 0:n)
    eval_bezier(τ::AbstractArray; n::Int, nodes::Points) = sum((i)->Bernstein(τ; i=i, n=n)' .* nodes[:, i+1], 0:n)

    """
        Bezier(nodes::Points;   δt::Float64=.01, knots_type::Symbol=:uniform, closed::Bool=false)

    Compute the Bezier curve given a set of `nodes` (d x N array of points)
    """
    function Bezier(nodes::Points; δt::Float64=.01, closed::Bool=false)::Points
        if closed
            nodes = [nodes nodes[:, 1]] # repeat first control point to make it loop
        end

        n = size(nodes, 2) - 1 # number of control points
        τ = Array(0:δt:1)  # paramter interval 

        return eval_bezier(τ; n=n, nodes)[:, 1:end-1]
    end

    function Bezier!(curve::Points, nodes::Points; δt::Float64=.01, closed::Bool=false)::Points
        curve = Bezier(nodes; δt=δt, closed=closed)
    end
end

