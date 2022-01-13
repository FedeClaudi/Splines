module Interpolation
    import StaticArrays: @MMatrix


    include("Geometry.jl")
    import .Geometry: Point, Points

    export lerp, PiecewiseLinear, PiecewiseLinear, BSpline, Bezier, Bezier!, RationalBezier, RationalBezier!

    const Knots = Vector{Float64}

    # ---------------------------------------------------------------------------- #
    #                                     CURVE                                    #
    # ---------------------------------------------------------------------------- #
    """
        Curve

    Stores the results of fitting a spline to a set of nodes.
    Contains the parameter values (`τ`) at which the fitted spline
    was evaluated at, the coordinates (`points`) of the spline and a 
    function (`f`) with signature `f(t::Number)::Vector{Number}` giving
    the position of the spline at a parameter value `t ∈ [0,1]`.
    """
    struct Curve
        name::String
        τ::Vector{Float64}   # parameter values used when constructing curve
        points::Points       # points along curve at param values
        func                 # function to evaluate the curve at a param value
    end

    """
        Make `Curve` callable such that `mycurve = Curve(...)` and then `mycurv(t)`
        results in evaluating the curve function at the parameter value `t`.
    """
    function (curve::Curve)(t)
        curve.func(t)
    end

    """
        Custom printing of curve object
    """
    function Base.show(io::IO, curve::Curve)
        print(io, "Curve: '$(curve.name)' $(size(curve.points))")
    end


    # ----------------------------- helper functions ----------------------------- #
    """
        Get the order of a curve (n-1) given a set of n nodes
    """
    ν(X::Points) = size(X, 2)-1


    """
        Creates an empty MMatrix with the same dimensions as the curve generate 
        by a spline function.
    """
    function init_curve(nodes::Points, δt::Float64)::Points
        ndim = size(nodes, 1)  # number of dimensions
        nt = Int(1/δt)  # number of parameter steps
        return @MMatrix zeros(ndim, nt)
    end

    """
        Alteranate method signature, for when the dimensions of the target curve
        are alraedy known.
    """
    init_curve(ndim::Int, nt::Int) = @MMatrix zeros(ndim, nt)


    """
        Common computations carried out before fitting any spline: getting a range
        over the parameter's interval, 'closing' the curve by repeating a node etc..
    """
    function prep_spline_parameters(nodes::Points, δt::Float64, closed::Bool)

        if closed
            nodes = [nodes nodes[:, 1]] # repeat first control point to make it loop
        end

        ndim = size(nodes, 1)
        n = ν(nodes)
        τ = Array(0:δt:1-δt)  # parameter range

        return nodes, ndim, n, τ
    end
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
    function PiecewiseLinear(nodes::Points; η::Int64=100, closed::Bool=false)::Curve
        # compute quantities
        N = size(nodes, 2)
        n_nodes = closed ? N : N - 1
        ηₜ = η * n_nodes
        return PiecewiseLinear!(init_curve(size(nodes, 1), ηₜ), nodes; η=η, closed=closed)
    end

    """
        PiecewiseLinear!(curve_points, nodes[; η=100])

    In-place piecewise linear interpolation between consecutive pairs of nodes called nodes, given a
    pre-allocated curve_points array (d x η*N Points array).
    """
    function PiecewiseLinear!(curve_points::Points, nodes::Points; η::Int64=100, closed::Bool=false)::Curve
        # compute quantities
        N = size(nodes, 2)
        n_nodes = closed ? N : N - 1
        P = range(0, 1, length=η)  # for interpolation
        
        # create curve_points by linear interpolation of each segment
        for n in range(1, length=n_nodes)
            if n == N
                k₀, k₁ = nodes[:, end], nodes[:, 1]
            else
                k₀, k₁ = nodes[:, n], nodes[:, n+1]
            end
            curve_points[:, η * (n - 1) + 1: η * n] = lerp(k₀, k₁, P)'
        end
    
        return Curve(
            "PiecewiseLinear",
            0:1,  # no params yet in piecewise linear
            curve_points,
            (t) -> print("No evaluating function for piecewise linear")
        )
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


    # ----------------------------- B-SPLINE FUNCTION ---------------------------- #
    """
        BSpline(X; d, [δt=0.01, knots_type=:uniform])

    Computes the b-spline of degree 'd' given a set of control nodes (d x N, of type ::Points).
    The paramtert 'δt' specifies how densly to sample the parameter interval τ=[0,1] (i.e. how many
    points in the bspline curve).

    Value of a d-dimensional b-spline defined by a set of nkots `k` at `t ∈ τ = [0,1]` given
    control points `X` (dxN array of points).

        `S(t) = ∑ᵢⁿ Nᵢ(t)Xᵢ`
    """
    function BSpline(nodes::Points; d::Int, δt::Float64=.01, knots_type::Symbol=:uniform, closed::Bool=false)::Curve
        BSpline!(init_curve(nodes, δt), nodes; d=d, δt=δt, knots_type=knots_type, closed=closed)
    end

    """
        BSpline(t; nodes, d, knots)

    Evaluate a BSpline function at a parameter value t
    """
    function BSpline(t::Number; nodes::Points, d::Int, knots::Knots)
        sum(N(τ; k=knots, i=i, d=d)' .* nodes[:, i+1], 0:ν(nodes))
    end

    """
        In place implementation of BSpline function. See `BSpline`
    """
    function BSpline!(curve_points::Points, nodes::Points; d::Int, δt::Float64=.01, knots_type::Symbol=:uniform, closed::Bool=false)::Curve
        if d == 1
            @warn "For b-splines with d=1, `PiecewiseLinear` offers a more efficient implementation"
        end
        
        # prep
        nodes, ndim, n, τ = prep_spline_parameters(nodes, δt, closed)
        k = eval(:($knots_type($n, $d)))  # initialize knots with the desired method

        # compute spline
        B(i) = N(τ; k=k, i=i, d=d)' .* nodes[:, i+1]
        curve_points = sum(B, 0:n)

        return Curve(
            "Bspline (d=$d)",
            τ,
            curve_points,
            t -> BSpline(t; nodes=nodes, d=d, knots=k)
        )
    end

    # ---------------------------------------------------------------------------- #
    #                                    Bezier                                    #
    # ---------------------------------------------------------------------------- #

    """
        Bernstein(t::Float64; i::Int, n::Int)

    Evaluate the Bernstein polynomial at parameter value `t` given the index `i` and the number
    of polynomials `n`.
    """
    Bernstein(τ::Vector{Float64}; i::Int, n::Int)::Vector{Float64} = @. binomial(n, i) * τ^i * (1 - τ)^(n-i) 
    Bernstein(t::Number; i::Int, n::Int)::Float64 = binomial(n, i) * t^i * (1 - t)^(n-i) 

    """
        Bezier(nodes::Points;   δt::Float64=.01, knots_type::Symbol=:uniform, closed::Bool=false)

    Compute the Bezier curve given a set of `nodes` (d x N array of points)
    """
    function Bezier(nodes::Points; δt::Float64=.01, closed::Bool=false)::Curve
        return Bezier!(init_curve(nodes, δt), nodes; δt=δt, closed=closed)
    end

    """
        Bezier(t, X)

    Evaluate a Bezier function at a given parameter value `t` and given a set
    of control nodes `X`
    """
    Bezier(t::Number; nodes::Points) = sum(
        (i)->Bernstein(t; i=i, n=ν(nodes))' .* nodes[:, i+1], 0:ν(nodes)
    )

    """
        In place Bezier curve computation, see `Bezier` for more details.
    """
    function Bezier!(curve_points::Points, nodes::Points; δt::Float64=.01, closed::Bool=false)::Curve
        # prepare parameters
        nodes, ndim, n, τ = prep_spline_parameters(nodes, δt, closed)

        # compute bezier curve
        B(i) = Bernstein(τ; i=i, n=n)' .* nodes[:, i+1]
        curve_points = sum(B, 0:n)

        return Curve(
            "Bezier",
            τ,
            curve_points,
            (t) -> Bezier(t; nodes=nodes)
        )
    end



    # ---------------------------------------------------------------------------- #
    #                                Rational Bezier                               #
    # ---------------------------------------------------------------------------- #

    """
        RationalBezier(nodes::Points, weights::Vector{Float64}; δt::Float64=.01, closed::Bool=false)::Points
    
    Computes a rational Bezier curve: https://en.wikipedia.org/wiki/B%C3%A9zier_curve#Rational_B%C3%A9zier_curves
    Similar to `Bezier` but with an additional 1xN weights vector.
    """
    function RationalBezier(nodes::Points, weights::Vector{Float64}; δt::Float64=.01, closed::Bool=false)::Curve
        return RationalBezier!(init_curve(nodes, δt), nodes, weights; δt=δt, closed=closed)
    end

    """
        RationalBezier(t, X)

    Evaluate a RationalBezier function at a given parameter value `t` and given a set
    of control nodes `X`
    """
    function RationalBezier(t::Number; nodes::Points, weights::Vector{Float64})
        n = ν(nodes) # number of control points

        numerator(i) = Bernstein(t; i=i, n=n)' .* nodes[:, i+1] .* weights[i+1]
        denominator(i) = Bernstein(t; i=i, n=n)' .* weights[i+1]
        return sum(numerator, 0:n) ./ sum(denominator, 0:n)
    end

    """
        In-place rational Bezier curve computation, see `RationalBezier`
    """
    function RationalBezier!(curve_points::Points, nodes::Points, weights::Vector{Float64}; δt::Float64=.01, closed::Bool=false)::Curve
        # prepare parameters
        nodes, ndim, n, τ = prep_spline_parameters(nodes, δt, closed)

        numerator(i) = Bernstein(τ; i=i, n=n)' .* nodes[:, i+1] .* weights[i+1]
        denominator(i) = Bernstein(τ; i=i, n=n)' .* weights[i+1]
        curve_points = sum(numerator, 0:n) ./ sum(denominator, 0:n)

        return Curve(
            "RationalBezier",
            τ,
            curve_points,
            (t) -> RationalBezier(t; nodes=nodes, weights=weights)
        )
    end


end

