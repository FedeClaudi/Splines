module Curves

    include("Utils.jl")
    include("Types.jl")
    include("Polynomials.jl")

    using .Types
    import .Utils: ν, init_curve, prep_spline_parameters, uniform, periodic
    import .Polynomials: N, bernstein

    export BSpline, Bezier, Bezier!, RationalBezier, RationalBezier!

    # ---------------------------------------------------------------------------- #
    #                                   B-SPLINE                                   #
    # ---------------------------------------------------------------------------- #
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
        # prep
        nodes, ndim, n, τ = prep_spline_parameters(nodes, δt, closed)
        k = eval(:($knots_type($n, $d)))  # initialize knots with the desired method

        # compute spline
        B(i) = N(τ; k=k, i=i, d=d)' .* nodes[:, i+1]
        curve_points += sum(B, 0:n)

        return Curve(
            "Bspline (d=$d)",
            nodes,
            τ,
            curve_points,
            t -> BSpline(t; nodes=nodes, d=d, knots=k)
        )
    end

    # ---------------------------------------------------------------------------- #
    #                                    Bezier                                    #
    # ---------------------------------------------------------------------------- #


    """
        Bezier(nodes::Points;   δt::Float64=.01, knots_type::Symbol=:uniform, closed::Bool=false)

    Compute the Bezier curve given a set of `nodes` (d x N array of points)
    """
    Bezier(nodes::Points; δt::Float64=.01, closed::Bool=false)::Curve = Bezier!(
                                            init_curve(nodes, δt), nodes; δt=δt, closed=closed
                                        )

    """
        Bezier(t, X)

    Evaluate a Bezier function at a given parameter value `t` and given a set
    of control nodes `X`
    """
    Bezier(t::Number; nodes::Points) = sum(
        (i)->bernstein(t; i=i, n=ν(nodes))' .* nodes[:, i+1], 0:ν(nodes)
    )

    """
        In place Bezier curve computation, see `Bezier` for more details.
    """
    function Bezier!(curve_points::Points, nodes::Points; δt::Float64=.01, closed::Bool=false)::Curve
        # prepare parameters
        nodes, ndim, n, τ = prep_spline_parameters(nodes, δt, closed)

        # compute bezier curve
        B(i) = bernstein(τ; i=i, n=n)' .* nodes[:, i+1]
        curve_points = sum(B, 0:n)

        return Curve(
            "Bezier",
            nodes,
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

        numerator(i) = bernstein(t; i=i, n=n)' .* nodes[:, i+1] .* weights[i+1]
        denominator(i) = bernstein(t; i=i, n=n)' .* weights[i+1]
        return sum(numerator, 0:n) ./ sum(denominator, 0:n)
    end

    """
        In-place rational Bezier curve computation, see `RationalBezier`
    """
    function RationalBezier!(curve_points::Points, nodes::Points, weights::Vector{Float64}; δt::Float64=.01, closed::Bool=false)::Curve
        # prepare parameters
        nodes, ndim, n, τ = prep_spline_parameters(nodes, δt, closed)

        numerator(i) = bernstein(τ; i=i, n=n)' .* nodes[:, i+1] .* weights[i+1]
        denominator(i) = bernstein(τ; i=i, n=n)' .* weights[i+1]
        curve_points = sum(numerator, 0:n) ./ sum(denominator, 0:n)

        return Curve(
            "RationalBezier",
            nodes,
            τ,
            curve_points,
            (t) -> RationalBezier(t; nodes=nodes, weights=weights)
        )
    end


end

