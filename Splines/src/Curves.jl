module Curves

    include("Utils.jl")
    include("Types.jl")
    include("Polynomials.jl")

    using .Types
    import .Utils: ν, init_curve, prep_spline_parameters, uniform, periodic
    import .Polynomials: N, bernstein

    export bspline, bezier, bezier!, rational_bezier, rational_bezier!, rational_bspline, rational_bspline!

    # ---------------------------------------------------------------------------- #
    #                                   B-SPLINE                                   #
    # ---------------------------------------------------------------------------- #
    """
        bspline(X; d, [δt=0.01, knots_type=:uniform])

    Computes the b-spline of degree 'd' given a set of control nodes (d x N, of type ::Points).
    The paramtert 'δt' specifies how densly to sample the parameter interval τ=[0,1] (i.e. how many
    points in the bspline curve).

    Value of a d-dimensional b-spline defined by a set of nkots `k` at `t ∈ τ = [0,1]` given
    control points `X` (dxN array of points).

        `S(t) = ∑ᵢⁿ Nᵢ(t)Xᵢ`
    """
    function bspline(nodes::Points; d::Int, δt::Float64=.01, knots_type::Symbol=:uniform, closed::Bool=false)::Curve
        bspline!(init_curve(nodes, δt), nodes; d=d, δt=δt, knots_type=knots_type, closed=closed)
    end

    """
        bspline(t; nodes, d, knots)

    Evaluate a bspline function at a parameter value t
    """
    function bspline(t::Number; nodes::Points, d::Int, knots::Knots)
        sum((i) -> N(t; k=knots, i=i, d=d)' .* nodes[:, i+1], 0:ν(nodes))
    end

    """
        In place implementation of bspline function. See `bspline`
    """
    function bspline!(curve_points::Points, nodes::Points; d::Int, δt::Float64=.01, knots_type::Symbol=:uniform, closed::Bool=false)::Curve
        # prep
        nodes, n, τ = prep_spline_parameters(nodes, δt, closed)
        k = eval(:($knots_type($n, $d)))  # initialize knots with the desired method

        # compute spline
        B(i) = N(τ; k=k, i=i, d=d)' .* nodes[:, i+1]
        curve_points += sum(B, 0:n)

        return Curve(
            "Bspline (d=$d)",
            nodes,
            τ,
            curve_points,
            t -> bspline(t; nodes=nodes, d=d, knots=k)
        )
    end

    

    # ---------------------------------------------------------------------------- #
    #                                  NURB curve                                  #
    # ---------------------------------------------------------------------------- #

    function rational_bspline(
                    nodes::Points, 
                    weights::Points; 
                    d::Int, 
                    δt::Float64=.01, 
                    knots_type::Symbol=:uniform, 
                    closed::Bool=false
                )::Curve
        rational_bspline!(init_curve(nodes, δt), nodes, weights; d=d, δt=δt, knots_type=knots_type, closed=closed)
    end

    function rational_bspline(t::Number; nodes::Points, weights::Points, d::Int, knots::Knots)
        num = sum((i) -> N(t; k=knots, i=i, d=d)' .* nodes[:, i+1] .* weights[i+1], 0:ν(nodes))
        den = sum((i) -> N(t; k=knots, i=i, d=d)' .* weights[i+1], 0:ν(nodes)) + eps()
        return num /den
    end

    function rational_bspline!(
            curve_points::Points,
            nodes::Points,
            weights::Points;
            d::Int,
            δt::Float64=.01,
            knots_type::Symbol=:uniform,
            closed::Bool=false
        )::Curve
        # prep
        nodes, n, τ = prep_spline_parameters(nodes, δt, closed)
        k = eval(:($knots_type($n, $d)))  # initialize knots with the desired method

        # compute spline
        Bₙ(i) = N(τ; k=k, i=i, d=d)' .* nodes[:, i+1] .* weights[i+1]
        Bₘ(i) =  N(τ; k=k, i=i, d=d)' .* weights[i+1]
        curve_points += sum(Bₙ, 0:n) ./ sum(Bₘ, 0:n)

        return Curve(
            "NURB curve (d=$d)",
            nodes,
            τ,
            curve_points,
            t -> rational_bspline(t; nodes=nodes, weights=weights, d=d, knots=k)
        )
    end


    # ---------------------------------------------------------------------------- #
    #                                    bezier                                    #
    # ---------------------------------------------------------------------------- #


    """
        bezier(nodes::Points;   δt::Float64=.01, knots_type::Symbol=:uniform, closed::Bool=false)

    Compute the bezier curve given a set of `nodes` (d x N array of points)
    """
    bezier(nodes::Points; δt::Float64=.01, closed::Bool=false)::Curve = bezier!(
                                            init_curve(nodes, δt), nodes; δt=δt, closed=closed
                                        )

    """
        bezier(t, X)

    Evaluate a bezier function at a given parameter value `t` and given a set
    of control nodes `X`
    """
    bezier(t::Number; nodes::Points) = sum(
        (i)->bernstein(t; i=i, n=ν(nodes))' .* nodes[:, i+1], 0:ν(nodes)
    )

    """
        In place bezier curve computation, see `bezier` for more details.
    """
    function bezier!(curve_points::Points, nodes::Points; δt::Float64=.01, closed::Bool=false)::Curve
        # prepare parameters
        nodes, n, τ = prep_spline_parameters(nodes, δt, closed)

        # compute bezier curve
        B(i) = bernstein(τ; i=i, n=n)' .* nodes[:, i+1]
        curve_points = sum(B, 0:n)

        return Curve(
            "bezier",
            nodes,
            τ,
            curve_points,
            (t) -> bezier(t; nodes=nodes)
        )
    end


    # ---------------------------------------------------------------------------- #
    #                                Rational bezier                               #
    # ---------------------------------------------------------------------------- #

    """
        rational_bezier(nodes::Points, weights::Vector{Float64}; δt::Float64=.01, closed::Bool=false)::Points
    
    Computes a rational bezier curve: https://en.wikipedia.org/wiki/B%C3%A9zier_curve#Rational_B%C3%A9zier_curves
    Similar to `bezier` but with an additional 1xN weights vector.
    """
    function rational_bezier(nodes::Points, weights::Vector{Float64}; δt::Float64=.01, closed::Bool=false)::Curve
        return rational_bezier!(init_curve(nodes, δt), nodes, weights; δt=δt, closed=closed)
    end

    """
        rational_bezier(t, X)

    Evaluate a rational_bezier function at a given parameter value `t` and given a set
    of control nodes `X`
    """
    function rational_bezier(t::Number; nodes::Points, weights::Vector{Float64})
        n = ν(nodes) # number of control points

        numerator(i) = bernstein(t; i=i, n=n)' .* nodes[:, i+1] .* weights[i+1]
        denominator(i) = bernstein(t; i=i, n=n)' .* weights[i+1]
        return sum(numerator, 0:n) ./ sum(denominator, 0:n)
    end

    """
        In-place rational bezier curve computation, see `rational_bezier`
    """
    function rational_bezier!(curve_points::Points, nodes::Points, weights::Vector{Float64}; δt::Float64=.01, closed::Bool=false)::Curve
        # prepare parameters
        nodes, n, τ = prep_spline_parameters(nodes, δt, closed)

        numerator(i) = bernstein(τ; i=i, n=n)' .* nodes[:, i+1] .* weights[i+1]
        denominator(i) = bernstein(τ; i=i, n=n)' .* weights[i+1]
        curve_points = sum(numerator, 0:n) ./ sum(denominator, 0:n)

        return Curve(
            "rational_bezier",
            nodes,
            τ,
            curve_points,
            (t) -> rational_bezier(t; nodes=nodes, weights=weights)
        )
    end


end

