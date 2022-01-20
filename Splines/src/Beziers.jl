module Beziers
    include("Utils.jl")

    
    import ..Spline: AbstractBezier
    import ..Polynomials: bernstein
    import ..Utils: get_curve_parameters, ∑

    export Bezier



    """
        Bezier

    Struct of type Bezier <: AbstractBezier.
    Contains all information about an N-dimensional Bezier surface.

    For `deg=1` we have a Bezier curve and:
        `nodes`: `dxN` array of nodes coordinates
        `η`: dimensionality of the nodes (`N-1`)
        `coordinates`: d x K points with the coordinates of the Bezier curve.

    for `deg>1` the type represents `deg` dimensional Bezier surfaces
    """
    struct Bezier <: AbstractBezier
        nodes::AbstractArray  # coordinates of nodes (d x N (x M ...))
        coordinates::AbstractArray   # coordinates of points along the spline 
        weights::AbstractArray
        η  # shape of the nodes -1  | degree (linera, quadratic, cubic..)
        d::Int  # dimensionality of the Euclidean space
        N::Int  # spline dimensionality: 1 - curve, 2 - surface...
        inner_fn  # a collection of lower dimensionality bezier curves
    end
    

    """
        Bezier(nodes; δt=0.05)

    Constructor for a degree 1 (`n=1`) `Bezier` curve defined by a set of nodes
    """
    function Bezier(nodes; δt=0.025)::Bezier
        τ, N, η, d = get_curve_parameters(nodes, δt)
        weights = ones(size(nodes)...)

        # compute coordinates
        B(i) = bernstein(τ; i=i, n=η)' .* nodes[:, i+1]
        coordinates = ∑(B, 0:η)

        return Bezier(nodes, coordinates, weights, η, d, N, nothing)
    end


    struct RationalBezier <: AbstractBezier
        nodes::AbstractArray  # coordinates of nodes (d x N (x M ...))
        coordinates::AbstractArray   # coordinates of points along the spline 
        weights::AbstractArray  # of shape d x N x M...
        η  # shape of the nodes -1  | degree (linera, quadratic, cubic..)
        d::Int  # dimensionality of the Euclidean space
        N::Int  # spline dimensionality: 1 - curve, 2 - surface...
        inner_fn  # a collection of lower dimensionality bezier curves
    end

    function RationalBezier(nodes, weights; δt = 0.025)::RationalBezier
        τ, N, η, d = get_curve_parameters(nodes, δt)

        numerator(i) = bernstein(τ; i=i, n=η)' .* nodes[:, i+1] .* weights[i+1]
        denominator(i) = bernstein(τ; i=i, n=η)' .* weights[i+1]
        coordinates = ∑(numerator, 0:η) ./ ∑(denominator, 0:η)

        return RationalBezier(
            nodes, coordinates, weights, η, d, N, nothing
        )
    end

    """
        Bezier(t) or RationalBezier(t)

    Turn an instance of `Bezier` into a callable.
    Evaluates a bezier curve (N=1) at parameter value `t ∈ [0,1]`
    """
    function (bez:: AbstractBezier)(t)
        numerator(i) = bernstein(t; i=i, n=bez.η)' .* nodes[:, i+1] .* weights[i+1]
        denominator(i) = bernstein(t; i=i, n=bez.η)' .* weights[i+1]

        ∑(numerator, 0:bez.evaluates) ./ ∑(denominator, 0:bez.evaluates)
    end







end