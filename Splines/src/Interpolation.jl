module Interpolation
    include("Geometry.jl")
    import .Geometry: Point, Points

    export lerp, PiecewiseLinear, PiecewiseLinear
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
        PiecewiseLinear(points[; η=100])

    Piecewise linear interpolation between consecutive pairs of points (d x N, `Points`) ∈ ℝⁿ to return
    a curve connecting them. The curve is specified by (η*N)-many points (d x η*N Points array), the 
    parameter η specifies the number of points in each linear segment.
    """
    function PiecewiseLinear(points::Points; η::Int64=100, closed::Bool=false)::Points
        # compute quantities
        N = size(points, 2)
        n_points = closed ? N : N - 1
        ηₜ = η * n_points
        curve = zeros(2, ηₜ)
        return PiecewiseLinear!(curve, points; η=η, closed=closed)
    end


    """
        PiecewiseLinear!(curve, points[; η=100])

    Piecewise linear interpolation between consecutive pairs of points, given a
    pre-allocated curve array (d x η*N Points array).
    """
    function PiecewiseLinear!(curve::Points, points::Points; η::Int64=100, closed::Bool=false)::Points
        # compute quantities
        N = size(points, 2)
        n_points = closed ? N : N - 1
        P = range(0, 1, length=η)  # for interpolation

        # @info "PiecewiseLinear! with curve $(size(curve)), points $(size(points)); η=$η, closed=$closed | N=$N"
        
        # create curve by linear interpolation of each segment
        for n in range(1, length=n_points)
            if n == N
                @info "ops"
                k₀, k₁ = points[:, end], points[:, 1]
            else
                k₀, k₁ = points[:, n], points[:, n+1]
            end
            curve[:, η * (n - 1) + 1: η * n] = lerp(k₀, k₁, P)'
        end
    
        return curve
    end

end