module Interpolation
    include("Geometry.jl")
    import .Geometry: Point, Points

    export lerp, PiecewiseLinear
    """
        lerp(x₀, x₁, p)

    Linear interpolation between two values, given an interpolation factor
    `p ∈ [0, 1]`. `x₀` and `x₁` can be number of `Point`, but they need to have the same
    Type and shape (if `Point`).
    """

    lerp(x₀::Point, x₁::Point, p::Float64)::Point = x₀ * (1 - p) + x₁ * p
    lerp(x₀::Number, x₁::Number, p::Float64)::Number = x₀ * (1 - p) + x₁ * p

    """
        PiecewiseLinear(points[; η=100])

    Piecewise linear interpolation between consecutive pairs of points (d x N, `Points`) ∈ ℝⁿ to return
    a curve connecting them. The curve is specified by (η*N)-many points (d x η*N Points array), the 
    parameter η specifies the number of points in each linear segment.
    """
    function PiecewiseLinear(points::Points; η::Int64=100)::Points
        # compute quantities
        n_points = size(points, 2) - 1
        ηₜ = η * n_points
        P = range(0, 1, length=η)  # for interpolation
        
        # create curve by linear interpolation of each segment
        curve = zeros(2, ηₜ)
        @info size(curve)
        for n in range(1, length=n_points)
            k₀, k₁ = points[:, n], points[:, n+1]
            for (nₚ, p) in enumerate(P)
                curve[:, nₚ + η * (n-1)] = lerp(k₀, k₁, p)
            end
        end
    
        return curve
    end
end