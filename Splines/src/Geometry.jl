module Geometry

    export Point, Points, distances, curve_length

    const Point = Vector{Float64}
    const Points = AbstractArray{Float64}


    # --------------------------------- DISTANCES -------------------------------- #
    """
        distances(points::Points)::Points

    Given a set of {d x N} points it computes the distance
    between each consecutive pair of points and returns an
    (N-1)-vector of distance values.
    The distance is computed as euclidean distance:
    `d = √(∑Δxᵢ)`
    where `Δxᵢ` is the `diff` if the i-th dimension of the points array.

    See also [`curve_length`](@ref),
    """
    distances(points::Points)::Vector = vec(sum(
        sqrt.(
            diff(points, dims=2).^2
            ),
        dims=1
    ))

    """
        distances(points::Points, point::Point)::Points

    When an additional `point` is given, the distances represent the distance 
    of each point in `points` to the reference `point`
    """
    function distances(points::Points, point::Point)::Vector
        return vec(sqrt.(sum((points .- point).^2, dims=1)))
    end

    """
        curve_length(points)

    Given a set of {d x N} points, returns the length of the curve
    going through them
    """
    curve_length(points::Points) = sum(distances(points))


end

