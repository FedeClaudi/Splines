module Geometry
    include("Utils.jl")
    import .Utils: dropCol

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


    # ----------------------------------- misc ----------------------------------- #    

    """
        first(points)

    Just one, the index of the first element
    """
    first(points::Points)::Int = 1

    """
        smallest(points)

    Index the point with smallest magnitude
    """
    smallest(points::Points)::Int = argmin(sum(points.^2, dims=1)')[1]


    """
        sort_points(points::Points)

    Sorts a set of `Points`, starts with a single Point, gets the next
    closest one and so on.
    """
    function sort_points(points::Points; selection_method::Symbol=:first)::Points
        #initialize
        sorted = Points(zeros(size(points)))
        
        # select a staring point
        sel_idx = eval(:($selection_method($points)))
        points = [points[:, sel_idx] dropCol(points, sel_idx)]
        sorted[:, 1] = points[:, 1]
        remaining = points[:, 2:end]

        # fill in
        while size(remaining, 2) > 0
            # get next closest
            next_idx = argmin(distances(remaining, sorted[:, end-size(remaining, 2)]))
            sorted[:, end-size(remaining, 2)+1] = remaining[:, next_idx]
            # remove from remaining
            remaining = remaining[:, 1:end .!= next_idx]
        end
    
        return sorted
    end
    
end

