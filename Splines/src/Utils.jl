
module Utils

    include("Geometry.jl")
    import .Geometry: Points, distances

    export dropCol


    """
        dropCol(X::AbstractArray, idx::Int)

    Drops a column in an array, specified by an `idx` value
    """
    dropCol(X::AbstractArray, idx::Int)::AbstractArray = X[:, filter(x->x!=idx, 1:end)]


    """
    Drops a set of columns specified by an array of indices
    """
    dropCol(X::AbstractArray, idx::Vector{Int})::AbstractArray = X[:, filter(x->!(x ∈ idx), 1:end)]



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