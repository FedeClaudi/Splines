
module Utils
    import StaticArrays: @MMatrix

    include("Types.jl")
    include("Geometry.jl")

    import .Types: Points, asPoints
    import .Geometry: distances


    # ------------------------- curve creation helper fns ------------------------ #
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



    # ---------------------------------- sorting --------------------------------- #
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
        sorted = asPoints(zeros(size(points)))
        
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