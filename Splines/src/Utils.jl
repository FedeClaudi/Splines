
module Utils
    import StaticArrays: @MMatrix

    include("Types.jl")
    include("Geometry.jl")

    import .Types: Points, asPoints, Knots
    import .Geometry: distances

    # -------------------------------- coordinates ------------------------------- #
    """
        x(X)

    Select the first coordinate from an array d x N ( x M)
    """
    x(X::Points) = ndims(X)==2 ? X[1, :] : X[1, :, :]

    """
        y(X)
        
    Select the second coordinate from an array d x N ( x M)
    """
    y(X::Points) = ndims(X)==2 ? X[2, :] : X[2, :, :]

    """
        z(X)
        
    Select the third coordinate from an array d x N ( x M)
    """
    z(X::Points) = ndims(X)==2 ? X[3, :] : X[3, :, :]

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