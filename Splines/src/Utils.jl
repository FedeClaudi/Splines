module Utils

    using EllipsisNotation

    export dimensions, ∑

    # shorthand for sum operation
    ∑(fn, values) = sum(fn, values)


    """
    Gets the dimensions of an array of shape d x N x M ... ignoring
    the first dimension that here represents the 3 cartesian coordinates.
    E.g., d x N x M -> N x M
    """
    dimensions(X::AbstractArray) =  size(X)[2:end]


    """
    Get the number of nodes along each dimension (given d x N x M input array
    it returns (N-1)x(M-1)).
    """
    nnodes(X::AbstractArray) = ndims(X)==2 ? size(X, 2)-1 : dimensions(X) .- 1


    """
        x(X)

    Select the first coordinate from an array d x N ( x M)
    """
    x(X::AbstractArray) = X[1, ..]

    """
        y(X)
        
    Select the second coordinate from an array d x N ( x M)
    """
    y(X::AbstractArray) = X[2, ..]

    """
        z(X)
        
    Select the third coordinate from an array d x N ( x M)
    """
    z(X::AbstractArray) = X[3, ..]

    """
    Prepares standard parameters used in the creation of 1D splines
    """
    function get_curve_parameters(nodes::AbstractArray, δt)
        # parameter interval
        τ = 0:δt:1

        # dimesions
        N = 1
        η = size(nodes, 2) - 1
        d = size(nodes, 1)

        return τ, N, η, d
    end
end