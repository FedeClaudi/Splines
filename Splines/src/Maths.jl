
module Maths
    using EllipsisNotation

    export ∑, ⊗, dimensions

    # shorthand for sum operation
    ∑(fn, values) = sum(fn, values)
    

    """
        ⊗(x::AbstractArray, y::AbstractArray)
    
    Tensor product between two arrays.
    Given x as a d x N x M... and y as a d x P x Q...
    it returns an array of shape d x N x M .. x P x Q x ...
    with the product of the two input arrays.

    This is done by iterating over the dimension d using Einstein summatin.
    Then each N x M... and P x Q... array is repeated/permuted to 
    to create two N x M... x P x Q... arrays which are then summed.
    """
    function ⊗(x::AbstractArray, y::AbstractArray)::AbstractArray
        # initialize an array of the appropriate shape
        N = zeros(3, dimensions(x) ..., dimensions(y)...)
        
        for Ix in CartesianIndices(x[1, ..]), Iy in CartesianIndices(y[1, ..])
            N[:, Ix, Iy] = x[:, Ix] .+ y[:, Iy]
        end
        return N
    end


    # ------------------------------ Misc.utilities ------------------------------ #

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
end