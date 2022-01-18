
module Maths
    using EllipsisNotation
    using Einsum

    export ∑, ⊗, ϕ!, dimensions

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
        if ndims(x) == ndims(y) == 1
            return x * y'
        else
            # compute product for each Euclidean dimension
            n_x = ndims(x) - 1
            n_y = ndims(y) - 1
            dimsx = dimensions(x) 
            dimsy = dimensions(y)

            # initialize an array of the appropriate shape
            N = zeros(3, dimsx..., dimsy...)
            
            # prepare conuts/indices for array repeats
            counts_x = Int.([ones(n_x)... dimsy...])
            counts_y = Int.([ones(n_y)... dimsx...])
            dims_permutation = [collect(2:n_x+n_y)... 1]

            # define two utility functions 
            α(x) = repeat(x, counts_x...)
            β(y) = permutedims(repeat(y, counts_y...), dims_permutation)

            # @einsum N[d, ..] := α(x[d, ..]) .+ β(y[d, ..])
            ϕ!(N, (d) -> α(x[d, ..]) .+ β(y[d, ..]))

            return N
        end
    end


    """
        ϕ!(X::AbstractArray, fn)

    Apply a function `fn` along the first dimension of a given array `X`.
    `fn` should accept only 1 argument, the index of the dimension slice.
    """
    function ϕ!(X, fn)
        for d in 1:size(X, 1)
            X[d, ..] .= fn(d)
        end
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