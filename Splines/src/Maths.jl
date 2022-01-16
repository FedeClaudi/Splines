
module Maths
    using EllipsisNotation

    export ∑, ⊗, ϕ!

    # shorthand for sum operation
    ∑(fn, values) = sum(fn, values)
    

    """
    Tensor product between two arrays.
    """
    function ⊗(x::AbstractArray, y::AbstractArray)::AbstractArray
        if ndims(x) == ndims(y) == 1
            return x * y'
        else
            if ndims(x) > 2 || ndims(y) > 2
                @warn "This is only tested with ndims=2"
            end

            # initialize an array of the appropriate shape
            n, m = size(x, 2), size(y, 2)
            N = zeros(3, n, m)
            
            # compute product for each Euclidean dimension
            ϕ!(N, (d)->repeat(x[d, :, :], 1, m) .+ y[d, :, :]')

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
        Get the number of nodes (given either a d x N or d x N x M input array).
    """
    nnodes(X::AbstractArray) = ndims(X)==2 ? size(X, 2) - 1 : size(X)[2:end] .- 1


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