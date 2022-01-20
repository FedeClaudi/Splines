
module Maths

    using EllipsisNotation

    import ..Spline: AbstractBezier

    export ⊗, dimensions


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
            return x .* y'
        else
            # initialize an array of the appropriate shape
            N = zeros(3, dimensions(x) ..., dimensions(y)...)
            
            for Ix in CartesianIndices(x[1, ..]), Iy in CartesianIndices(y[1, ..])
                N[:, Ix, Iy] = x[:, Ix] .+ y[:, Iy]
            end
        end
        return N
    end


    """
        ⊗(b1::AbstractBezier, b2::AbstractBezier)::AbstractBezier

    Tensor product of two beziers
    """
    function ⊗(b1::T, b2::T) where T <: AbstractBezier
        # get nodes
        nodes = b1.nodes ⊗ b2.nodes
        n, m = nnodes(nodes)

        # get weights
        weights = b1.weights ⊗ b2.weights
        
        # get coordinates
        coordinates = b1.coordinates ⊗ b2.coordinates

        # concatenate inner functions
        β1 = b1.inner_fn == nothing ? [b1] : b1.inner_fn
        β2 = b2.inner_fn == nothing ? [b2] : b2.inner_fn
                
        return T(nodes, coordinates, weights, (n, m), b1.d, b1.N+b2.N, hcat(β1, β2))
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