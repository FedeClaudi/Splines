
module Utils
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


end