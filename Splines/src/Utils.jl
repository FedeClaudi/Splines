module Utils

    # shorthand for sum operation
    ∑(fn, values) = sum(fn, values)

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