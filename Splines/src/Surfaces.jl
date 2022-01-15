module Surfaces
    include("Types.jl")
    include("Polynomials.jl")

    using .Types
    import .Polynomials: bernstein

    export BezierSurface

    function BezierSurface(nodes::Points; δt = .01)::Surface
        # get parameters
        n, m  = size(nodes)[1+1:end] .- 1
        
        τ = 0:δt:1
        T = length(τ)

        # initialize array
        surface_coords = zeros(3, T, T)

        # compute surface points coordinates values
        for i in 0:n, j in 0:m
            # compute matrix of bernstein values    
            β = bernstein(τ; i=j, n=m) * bernstein(τ; i=i, n=n)'

            # multiply β for each dimension of the control nodes
            for d in 1:3
                surface_coords[d, :, :] .+=  β .* nodes[d, i+1, j+1]
            end
        end

        return Surface(
            "Bezier Surface",
            nodes,
            τ .* τ',
            surface_coords
        )
    end


end