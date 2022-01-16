module Surfaces
    include("Types.jl")
    include("Polynomials.jl")
    include("Utils.jl")

    using .Types
    import .Polynomials: bernstein, N
    import .Utils: prep_surface_parameters, uniform, periodic

    export bezier_surface, bspline_surface

    function bezier_surface(nodes::Points; δt = 0.01)::Surface
        # get params
        n, m, τ, T = prep_surface_parameters(nodes, δt)

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
            "bezier Surface",
            nodes,
            τ .* τ',
            surface_coords
        )
    end



    function bspline_surface(nodes::Points; δt = 0.01, d1::Int=3, d2::Int=3, knots_type::Symbol=:uniform)
        # get params
        n, m, τ, T = prep_surface_parameters(nodes, δt)
        k_n = eval(:($knots_type($n, $d1)))
        k_m = eval(:($knots_type($m, $d2)))
        
        # initialize array
        surface_coords = zeros(3, length(τ), length(τ));

        # compute surfface points
        for i in 0:n, j in 0:m
            β = N(τ; k=k_n, i=i, d=d1) * N(τ; k=k_m, i=j, d=d2)'
            for dim in 1:3
                surface_coords[dim, :, :] .+= β .* nodes[dim, i+1, j+1]
            end
        end

        return Surface(
            "b-spline Surface",
            nodes,
            τ .* τ',
            surface_coords
        )
    end

end