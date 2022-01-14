import Pkg
Pkg.activate("Splines")

using PlotlyJS

import Splines.Polynomials.Bernstein

# define nodes grid
n = 8
m = 8

x = (1:n) ./ n .* ones(m)'  # n x m
y = (1:m)' ./ m .* ones(n)  # n x m

ϕ₁, ϕ₂ = 2, 1  # sines frequency
z = (sin.(ϕ₁ .* range(0, stop=pi, length=m)) .* sin.(ϕ₂ .* range(0, stop=pi, length=n))')'  # n x m

X = cat(x, y, z, dims=3)  # n x m x 3
X = permutedims(X, [3, 1, 2]) # 3 x n x m




function BezierSurface(nodes::Array{Float64, 3})::Array{Float64, 3}
    # get parameters
    n, m  = size(nodes)[1+1:end] .- 1
    δt = .01
    τ = Array(0:δt:1)
    T = length(τ)

    # initialize array
    surface_coords = zeros(3, T, T)

    # compute surface points coordinates values
    for i in 0:n
        for j in 0:m
            # compute matrix of Bernstein values
            β = Bernstein(τ; i=j, n=m) * Bernstein(τ; i=i, n=n)'

            # multiply β for each dimension of the control node
            p = nodes[:, i+1, j+1]
            for d in 1:3
                surface_coords[d, :, :] .+=  β .* p[d]
            end
        end
    end

    return surface_coords
end


@time BezierSurface(X)
surface_coords = @time BezierSurface(X)

# plot
@info "Plotting"
display(
    plot([
        scatter(
            x=vec(X[1, :, :]), 
            y=vec(X[2, :, :]), 
            z=vec(X[3, :, :]), 
            mode="markers", type="scatter3d",
            marker=attr(
                    size=4,
                    color="red",
                    edgecolor="black",
                    opacity=1,
                ),
            name = "nodes"
            ),

            surface(
                x=surface_coords[1, :, :], 
                y=surface_coords[2, :, :], 
                z=surface_coords[3, :, :], 
                name="Bezier surface",
                showscale=false
                ),

            scatter(
                x=vec(surface_coords[1, :, :]), 
                y=vec(surface_coords[2, :, :]), 
                z=vec(surface_coords[3, :, :]), 
                mode="markers", type="scatter3d",
                marker=attr(
                        size=1,
                        color="black",
                        edgecolor="black",
                        opacity=1,
                    ),
                name="Bezier surface",
                )
    ]))
@info "Done"