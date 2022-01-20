# ? activate Splines
import Pkg
Pkg.activate("Splines")

using PlotlyJS
using Revise
Revise.revise()

using Splines

n, m, p = 8, 12, 10
nodes1 = cat([0:n-1 zeros(n) sin.(range(0, stop=π, length=n))], dims=2)'
nodes2 = cat([zeros(m) 0:m-1 sin.(range(0, stop=π, length=m))], dims=2)'
nodes3 = cat([0:p-1 .05 .* 0:(p-1) .2 .* (0:-1:(-p+1))], dims=2)'

# create two bezier curves
@info "N=1"
b₁ = Bezier(nodes1)
b₂ = Bezier(nodes2)
b₃ = Bezier(nodes3)


# create bezier surface
@info "N=2"
b2n = b₁ ⊗ b₂

# create a Bezier volume!
# @info "N=3"
# b3n =  b2n ⊗ b₃


# plots = []
# N = dimensions(b3n.coordinates)[end]
# for i in 1:N
#     push!(plots, plot_surface(b3n.coordinates[:, :, :, i]; with_points=false, opacity=1 - N/i))
# end

@info "plots"
display(
    plot([
        plot_nodes(nodes1, color="red"),
        plot_curve(b₁.coordinates, color="red"),
        plot_nodes(nodes2, color="blue"),
        plot_curve(b₂.coordinates, color="blue"),
        # plot_curve(b₃.coordinates, color="green"),
        plot_surface(b2n.coordinates, opacity=1)...,
        # plots...
]))

# display(plot(volume(
#     x=vec(b3n.coordinates[1, :, :, :]),
#     y=vec(b3n.coordinates[2, :, :, :]),
#     z=vec(b3n.coordinates[3, :, :, :]),
#     type="scatter3d", 
#     color=vec(b3n.coordinates[3, :, :, :])
#     )))
