# ? activate Splines
import Pkg
Pkg.activate("Splines")

using PlotlyJS
using Revise
Revise.revise()

using Splines

n, m = 8, 12
nodes1 = cat([0:n-1 zeros(n) sin.(range(0, stop=π, length=n))], dims=2)'
nodes2 = cat([zeros(m) 0:m-1 sin.(range(0, stop=π, length=m))], dims=2)'

# create two bezier curves
b1 = @time Bezier(nodes1)
b2 = @time Bezier(nodes2)


# create bezier surface
bb = b1 ⊗ b2



display(
    plot([
        plot_nodes(nodes1, color="red"),
        plot_curve(b1.coordinates, color="red"),
        plot_nodes(nodes2, color="blue"),
        plot_curve(b2.coordinates, color="blue"),
        plot_surface(bb.coordinates)...
]))