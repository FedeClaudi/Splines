# ? activate Splines
import Pkg
Pkg.activate("Splines")

using PlotlyJS
using Revise
Revise.revise()

using Splines

n, m = 8, 12
nodes1 = cat([0:n-1 zeros(n) sin.(range(0, stop=π, length=n))], dims=2)'
nodes2 = cat([zeros(m) 4 .* (0:m-1) sin.(range(0, stop=π, length=m))], dims=2)'


weights1 = ones(size(nodes1, 2))
weights2 = ones(size(nodes2, 2))
weights1[5] = 3
weights2[5] = -.1

# create two bezier curves
@info "N=1"
b₁ = BSpline(nodes1, weights1, 3)
b₂ = BSpline(nodes2, weights2, 3)


# create bezier surface
@info "N=2"
b2n = b₁ ⊗ b₂


@info "plots"
display(
    plot([
        plot_nodes(nodes1, color="red"),
        plot_curve(b₁.coordinates, color="red"),
        plot_nodes(nodes2, color="blue"),
        plot_curve(b₂.coordinates, color="blue"),
        plot_surface(b2n.coordinates, opacity=1)...,
]))
