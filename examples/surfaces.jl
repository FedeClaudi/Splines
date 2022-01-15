import Pkg
Pkg.activate("Splines")

using PlotlyJS
using BenchmarkTools

using Revise
Revise.revise()

import Splines: bezier_surface, bspline_surface
import Splines.Visuals: plot_surface


# define nodes grid
n = 8
m = 8

x̂ = (1:n) ./ n .* ones(m)'  # n x m
ŷ = (1:m)' ./ m .* ones(n)  # n x m

ϕ₁, ϕ₂ = 1, 1  # sines frequency
ẑ = (sin.(ϕ₁ .* range(0, stop=pi, length=m)) .* sin.(ϕ₂ .* range(0, stop=pi, length=n))')'  # n x m

X = cat(x̂, ŷ, ẑ, dims=3)  # n x m x 3
X = permutedims(X, [3, 1, 2]) # 3 x n x m

# compute surfaces
bezier_surf = bezier_surface(X; δt=.025)
bspline_surf = bspline_surface(X; d1=5, d2=1, δt=.025)


# plot surfaces
fig = make_subplots(
    rows=2, cols=1,
    specs=fill(Spec(kind="scene"), 2, 1),
    subplot_titles=["Bezier" "B Spline"],
)



map((tr) -> add_trace!(
        fig,
        tr,
        row=1, col=1
    ), plot_surface(bezier_surf))

map((tr) -> add_trace!(
        fig,
        tr,
        row=2, col=1
    ), plot_surface(bspline_surf))



display(fig)
@info "done!"