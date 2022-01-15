import Pkg
Pkg.activate("Splines")

using PlotlyJS
using BenchmarkTools

using Revise
Revise.revise()

import Splines: BezierSurface
import Splines.Visuals: plot_surface

# define nodes grid
n = 4
m = 4

x = (1:n) ./ n .* ones(m)'  # n x m
y = (1:m)' ./ m .* ones(n)  # n x m

ϕ₁, ϕ₂ = 1, 1  # sines frequency
z = (sin.(ϕ₁ .* range(0, stop=pi, length=m)) .* sin.(ϕ₂ .* range(0, stop=pi, length=n))')'  # n x m

X = cat(x, y, z, dims=3)  # n x m x 3
X = permutedims(X, [3, 1, 2]) # 3 x n x m


bezier = BezierSurface(X)


display(
    plot(
        plot_surface(bezier)
    )
)