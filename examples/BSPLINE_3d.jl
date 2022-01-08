# ? activate Splines
import Pkg
Pkg.activate("Splines")

using Plots
using Revise
Revise.revise()

import Splines: BSpline, PiecewiseLinear, PiecewiseLinear!


"""
    This script shows how to construct a b-spline of a given degree given a set of control points.
"""

# define control points
X = [[3, 1, 0] [2.5, 4, .2] [0, 6, .1] [-2.5, 4, .2] [-1, 0, -.1] [-2.5, -4, 0] [0, -1, -.1] [2.5, -4, -.05] [3, -1, 0]]

# compute spline curve of 3d degree and firsdt degree
spline = BSpline(X, d=3);
spline_linear = BSpline(X; d=1)

# fit a piecewise linear curve for comparison
pwl = PiecewiseLinear(X)

# plot
plt = scatter(X[1, :], X[2, :], X[3, :], label="X", ms=8, color="black", ylim=[-5, 5], camera=(30, 80))

plot!(spline[1, :], spline[2, :], spline[3, :], label="d=3", color="blue", lw=3)
plot!(spline_linear[1, :], spline_linear[2, :], spline_linear[3, :], label="d=1", color="red", lw=2)
plot!(pwl[1, :], pwl[2, :], pwl[3, :], label="pwl", color="#000000", lw=6, alpha=.4)

display(plt)
