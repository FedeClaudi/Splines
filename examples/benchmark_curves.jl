# ? activate Splines
import Pkg
Pkg.activate("Splines")

using BenchmarkTools
using Revise
Revise.revise()

import Splines: bspline, bezier, rational_bezier, rational_bspline


"""
    This example shows how to create different kind of curves
    given a set of control points.

    Curves:
        - bspline
        - bezier
"""

# define control points
X = [[3, 1, 0] [2.5, 3, .2] [0, 4, .6] [-2.5, 3, 1] [-1, 0, 1.4] [-2.5, -2, 2] [0, -1, 3.2] [2.5, -3, 4.3] [3, -4, 5] [1, -5, 7] [-3, -1, 8] [-1.2, 1, 8.9]]

# call each functiuon first to ensure it's compiled
bspline(X, d=3)
bezier(X)
weights = ones(size(X, 2))
weights[3] = -1
rational_bezier(X, weights)

# benchmark
# @info "Benchmarking B spline"
# display(@benchmark bspline(X, d=3))

# @info "Benchmarking Rational bspline"
# display(@benchmark rational_bspline(X, weights; d=3))


# @info "Benchmarking bezier"
# display(@benchmark bezier(X))

# @info "Benchmarking Rational bezier"
# display(@benchmark bezier(X))


# time
@info "Benchmarking B spline"
@btime bspline(X, d=3); print("\n")

@info "Benchmarking rational B spline"
@btime rational_bspline(X, weights; d=3); print("\n")

@info "Benchmarking bezier"
@btime bezier(X); print("\n")

@info "Benchmarking Rational bezier"
@btime bezier(X); print("\n")
