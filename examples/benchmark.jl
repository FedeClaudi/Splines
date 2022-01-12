# ? activate Splines
import Pkg
Pkg.activate("Splines")

using BenchmarkTools
using Revise
Revise.revise()

import Splines: BSpline, PiecewiseLinear, Bezier, RationalBezier


"""
    This example shows how to create different kind of curves
    given a set of control points.

    Curves:
        - PiecewiseLinear
        - BSpline
        - Bezier
"""

# define control points
X = [[3, 1, 0] [2.5, 3, .2] [0, 4, .6] [-2.5, 3, 1] [-1, 0, 1.4] [-2.5, -2, 2] [0, -1, 3.2] [2.5, -3, 4.3] [3, -4, 5] [1, -5, 7] [-3, -1, 8] [-1.2, 1, 8.9]]

# spline curve of 3d degree
spline = BSpline(X, d=3)

# piecewise linear curve
pwl = PiecewiseLinear(X)

# Bezier curve
bezier = Bezier(X)

# rational bezier curve
weights = ones(size(X, 2))
weights[3] = -1
weights[5] = 2
rbezier = RationalBezier(X, weights)

# evaluate performance
@info "Benchmarking piecewise linear"
display(@benchmark PiecewiseLinear(X))

@info "Benchmarking B spline"
display(@benchmark BSpline(X, d=3))

@info "Benchmarking Bezier"
display(@benchmark Bezier(X))

@info "Benchmarking Rational Bezier"
display(@benchmark Bezier(X))