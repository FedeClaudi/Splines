# ? activate Splines
import Pkg
Pkg.activate("Splines")

using PlotlyJS
using Revise
Revise.revise()

import Splines: BSpline, Bezier, RationalBezier
import Splines.Visuals: plot_curve, plot_nodes


"""
    This example shows how to create different kind of curves
    given a set of control points.

    Curves:
        - BSpline
        - Bezier
"""

# define control points
X = [[3, 1, 0] [2.5, 3, .2] [0, 4, .6] [-2.5, 3, 1] [-1, 0, 1.4] [-2.5, -2, 2] [0, -1, 3.2] [2.5, -3, 4.3] [3, -4, 5] [1, -5, 7] [-3, -1, 8] [-1.2, 1, 8.9]]

# spline curve of 3d degree
spline = BSpline(X, d=3)

# Bezier curve
bezier = Bezier(X)

# rational bezier curve
weights = ones(size(X, 2))
weights[3] = -1
weights[5] = 2
weights[end-2] = 2
weights[end-1] = -1
rbezier = RationalBezier(X, weights)

# plot
display(
    plot(
        [   
            plot_nodes(X),
            plot_curve(spline; color="#D81B60"),
            plot_curve(bezier; color="#F4511E"),
            plot_curve(rbezier; color="#blue"),
    ],
))
@info "Done!"


