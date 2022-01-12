# ? activate Splines
import Pkg
Pkg.activate("Splines")

using PlotlyJS
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
weights[end-2] = 2
weights[end-1] = -1
rbezier = RationalBezier(X, weights)

# plot
display(plot([
    scatter(x=X[1, :], y=X[2, :], z=X[3, :], mode="markers", type="scatter3d",     
        marker=attr(
            size=8,
            color="black",
            edgecolor="black",
            opacity=1,
        ), name="nodes"
    ),
    scatter(x=spline.points[1, :], y=spline.points[2, :], z=spline.points[3, :], mode="lines", type="scatter3d",
            line = attr(color="#D81B60", width=12), name="spline"),

    scatter(x=pwl.points[1, :], y=pwl.points[2, :], z=pwl.points[3, :], mode="lines", type="scatter3d", 
            line = attr(color="#546E7A", width=12), name="PWL"),

    scatter(x=bezier.points[1, :], y=bezier.points[2, :], z=bezier.points[3, :], mode="lines", type="scatter3d",
            line = attr(color="#F4511E", width=12), name="bezier"),

    scatter(x=rbezier.points[1, :], y=rbezier.points[2, :], z=rbezier.points[3, :], mode="lines", type="scatter3d",
            line = attr(color="blue", width=12), name="Rational bezier"),
],  Layout(
    scene=attr(        

        xaxis=attr(
            nticks=3,

        ),
        yaxis=attr(
            nticks=3,

        ),
        zaxis=attr(
            nticks=3,

        ),
    ),
)))
@info "Done!"


