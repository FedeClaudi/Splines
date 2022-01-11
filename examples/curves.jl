# ? activate Splines
import Pkg
Pkg.activate("Splines")

using PlotlyJS
using Revise
Revise.revise()

import Splines: BSpline, PiecewiseLinear, Bezier


"""
    This example shows how to create different kind of curves
    given a set of control points.

    Curves:
        - PiecewiseLinear
        - BSpline
        - Bezier
"""

# define control points
X = [[3, 1, 0] [2.5, 3, .2] [0, 4, .6] [-2.5, 3, 1] [-1, 0, 1.4] [-2.5, -2, 2] [0, -1, 2.2] [2.5, -3, 3] [3, -1, 3.5]]

# spline curve of 3d degree
spline = BSpline(X, d=3)

# piecewise linear curve
pwl = PiecewiseLinear(X)

# Bezier curve
bezier = Bezier(X)

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
    scatter(x=spline[1, :], y=spline[2, :], z=spline[3, :], mode="lines", type="scatter3d",
            line = attr(color="#D81B60", width=12), name="spline"),

    scatter(x=pwl[1, :], y=pwl[2, :], z=pwl[3, :], mode="lines", type="scatter3d", 
            line = attr(color="#546E7A", width=12), name="PWL"),

    scatter(x=bezier[1, :], y=bezier[2, :], z=bezier[3, :], mode="lines", type="scatter3d",
            line = attr(color="#F4511E", width=12), name="bezier"),
],  Layout(
    xaxis=attr(
        nticks=4,
        range=[0, 255],
        axis=([], false),
        name="R",
    ),
    yaxis=attr(
        nticks=4,
        range=[0, 255],
        axis=([], false),
    ),
    zaxis=attr(
        nticks=4,
        range=[0, 255],
        axis=([], false),
    ),
)))
@info "Done!"
