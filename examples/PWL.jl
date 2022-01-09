# ? activate Splines
import Pkg
Pkg.activate("Splines")

using PlotlyJS
using Revise
Revise.revise()

import Splines: PiecewiseLinear

"""
    This script shows how to construct a b-spline of a given degree given a set of control points.
"""

# define control points
X = [[3, 1, 0] [2.5, 4, .2] [0, 6, .1] [-2.5, 4, .2] [-1, 0, -.1] [-2.5, -4, 0] [0, -1, -.1] [2.5, -4, -.05] [3, -1, 0]]



# fit a piecewise linear curve for comparison
pwl = PiecewiseLinear(X; closed=true)

# plot
display(plot([
    scatter(x=X[1, :], y=X[2, :], z=X[3, :], mode="markers", type="scatter3d",     
        marker=attr(
            size=8,
            color="black",                # set color to an array/list of desired values
            edgecolor="black",
            opacity=1,
        ), name="nodes"
    ),
    scatter(x=pwl[1, :], y=pwl[2, :], z=pwl[3, :], mode="lines", type="scatter3d", color="green", 
            line = attr(color="royalblue", width=4), name="curve"),
]))
@info "Done!"
