# ? activate Splines
import Pkg
Pkg.activate("Splines")

using PlotlyJS
using Revise
Revise.revise()

import Splines: bspline, bezier, rational_bezier, NURBS_curve
import Splines.Visuals: plot_curve, plot_nodes


"""
    This example shows how to create different kind of curves
    given a set of control points.

    Curves:
        - bspline
        - bezier
"""

# define control points
X = [[3, 1, 0] [2.5, 3, .2] [0, 4, .6] [-2.5, 3, 1] [-1, 0, 1.4] [-2.5, -2, 2] [0, -1, 3.2] [2.5, -3, 4.3] [3, -4, 5] [1, -5, 7] [-3, -1, 8] [-1.2, 1, 8.9]]

# ------------------------------- create curves ------------------------------ #
# spline curve of 3d degree
spline = bspline(X, d=2)

# bezier curve
bezier_curve = bezier(X)

# rational b-spline curve
weights = ones(size(X, 2)) 
weights[1:2:end] .= 2
weights[1+1:2:end-1] .= -.1

rbspline = NURBS_curve(X, weights; d=3)

# rational bezier curve
weights = ones(size(X, 2))
weights[3] = -1
weights[5] = 2
weights[end-2] = 2
weights[end-1] = -1

rbezier = rational_bezier(X, weights)

# ------------------------------ what's a curve ------------------------------ #
"""
    Curve generating function return a Curve object:
"""

println(spline)

"""
    which you can use to compute the value of the spline at a parameter value
    with `t ∈ [0, 1]`
"""

println(spline(.3))
half_way = spline(.5)

# ----------------------------------- plot ----------------------------------- #
display(
    plot(
        [   
            plot_nodes(X),
            plot_curve(spline; color="#D81B60"),
            plot_curve(bezier_curve; color="#F4511E"),
            plot_curve(rbspline; color="green", width=3),
            plot_curve(rbezier; color="blue", width=3),
            scatter(x=[half_way[1]], y=[half_way[2]], z=[half_way[3]], marker=attr(size=14, color="#D81B60"), type="scatter3d")
    ],
))
@info "Done!"


