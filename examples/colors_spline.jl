# ? activate Splines
import Pkg
Pkg.activate("Splines")

using PlotlyJS
using Revise

import Splines: BSpline, ⊗
import Splines.Utils: x, y, z

"""
    Fitting splines in R³ to create interpolations in RGB color space.
"""

Revise.revise()


t = 1:6*π
nodes = Array{Float64}(undef, 3, length(t))
nodes[1, :] = (sin.(t) .+ 1) .* 255/2
nodes[2, :] = (cos.(t) .+ 1) .* 255/2
nodes[3, :] = Array(range(0, stop=255, length=size(t, 1)))

r_nodes = [[0, 0, 0] [150, 0, 0] [255, 0, 0] [255, 0, 0]]
b_nodes = [[0, 0, 0] [0, 150, 0] [0, 255, 0] [0, 255, 0]]
g_nodes = [[0, 0, 0] [0, 0, 150] [0, 0, 255] [0, 0, 255]]

lines = []

# plot RGB lines
for points in [r_nodes, b_nodes, g_nodes]

    spline = BSpline(points, 1; δt=0.005)

    plotted_line = scatter(
        x=spline.coordinates[1, :],
        y=spline.coordinates[2, :],
        z=spline.coordinates[3, :],
        mode="lines",
        type="scatter3d",
        line = attr(color=spline.coordinates, width=5), name=nothing)
    push!(lines, plotted_line)
end

# plot splines of different degrees

for degree in (1, 3, 5, 7,)
    spline = BSpline(nodes, degree; δt=0.005)
    plotted_line = scatter(
        x=spline.coordinates[1, :],
        y=spline.coordinates[2, :],
        z=spline.coordinates[3, :],
        mode="lines",
        type="scatter3d",
        line = attr(color=spline.coordinates, width=20), name=nothing)
    push!(lines, plotted_line)
end

# plot
layout = Layout(
    showlegend=false,
    scene=attr(        
            xaxis_title="R",
            yaxis_title="G",
            zaxis_title="B",
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
    ),
)

display(plot([lines...], layout))
@info "Done!"