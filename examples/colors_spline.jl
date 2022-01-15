
# ? activate Splines
import Pkg
Pkg.activate("Splines")

using PlotlyJS
using Revise

import Splines: bspline

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
for points in [nodes, r_nodes, b_nodes, g_nodes]
    spline = bspline(points, d=3; δt=0.005)
    plotted_line = scatter(
        x=spline.points[1, :],
        y=spline.points[2, :],
        z=spline.points[3, :],
        mode="markers",
        type="scatter3d",
        line = attr(color=spline.points, width=20), name=nothing)
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