using Plots
using Revise
Revise.revise()

# ? activate Splines
import Pkg
Pkg.activate("Splines")


import Splines
import Splines.Fit: fitPWL as fit

"""
    Fit a piecewise linear curve to a circle in 3D Euclidean space
"""

gr()

# initialize raw data points
data = Splines.TestData.circle3D(σ=.1, δ=.001)

# fit
knots_init, knots_optim, curve = fit(data; η=20, closed=true, α=1.0)


# # plot stuff
plt = scatter3d(data[1, :], data[2, :], data[3, :], label="data", color="white", camera=(30, 70))
plot3d!(curve[1, :], curve[2, :], curve[3, :], label="curve", color="green", lw=4)

scatter3d!(knots_optim[1, :], knots_optim[2, :], knots_optim[3, :], label="optim. knoblackts", color="red", ms=10)
scatter3d!(knots_init[1, :], knots_init[2, :], knots_init[3, :], label="init. knots", color="black", ms=7)
scatter3d!([knots_init[1, 1]], [knots_init[2, 1]], [knots_init[3, 1]], label=nothing, color="blue", ms=5)


display(plt)
@info "Done, happy days!"