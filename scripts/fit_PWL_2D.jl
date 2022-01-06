using Plots


# ? activate Splines
import Pkg
Pkg.activate("Splines")

using Revise
Revise.revise()

import Splines
import Splines.Fit: fitPWL as fit

"""
    Fit a piecewise linear curve to a circle in the 2D plane
"""

gr()

# initialize raw data points
data = Splines.TestData.sine(δ=0.001)

# fit
knots_init, knots_optim, curve = fit(data; η=8, closed=false, α=10.0, β=1.0)
@info "Timed fitting"
@time fit(data; η=12, closed=false)

# plot stuff
plt = scatter(data[1, :], data[2, :], label="data", color="white")
plot!(curve[1, :], curve[2, :], label="curve", color="green", lw=4)

scatter!(knots_optim[1, :], knots_optim[2, :], label="optim. knoblackts", color="red", ms=10)
scatter!(knots_init[1, :], knots_init[2, :], label="init. knots", color="black", ms=7)
scatter!([knots_init[1, 1]], [knots_init[2, 1]], label=nothing, color="blue", ms=5)


display(plt);
@info "Done, happy days!"