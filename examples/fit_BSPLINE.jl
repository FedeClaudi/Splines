
# ? activate Splines
import Pkg
Pkg.activate("Splines")


using Plots
using Revise
Revise.revise()


import Splines
import Splines: fit

"""
    Fit a b-spline curve to a point-cloud in 3D Euclidean space.

    The curve is defined by a set of 'nodes', point in 3D space. Thes nodes
    are initialized as the centroids of clusters identified through kMeans 
    clustering of the data points and their position is adjusted to minimize
    an error (lenght of the resulting curve + fit to the data).
    The curve is then done through piecewise interpolation of consecutive pairs of data.
"""

gr()

# initialize raw data points
data = Splines.TestData.circle3D(σ=.1, δ=.001)

# fit
nodes_init, nodes_optim, curve = fit(
            data, 
            :BSline;  # type of curve to fit
            n=10, # number of nodes
            d=3,  # dimensionality of the spline
            )


# plot results
# points and curve
plt = scatter3d(data[1, :], data[2, :], data[3, :], label="data", color="white", camera=(30, 70))
plot3d!(curve[1, :], curve[2, :], curve[3, :], label="curve", color="green", lw=4)

# nodes position
scatter3d!(nodes_optim[1, :], nodes_optim[2, :], nodes_optim[3, :], label="optim. knoblackts", color="red", ms=10)
scatter3d!(nodes_init[1, :], nodes_init[2, :], nodes_init[3, :], label="init. nodes", color="black", ms=7)

display(plt)
@info "Done, happy days!"