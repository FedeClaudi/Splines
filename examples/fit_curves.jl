# ? activate Splines
import Pkg
Pkg.activate("Splines")


using PlotlyJS
using Revise
Revise.revise()


import Splines
import Splines: fit

"""
    Fit a curve to a point-cloud in 3D Euclidean space.

    The curve is defined by a set of 'nodes', point in 3D space. Thes nodes
    are initialized as the centroids of clusters identified through kMeans 
    clustering of the data points and their position is adjusted to minimize
    an error (lenght of the resulting curve + fit to the data).
    The curve is then done through piecewise interpolation of consecutive pairs of data.
"""
ENV["JULIA_DEBUG"]="all"


# initialize raw data points
data = Splines.TestData.circle3D(σ=.1, δ=.001)

# fit
nodes_init, nodes_optim, curve, opt_res = fit(
            data, 
            :bezier;  # type of curve to fit: :PiecewiseLinear, :bspline, :bezier
            n=20, # number of nodes
            closed=true,
            α=0.1, β=5.0,
            )

print(opt_res)


# plot results
Splines.plot_fit_results(
    data, curve, nodes_optim, nodes_init
)
@info "Done, happy days!"