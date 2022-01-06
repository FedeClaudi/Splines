using Plots
import MultivariateStats: fit as fit_PCA
import MultivariateStats: PCA, transform

using Revise
Revise.revise()

# ? activate Splines
import Pkg
Pkg.activate("Splines")


import Splines
import Splines.Fit: fitPWL as fit

"""
    Fit a piecewise linear curve to a circle in 4D Euclidean space
"""

gr()

# initialize raw data points
data = Splines.TestData.circle4D(σ=.1, δ=.001)

# fit
knots_init, knots_optim, curve = @time fit(data; η=20, closed=true, α=1.0)


# project to 2D with PCA
M = fit_PCA(PCA, data, maxoutdim=2)

ŷ = transform(M, data)
k̂ᵢ = transform(M, knots_init)
k̂ₒ = transform(M, knots_optim)
c = transform(M, curve)


# plot stuff
plt = scatter(ŷ[1, :], ŷ[2, :], label="data", color="white", camera=(30, 70))
plot!(c[1, :], c[2, :], label="curve", color="green", lw=4)

scatter!(k̂ₒ[1, :], k̂ₒ[2, :], label="optim. knoblackts", color="red", ms=10)
scatter!(k̂ᵢ[1, :], k̂ᵢ[2, :], label="init. knots", color="black", ms=7)
scatter!([k̂ᵢ[1, 1]], [k̂ᵢ[2, 1]], label=nothing, color="blue", ms=5)


display(plt)
@info "Done, happy days!"