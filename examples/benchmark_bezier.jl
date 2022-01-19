# ? activate Splines
import Pkg
Pkg.activate("Splines")

using Revise
Revise.revise()

using Splines

n, m, p = 8, 12, 10
nodes1 = cat([0:n-1 zeros(n) sin.(range(0, stop=π, length=n))], dims=2)'
nodes2 = cat([zeros(m) 0:m-1 sin.(range(0, stop=π, length=m))], dims=2)'
nodes3 = cat([0:p-1 .05 .* 0:(p-1) .2 .* (0:-1:(-p+1))], dims=2)'

# create two bezier curves
@info "N=1"
b₁ = Bezier(nodes1, δt=0.05)
b₂ = @time Bezier(nodes2, δt=0.05)
b₃ = Bezier(nodes3, δt=0.05)


# create bezier surface
@info "N=2"
b2n = b₁ ⊗ b₂
@time b₁ ⊗ b₂;

# create a Bezier volume!
@info "N=3"
b3n =  b2n ⊗ b₃
@time b2n ⊗ b₃;

@info "N=5"
b5n = @time b2n ⊗ b2n




