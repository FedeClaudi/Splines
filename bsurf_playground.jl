import Pkg
Pkg.activate("Splines")

using PlotlyJS
using Revise
Revise.revise()


import Splines.Utils: ν
import Splines.Polynomials: bernstein
import Splines.Visuals: plot_nodes, plot_curve

struct Bezier
    nodes
    τ
    coordinates
    dim
    n
    nodes_flat
    inner
end

function Base.show(io::IO, b::Bezier)
    print(io, "$(b.dim) dimensional Bezier")
end

function Bezier(nodes; δt=0.05)  # constructor for 1 dimensional bezier
    τ = 0:δt:1
    coordinates = sum(
        (i)->bernstein(τ; i=i, n=ν(nodes))' .* nodes[:, i+1], 0:ν(nodes)
    )

    inner(t, i) = bernstein(t; i=i, n=ν(nodes))
    return Bezier(nodes, τ, coordinates, ndims(nodes) - 1, ν(nodes), nodes, inner)
end



function ⊗(n1::AbstractArray, n2::AbstractArray)::AbstractArray
    n, m = size(n1, 2), size(n2, 2)
    N = zeros(3, n, m)
    for d in 1:3
        N[d, :, :] = repeat(n1[d, :, :], 1, m) .+ n2[d, :, :]'
    end
    return N
end


function ⊗(b1::Bezier, b2::Bezier)::Bezier
    # stack nodes by interpolating between original sets
    nodes = b1.nodes ⊗ b2.nodes
    n, m = size(nodes)[2:end] .-1

    # compute coordinates
    coordinates = b1.coordinates ⊗ b2.coordinates
    # coordinates = zeros(3, 101, 101)
    # τ = 0:.01:1

    # for i in 0:n, j in 0:m
    #     β = b1.inner(τ, i) .* b2.inner(τ, j)'
    #     for d in 1:3
    #         coordinates[d, :, :] .+= β .* nodes[d, i+1, j+1]
    #     end
    # end
    # compute propertie
    nodes_flat = reshape(nodes, 3, *(size(nodes)[2:end]...))

    inner(τ, ψ, i, j, n, m) = b1.inner(τ, i, n) * b2.inner(ψ, j, m)
    return Bezier(nodes, b1.τ .* b2.τ', coordinates, b1.dim + b2.dim, ν(nodes), nodes_flat, inner)
end




@info "Computing"
n, m = 8, 8


# cupula
N1 = [
        cat([0:n-1 zeros(n) sin.(range(0, stop=π, length=n))], dims=2)',
        cat([zeros(m) 0:m-1 sin.(range(0, stop=π, length=n))], dims=2)'
    ]


# angled plane
N2 = [
        cat([1:n zeros(n) zeros(n)], dims=2)',
        cat([zeros(m) 1:m 1:m], dims=2)'
     ]

# flat plane
N3 = [
        cat([0:n-1 ones(n).-1 zeros(n)], dims=2)',
        cat([ones(m).-1 0:m-1 zeros(m)], dims=2)'
     ]


# funny shape
N4 = [
        [[-1, -3, 1] [1, 3, -1] [2, 6, -3] [3, 2, -5] [4, -1, 1] [5, 1, 3] [6, 3, 7]],
        [1:7 1:7 1:7]'
     ]


fig = make_subplots(
    rows=2, cols=2,
    specs=fill(Spec(kind="scene"), 2, 2),
)
rows = [1, 1, 2, 2]
cols = [1, 2, 1, 2]

for (n, (nodes_1, nodes_2)) in enumerate([N1, N2, N3, N4])
    b1 = Bezier(nodes_1)
    b2 = Bezier(nodes_2)

    bb = b1 ⊗ b2

    for (nodes, color) in zip((b1.nodes, b2.nodes, bb.nodes), ("red", "blue", "black"))
        add_trace!(fig,
        plot_nodes(nodes; name=nothing, color=color, knot_size=3), row=rows[n], col=cols[n]
        )
    end


    for (b, color) in zip((b1, b2), ("red", "blue"))
        add_trace!(fig,
        
            scatter3d(
                x=b.coordinates[1, :],
                y=b.coordinates[2, :],
                z=b.coordinates[3, :],
                mode="lines",            
                line=attr(width=5, color=color),
                name=nothing,
                ), row=rows[n], col=cols[n]
        )
    end

    add_trace!(
        fig,
        surface(
            x=bb.coordinates[1, :, :], 
            y=bb.coordinates[2, :, :], 
            z=bb.coordinates[3, :, :], 
            showscale=false, 
            name=nothing,
        ),
        row=rows[n], col=cols[n]
    )

    add_trace!(fig,
        plot_nodes(bb.coordinates; name=nothing, color="black", knot_size=1), row=rows[n], col=cols[n]
    )

end

display(fig)
@info "Done"

# b1 = Bezier(N1[1])
# b2 = Bezier(N1[2])

# bb = b1 ⊗ b2
# display(
#     plot([
#         # plot_nodes(b1.nodes; name=nothing, color="red", knot_size=3),
#         # plot_nodes(b2.nodes; name=nothing, color="blue", knot_size=3),
#         scatter3d(
#             x=b1.coordinates[1, :],
#             y=b1.coordinates[2, :],
#             z=b1.coordinates[3, :],
#             mode="lines",            
#             line=attr(width=30, color="red"),
#             name=nothing,
#         ),
#         scatter3d(
#             x=b2.coordinates[1, :],
#             y=b2.coordinates[2, :],
#             z=b2.coordinates[3, :],
#             mode="lines",            
#             line=attr(width=30, color="blue"),
#             name=nothing,
#         ),
#         surface(
#             x=bb.coordinates[1, :, :], 
#             y=bb.coordinates[2, :, :], 
#             z=bb.coordinates[3, :, :], 
#             showscale=false, 
#             name=nothing,
#         ),
#         plot_nodes(bb.coordinates; name=nothing, color="black", knot_size=1),
#     ], 
#     )
# )