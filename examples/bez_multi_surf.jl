# ? activate Splines
import Pkg
Pkg.activate("Splines")

using PlotlyJS
using Revise
Revise.revise()

using Splines

@info "Computing"
n, m = 8, 8


# cupula
N1 = [
        cat([0:n-1 zeros(n) 4 .* sin.(range(0, stop=π, length=n))], dims=2)',
        cat([zeros(m) 0:m-1 4 .* sin.(range(0, stop=π, length=n))], dims=2)'
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
    # compute beziers
    b1 = Bezier(nodes_1)
    b2 = Bezier(nodes_2)
    bb = b1 ⊗ b2

    # compute and plot points on the bezier curve
    points = hcat(b1(.5), b2(.5), bb(.5, .5))
    add_trace!(fig,
    plot_nodes(points, color="black", node_size=8), row=rows[n], col=cols[n]
    )

    # plot knots
    for (nodes, color) in zip((b1.nodes, b2.nodes, bb.nodes), ("red", "blue", "black"))
        add_trace!(fig,
            plot_nodes(nodes; color=color, node_size=1), row=rows[n], col=cols[n]
        )
    end

    # plot curves
    for (b, color) in zip((b1, b2), ("red", "blue"))
        add_trace!(fig,
            plot_curve(b.coordinates, color=color), row=rows[n], col=cols[n])
    end

    # plot surface
    add_trace!(
        fig,
        plot_surface(bb.coordinates, opacity=1)[1],
        row=rows[n], col=cols[n]
    )


end

display(fig)
@info "Done"