"""
Module for fitting splines/interpolated curves to data
"""
module Fit
    import Clustering: kmeans, assignments, kmeans!
    import Optim: optimize

    include("Geometry.jl")
    include("Interpolation.jl")
    include("Utils.jl")

    import .Geometry: Point, Points
    import .Geometry as gm
    import .Interpolation: PiecewiseLinear, PiecewiseLinear!
    import .Utils: sort_points

    export fitPWL

    """
        Computes the cost given a set of nodes and data points.
        The loss is proportional to the length of the curve (generated
        by interpolation through the nodes) and the distance to all data
        points.

        # Arguments
            - α : weight for the length cost
            - β : weight for the distance error cost
    """
    function cost(nodes::Points, labelled_data::Vector{Matrix{Float64}}, curve::Points; α=1.0, β=1.0, closed::Bool=false)::Float64
        # compute length of curve
        length_cost = gm.curve_length(PiecewiseLinear!(curve, nodes; closed=closed))

        # compute distance between knot position and datapoints assigned to it by clustering
        distance_cost = 0.0
        for n in range(1, stop=size(nodes, 2))
            distance_cost += sum(gm.distances(labelled_data[n], nodes[:, n]))
        end

        return α * length_cost + β * distance_cost
    end

    """
        Fits a data array (d x N) with a PiecewiseLinear curve specified by a set of 
        η nodes. The position of the nodes is intialized with `kmeans` and optimized to 
        minimize a cost that is a function of curve length (scaled by `α`) and nodes-data 
        distance (scaled by `β`)
    """
    function fitPWL(data; η=12, α=1.0, β=1.0, closed::Bool=false, nodes_sorting_method::Symbol=:smallest)
        @info "Fitting $(size(data, 2)) data points with $η nodes - α=$α, β=$β"
        # initialize KNOTS to data
        clusters = kmeans(data, η)
        knots_init = sort_points(clusters.centers; selection_method=nodes_sorting_method)

        labels = assignments(kmeans!(data, knots_init))  # re-run because now nodes are sorted
        labelled_data = [data[:, findall(labels .== n)] for n in range(1, stop=η)]

        # initialize a curve array
        curve = PiecewiseLinear(knots_init; closed=closed)

        # optimize nodes position
        @debug "Optimizing nodes placement"
        𝐿(k) = cost(k, labelled_data, curve; α=α, β=β, closed=closed)
        knots_optim = sort_points(optimize(𝐿, knots_init, iterations=100, show_every=3).minimizer, selection_method=nodes_sorting_method)

        # create curve
        @debug "Fitting Piecewise linear to nodes"
        curve = PiecewiseLinear(knots_optim; closed=closed)

        return knots_init, knots_optim, curve
    end



end