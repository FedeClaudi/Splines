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
    import .Interpolation: PiecewiseLinear, PiecewiseLinear!, BSpline, BSpline!, Bezier, Bezier!
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
    function cost(
            nodes::Points,
            labelled_data::Vector{Matrix{Float64}},  # datapoints with cluster label assigned
            curve::Points,
            curve_fn!;
            α=1.0,
            β=1.0,
            kwargs...
        )::Float64

        # compute length of curve
        length_cost = gm.curve_length(curve_fn!(curve, nodes; kwargs...).points)

        # compute distance between knot position and datapoints assigned to it by clustering
        distance_cost = 0.0
        for n in range(1, stop=size(nodes, 2))
            distance_cost += sum(gm.distances(labelled_data[n], nodes[:, n]))
        end

        return α * length_cost + β * distance_cost
    end

    """
        Fits a data array (d x N) with an interpolated curve specified by a set of 
        n nodes. The position of the nodes is intialized with `kmeans` as the centroids of n-many clusters.
        The position is optimized to minimize a cost function of curve length (scaled by `α`) and nodes-data 
        distance (scaled by `β`).

        The interpolated curve can be:
            - PiecewiseLinear (`curve_fn=:PiecewiseLinear`)
            - B-spline (`curve_fn=:BSpline`)

        `kwargs...` gets passed to the curve function.
    """
    function fit(
            data::Points, 
            curve_fn::Symbol; 
            n::Int=12, 
            α::Float64=1.0, 
            β::Float64=1.0, 
            nodes_sorting_method::Symbol=:smallest,
            n_iter::Int=250,
            kwargs...
        )
        @info "Fitting $curve_fn | $(size(data, 2)) data points with $n nodes - α=$α, β=$β"
        @debug "Keyword arguments: $(["$(v[1]):$(v[2])" for v in kwargs])"

        # initialize KNOTS to data
        clusters = kmeans(data, n)
        nodes_init = sort_points(clusters.centers; selection_method=nodes_sorting_method)

        labels = assignments(kmeans!(data, nodes_init))  # re-run because now nodes are sorted
        labelled_data = [data[:, findall(labels .== n)] for n in range(1, stop=n)]

        # get callables for curve generating functions
        curve_fn! = eval(Symbol(curve_fn, "!"))  # from :Symbol to in-place version of fn
        curve_fn = eval(curve_fn)  # from ::Symbol to callable function

        # initialize a curve array
        curve = curve_fn(nodes_init; kwargs...).points

        # optimize nodes position
        @debug "Optimizing nodes placement"
        𝐿(k) = cost(k, labelled_data, curve, curve_fn!; α=α, β=β, kwargs...)
        opt_res = optimize(𝐿, nodes_init, iterations=n_iter)
        @debug opt_res
        nodes_optim = sort_points(opt_res.minimizer, selection_method=nodes_sorting_method)

        # create curve
        @debug "Creating curve from nodes"
        curve = curve_fn(nodes_optim; kwargs...)
        
        @info "Curve fitting complete"
        return nodes_init, nodes_optim, curve, opt_res
    end



end