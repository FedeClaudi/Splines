"""
Module for fitting splines/interpolated curves to data
"""
module Fit
    import Clustering: kmeans, assignments, kmeans!
    import Optim: optimize

    include("Geometry.jl")
    include("Interpolation.jl")

    import .Geometry: Point, Points
    import .Geometry as gm
    import .Interpolation: PiecewiseLinear, PiecewiseLinear!

    export fitPWL

    """
        Computes the cost given a set of knots and data points.
        The loss is proportional to the length of the curve (generated
        by interpolation through the knots) and the distance to all data
        points.

        # Arguments
            - α : weight for the length cost
            - β : weight for the distance error cost
    """
    function cost(knots::Points, data::Points, labels::Vector{Int64}, curve::Points; α=1.0, β=1.0, closed::Bool=false)::Float64
        # compute length of curve
        length_cost = gm.curve_length(PiecewiseLinear!(curve, knots; closed=closed))

        distance_cost = 0.0
        for n in range(1, stop=size(knots, 2))
            # compute distance between knot position and datapoints assigned to it by clustering
            distance_cost += sum(gm.distances(data[:, findall(labels .== n)], knots[:, n]))
        end

        return α * length_cost + β * distance_cost
    end

    """
        Fits a data array (d x N) with a PiecewiseLinear curve specified by a set of 
        η knots. The position of the knots is intialized with `kmeans` and optimized to 
        minimize a cost that is a function of curve length (scaled by `α`) and knots-data 
        distance (scaled by `β`)
    """
    function fitPWL(data; η=12, α=1.0, β=1.0, closed::Bool=false)
        # initialize KNOTS to data
        clusters = kmeans(data, η)
        knots_init = gm.sort_points(clusters.centers; selection_method=:smallest)
        labels = assignments(kmeans!(data, knots_init))  # re-run because now knots are sorted

        # initialize a curve array
        curve = PiecewiseLinear(knots_init; closed=closed)

        # optimize knots position
        @info "Optimizing knots placement"
        𝐿(k) = cost(k, data, labels, curve; α=α, β=β, closed=closed)
        knots_optim = gm.sort_points(optimize(𝐿, knots_init).minimizer, selection_method=:smallest)

        # create curve
        @info "Fitting Piecewise linear to knots"
        curve = PiecewiseLinear(knots_optim; closed=closed)

        return knots_init, knots_optim, curve
    end



end