module Visuals
    using PlotlyJS

    include("Utils.jl")
    include("Types.jl")

    import .Utils: x, y, z
    using .Types


    export plot_fit_results, plot_curve, plot_nodes, plot_surface
    
    """
        plot_curve(
            curve::Curve; 
            as_scatter::Bool=false, 
            color::String="black",
            width::Number=10,
            )

    Creates a scatter or line representation of a 1 
    dimensional curve in N dimensional space
    """
    function plot_curve(
            points::AbstractArray;
            as_scatter::Bool=false, 
            color::String="black",
            width::Number=10,
            name::String="curve"
            )

            # prep params
            N = size(curve.points, 1)
            if N == 3
                _type = "scatter3d"
            else
                _type = "scatter"
            end

            if as_scatter
                mode= "marker"
            else
                mode = "lines"
            end

            # plot curve
            scatter(
                x=x(points), 
                y=y(points), 
                z=z(points), 
                mode=mode, 
                type=_type,
                line = attr(color=color, width=width), 
                name=name
            )
    end



    function plot_nodes(
            nodes::AbstractArray; 
            color::String="black",
            name="nodes",
            knot_size::Number=8,
        )
        N = ndims(nodes)
        if size(nodes, 1) == 3
            _type = "scatter3d"
        else
            _type = "scatter"
        end

        scatter(
            x= (N==2) ? x(nodes) : vec(x(nodes)), 
            y= (N==2) ? y(nodes) : vec(y(nodes)), 
            z= (N==2) ? z(nodes) : vec(z(nodes)), 
            mode="markers", 
            type=_type,     
            marker=attr(
                    size=knot_size,
                    color=color,
                    opacity=1,
                ), 
            name=name
            )
    end


    function plot_surface(surf; with_points::Bool=true)
        srf = surface(
            x=x(surf.points), 
            y=y(surf.points), 
            z=z(surf.points), 
            name=surf.name,
            showscale=false
            )

        if with_points == false
            return srf
        end

        sct = scatter(
            x=vec(x(surf.points)), 
            y=vec(y(surf.points)), 
            z=vec(z(surf.points)), 
            mode="markers", 
            type="scatter3d",
            marker=attr(
                    size=1,
                    color="black",
                ),
            name=surf.name,
            )

        return [srf, sct]


    end



    function plot_fit_results(data, curve, nodes_optim, nodes_init)
        display(plot(
            [
            # plot data
            scatter(x=data[1, :], y=data[2, :], z=data[3, :], mode="markers", type="scatter3d",     
                marker=attr(
                    size=3,
                    color="grey",                # set color to an array/list of desired values
                    edgecolor="black",
                    opacity=.2,
                ), name="data"
            ),
        
            # plot curve
            scatter(x=curve.points[1, :], y=curve.points[2, :], z=curve.points[3, :], mode="lines", type="scatter3d", color="green", 
                    line = attr(color="royalblue", width=4), name="curve"),
        
            # plot nodes
            scatter(x=nodes_optim[1, :], y=nodes_optim[2, :], z=nodes_optim[3, :], mode="markers", type="scatter3d",     
                marker=attr(
                    size=6,
                    color="red",                # set color to an array/list of desired values
                    edgecolor="black",
                    opacity=.6,
                ), name="optim. nodes"
            ),
            scatter(x=nodes_init[1, :], y=nodes_init[2, :], z=nodes_init[3, :], mode="markers", type="scatter3d",     
                marker=attr(
                    size=4,
                    color="black",                # set color to an array/list of desired values
                    edgecolor="black",
                    opacity=1,
                ), name="init. nodes"
            ),
        
            ]
        ))
    end
end