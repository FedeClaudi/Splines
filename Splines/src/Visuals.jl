module Visuals

    using PlotlyJS

    import ..Maths: x, y, z

    export plot_nodes, plot_curve, plot_surface
    
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
            points;
            as_scatter::Bool=false, 
            color::String="black",
            width::Number=10,
            name::String="curve"
            )

            # prep params
            N = size(points, 1)
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

    function plot_surface(points; with_points::Bool=true, name::String="surface")
        srf = surface(
            x=x(points), 
            y=y(points), 
            z=z(points), 
            name=name,
            showscale=false
            )

        if with_points == false
            return srf
        end

        sct = scatter(
            x=vec(x(points)), 
            y=vec(y(points)), 
            z=vec(z(points)), 
            mode="markers", 
            type="scatter3d",
            marker=attr(
                    size=1,
                    color="black",
                ),
            name=name,
            )

        return [srf, sct]


    end
end