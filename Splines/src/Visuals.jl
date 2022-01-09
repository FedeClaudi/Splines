module Visuals
    using PlotlyJS

    export plot_fit_results
    
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
            scatter(x=curve[1, :], y=curve[2, :], z=curve[3, :], mode="lines", type="scatter3d", color="green", 
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