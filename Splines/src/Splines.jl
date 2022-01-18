module Splines
    include("Beziers.jl")
    include("Maths.jl")
    include("Polynomials.jl")
    include("Spline.jl")
    include("Visuals.jl")

    import .Beziers: Bezier, ⊗
    import .Maths: ⊗, dimensions
    import .Polynomials
    import .Spline
    import .Visuals: plot_curve, plot_nodes, plot_surface

    export plot_curve, plot_nodes, plot_surface 
    export ∑, ⊗, dimensions
    export Bezier

end