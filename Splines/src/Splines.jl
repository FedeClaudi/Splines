module Splines
    include("Spline.jl")
    include("Polynomials.jl")
    include("Utils.jl")

    include("Beziers.jl")
    include("Maths.jl")
    include("Visuals.jl")

    import .Beziers: Bezier, RationalBezier
    import .Maths: ⊗, dimensions
    import .Polynomials
    import .Spline: AbstractSpline, AbstractBezier
    import .Visuals: plot_curve, plot_nodes, plot_surface
    import .Utils: ∑

    export plot_curve, plot_nodes, plot_surface 
    export ∑, ⊗, dimensions
    export Bezier, RationalBezier
    export AbstractSpline, AbstractBezier

end