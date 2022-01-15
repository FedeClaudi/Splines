module Splines
    include("Types.jl")
    include("Geometry.jl")
    include("Data.jl")
    include("Curves.jl")
    include("Fit.jl")
    include("Visuals.jl")
    include("Polynomials.jl")
    include("Surfaces.jl")

    import .Types: Point, Points, Curve
    import .Geometry: distances, curve_length
    import .Curves: 
        BSpline, BSpline!,
        Bezier, Bezier!,
        RationalBezier, RationalBezier!
    import .TestData
    import .Fit: fit
    import .Visuals: plot_fit_results, plot_surface
    import .Polynomials
    import .Surfaces: BezierSurface

    export Point, Points, BSpline, BSpline!, Bezier, Bezier!, fit, plot_fit_results, BezierSurface

end 
