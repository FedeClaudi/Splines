module Splines
    include("Types.jl")
    include("Geometry.jl")
    include("Data.jl")
    include("Curves.jl")
    include("Fit.jl")
    include("Visuals.jl")
    include("Polynomials.jl")
    include("Surfaces.jl")
    include("Utils.jl")

    import .Types: Point, Points, Curve
    import .Geometry: distances, curve_length
    import .Curves: 
        bspline, bspline!,
        bezier, bezier!,
        rational_bezier, rational_bezier!,
        rational_bspline, rational_bspline!
    import .TestData
    import .Fit: fit
    import .Visuals: plot_fit_results, plot_surface
    import .Polynomials
    import .Surfaces: bezier_surface, bspline_surface
    import .Utils

    export Point, Points, bspline, bspline!, bezier, bezier!, fit, plot_fit_results, bezier_surface, NURBS_curve, NURBS_curve!

end 
