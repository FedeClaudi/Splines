module Splines
    include("Types.jl")
    include("Geometry.jl")
    include("Data.jl")
    include("Curves.jl")
    include("Fit.jl")
    include("Visuals.jl")
    include("Polynomials.jl")

    import .Types: Point, Points
    import .Geometry: distances, curve_length
    import .Curves: 
        BSpline, BSpline!,
        Bezier, Bezier!,
        RationalBezier, RationalBezier!
    import .TestData
    import .Fit: fit
    import .Visuals: plot_fit_results
    import .Polynomials

    export Point, Points, BSpline, BSpline!, Bezier, Bezier!, fit, plot_fit_results

end 
