module Splines
    include("Geometry.jl")
    include("Data.jl")
    include("Interpolation.jl")
    include("Fit.jl")
    include("Visuals.jl")

    import .Geometry: Point, Points, distances, curve_length
    import .Interpolation: PiecewiseLinear, PiecewiseLinear!, lerp, BSpline
    import .TestData
    import .Fit: fit
    import .Visuals: plot_fit_results

    export Point, Points, distances, curve_length, PiecewiseLinear, PiecewiseLinear!, lerp, BSpline, fit, plot_fit_results

end 
