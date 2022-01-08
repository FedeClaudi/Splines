module Splines
    include("Geometry.jl")
    include("Data.jl")
    include("Interpolation.jl")
    include("Fit.jl")

    import .Geometry: Point, Points, distances, curve_length
    import .Interpolation: PiecewiseLinear, PiecewiseLinear!, lerp, BSpline
    import .TestData
    import .Fit

    export Point, Points, distances, curve_length, PiecewiseLinear, PiecewiseLinear!, lerp, BSpline

end 
