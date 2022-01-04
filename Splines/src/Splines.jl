module Splines
    include("Geometry.jl")
    include("Interpolation.jl")

    import .Geometry: Point, Points, distances, curve_length
    import .Interpolation: PiecewiseLinear, lerp

    export Point, Points, distances, curve_length, PiecewiseLinear, lerp

end 
