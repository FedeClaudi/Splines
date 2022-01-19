module Spline


# ---------------------------------------------------------------------------- #
#                                Abstract types                                #
# ---------------------------------------------------------------------------- #
abstract type AbtractSpline 
end

function Base.show(io::IO, s::AbtractSpline)
    print(io, "Spline (nodes: $(s.η), N: $(s.N))")
end

abstract type AbstractNURBS <: AbtractSpline end

abstract type AbstractBspline <: AbstractNURBS end
abstract type AbstractBezier <: AbstractBspline end

abstract type AbstractRationalBspline <: AbstractNURBS end
abstract type AbstractRationalBezier <: AbstractRationalBspline end


end