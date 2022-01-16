module Spline


# ---------------------------------------------------------------------------- #
#                                Abstract types                                #
# ---------------------------------------------------------------------------- #
abstract type AbtractSpline 
end

abstract type AbstractNURBS <: AbtractSpline end

abstract type AbstractBspline <: AbstractNURBS end
abstract type AbstractBezier <: AbstractBspline end

abstract type AbstractRationalBspline <: AbstractNURBS end
abstract type AbstractRationalBezier <: AbstractRationalBspline end


end