module Spline

    export AbstractBezier

    # ---------------------------------------------------------------------------- #
    #                                Abstract types                                #
    # ---------------------------------------------------------------------------- #
    abstract type AbstractSpline end

    function Base.show(io::IO, s::AbstractSpline)
        print(io, "Spline (nodes: $(s.η), N: $(s.N))")
    end

    """
    Method for N>1 bezier surfaces. 
    It evaluates each N=1 curve making the surface up at the corresponding parameter value.
    """
    function (bez :: AbstractSpline)(params...)        
        if length(params) != bez.N
            error("Number of parameters doesnt match the dimensionality of the Bezier curve")
        end

        ∑((n)->bez.inner_fn[n](params[n]), 1:length(params))
    end


    abstract type AbstractBezier <: AbstractSpline end



    # abstract type AbstractNURBS <: AbstractSpline end

    # abstract type AbstractBezier <: AbstractBspline end

    # abstract type AbstractBspline <: AbstractNURBS end
    # abstract type AbstractBezier <: AbstractBspline end

    # abstract type AbstractRationalBspline <: AbstractNURBS end
    # abstract type AbstractRationalBezier <: AbstractRationalBspline end


end