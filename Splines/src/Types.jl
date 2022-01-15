module Types
    export Point, Points, asPoints, Knots, Curve, Surface

    const Point{T<:Real} = Vector{T}
    const Points{T<:Real} = AbstractArray{T}
    asPoints(X) = AbstractArray{Float64}(X)

    const Knots = Vector{Float64}

    # ---------------------------------------------------------------------------- #
    #                                     CURVE                                    #
    # ---------------------------------------------------------------------------- #
    """
        Curve

    Stores the results of fitting a spline to a set of nodes.
    Contains the parameter values (`τ`) at which the fitted spline
    was evaluated at, the coordinates (`points`) of the spline and a 
    function (`f`) with signature `f(t::Number)::Vector{Number}` giving
    the position of the spline at a parameter value `t ∈ [0,1]`.
    """
    struct Curve
        name::String
        nodes::Points        # the curve control nodes
        τ::Vector{Float64}   # parameter values used when constructing curve
        points::Points       # points along curve at param values
        func                 # function to evaluate the curve at a param value
    end

    """
        Make `Curve` callable such that `mycurve = Curve(...)` and then `mycurv(t)`
        results in evaluating the curve function at the parameter value `t`.
    """
    function (curve::Curve)(t)
        curve.func(t)
    end

    """
        Custom printing of curve object
    """
    function Base.show(io::IO, curve::Curve)
        print(io, "Curve: '$(curve.name)' $(size(curve.points))")
    end


    # ---------------------------------------------------------------------------- #
    #                                    SURFACE                                   #
    # ---------------------------------------------------------------------------- #

    """
        Surface

    Stores data produced from computing a spline given a set of knots.
    """
    struct Surface
        name:: String
        nodes::Points
        params::Matrix{Float64}  # τ .* τ'
        points::Array  # 3 x T x T
    end

    """
    Custom printing of Surface object
    """
    function Base.show(io::IO, surface::Surface)
        print(io, "Surface: '$(surface.name)' $(size(surface.points))")
    end
end