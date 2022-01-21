module Splines
    using EllipsisNotation

    include("Utils.jl")
    include("Polynomials.jl")
    include("Visuals.jl")

    import .Polynomials: bernstein, N_basis
    import .Visuals: plot_curve, plot_nodes, plot_surface
    import .Utils: get_curve_parameters, ∑, dimensions, nnodes

    export plot_curve, plot_nodes, plot_surface 
    export Bezier, BSpline, AbstractSpline, AbstractBezier
    export ⊗, ∑


    # ---------------------------------------------------------------------------- #
    #                                Abstract types                                #
    # ---------------------------------------------------------------------------- #
    abstract type AbstractSpline end

    abstract type AbstractBezier <: AbstractSpline end
    abstract type AbstractBSpline <: AbstractSpline end

    # ---------------------------------- Methods --------------------------------- #

    function Base.show(io::IO, s::AbstractSpline)
        print(io, "Spline (nodes: $(s.η), N: $(s.N))")
    end

    """
        Spline(t, v, ...)

    Evaluates an N dimensional spline (with N > 1) at a set of N paramters t, v, ...
    values with each paramter `t ∈ [0, 1]`
    """
    function (bez::AbstractSpline)(params...)
        if length(params) != bez.N
            error("Number of parameters doesnt match the dimensionality of the spline surface")
        end

        sum((n)->bez.inner_fn[n](params[n]), 1:length(params))
    end


    # ---------------------------------------------------------------------------- #
    #                                   BSPLINES                                   #
    # ---------------------------------------------------------------------------- #


    struct BSpline <: AbstractBSpline
        nodes::AbstractArray            # coordinates of nodes (d x N (x M ...))
        coordinates::AbstractArray      # coordinates of points along the spline 
        weights::AbstractArray          # weights, shape N x M ...
        knots::AbstractArray
        degree::Int
        η                               # shape of the nodes -1  | degree (linera, quadratic, cubic..)
        d::Int                          # dimensionality of the Euclidean space
        N::Int                          # spline dimensionality: 1 - curve, 2 - surface...
        inner_fn                        # a collection of lower dimensionality bezier curves
    end

    # --------------------------- knots initialization --------------------------- #
    """
        uniform(n, d)

    Defines d+n+1 uniform knots ∈[0,1] for b-spline interpolation. 
    For i ∈ 0:n+d+1:
    if i ≤ d: k=0
    if d+1 ≤ i ≤ n: k=(i-d)/(n+1-d)
    else: k = 1
    See: https://www.geometrictools.com/Documentation/BSplineCurveLeastSquaresFit.pdf
    """
    function uniform(n::Int, d::Int)::AbstractArray
        if d > n
            @warn "When using B-splines, the degree should be < the number of control nodes"
        end

        knots = zeros(n+d+1+1)  # second 1 is because start from 0
        for i in 0:n+d+1
            if i ≤ d
                knots[i+1] = 0
            elseif (d+1) ≤ i ≤ n
                knots[i+1] = (i-d)/(n+1-d)
            else
                knots[i+1] = 1
            end
        end
        return knots
    end


    """
        periodic(n, d)

    Defines d+n+1 periodic knots ∈[0,1] for b-spline interpolation. 
    See: https://www.geometrictools.com/Documentation/BSplineCurveLeastSquaresFit.pdf
    """
    function periodic(n::Int, d::Int)::AbstractArray
        if d > n
            @warn "When using B-splines, the degree should be < the number of control nodes"
        end
        map((i)->(i-d)/(n+1-d), 0:n+d+1)
    end

    # ------------------------------- Constructors ------------------------------- #
    """
        BSpline(nodes::AbstractArray; knots_type::Symbol=:uniform, δt=0.025)::BSpline

    Constructor for a (non rational) B-Spline curve (N=1).
    """
    function BSpline(nodes::AbstractArray, degree::Int; knots_type::Symbol=:uniform, δt=0.025)::BSpline
        weights = ones(size(nodes)...)
        return BSpline(nodes, weights, degree;  knots_type= knots_type, δt=δt)
    end

    """
        BSpline(nodes::AbstractArray, weights::AbstractArray; knots_type::Symbol=:uniform δt=0.025)::BSpline

    Constructor for a rational B-Spline curve (N=1).
    """
    function BSpline(nodes::AbstractArray, weights::AbstractArray, degree::Int; knots_type::Symbol=:uniform, δt=0.025)::BSpline
        # get knots
        η = size(nodes, 2) - 1
        knots = eval(:($knots_type($η, $degree)))

        return BSpline(nodes, weights, knots, degree;  δt= δt)
    end

    """
        BSpline(nodes::AbstractArray, weights::AbstractArray, knots::AbstractArray; knots_type::Symbol=:uniform δt=0.025)::BSpline

    Constructor for a rational B Spline given a set of knot
    """
    function BSpline(nodes::AbstractArray, weights::AbstractArray, knots::AbstractArray, degree::Int; δt=0.025)::BSpline
        # get params
        τ, N, η, d = get_curve_parameters(nodes, δt)
        
        # compute coordinates
        numerator(i) = N_basis(τ; k=knots, i=i, d=degree)' .* nodes[:, i+1] .* weights[i+1]
        denominator(i) =  N_basis(τ; k=knots, i=i, d=degree)' .* weights[i+1]

        coordinates = ∑(numerator, 0:η) ./ ∑(denominator, 0:η)

        return BSpline(
            nodes, coordinates, weights, knots, degree, η, d, N, nothing
        )
    end

    
    """
        BSpline(t)

    Function to evaluate  a B Spline curve at a parameter `t ∈ [0, 1]`
    """
    function (bspline::AbstractBSpline)(t) 
        numerator(i) = N_basis(τ; k=bspline.knots, i=i, d=bspline.degree)' .* bspline.nodes[:, i+1] .* bspline.weights[i+1]
        denominator(i) =  N_basis(τ; k=bspline.knots, i=i, d=bspline.degree)' .* bspline.weights[i+1]

        ∑(numerator, 0:bspline.η) ./ ∑(denominator, 0:bspline.η)

    end
        

    # ---------------------------------------------------------------------------- #
    #                                    BEZIER                                    #
    # ---------------------------------------------------------------------------- #

    """
        Bezier

    Struct of type Bezier <: AbstractBezier.
    Contains all information about an N-dimensional Bezier surface.

    For `deg=1` we have a Bezier curve and:
        `nodes`: `dxN` array of nodes coordinates
        `η`: dimensionality of the nodes (`N-1`)
        `coordinates`: d x K points with the coordinates of the Bezier curve.

    for `deg>1` the type represents `deg` dimensional Bezier surfaces
    """
    struct Bezier <: AbstractBezier
        nodes::AbstractArray            # coordinates of nodes (d x N (x M ...))
        coordinates::AbstractArray      # coordinates of points along the spline 
        weights::AbstractArray          # weights, shape N x M ...
        η                               # shape of the nodes -1  | degree (linera, quadratic, cubic..)
        d::Int                          # dimensionality of the Euclidean space
        N::Int                          # spline dimensionality: 1 - curve, 2 - surface...
        inner_fn                        # a collection of lower dimensionality bezier curves
    end
    
    # ------------------------------- constructors ------------------------------- #
    """
        Bezier(nodes; δt=0.05)

    Constructor for a degree 1 (`n=1`) `Bezier` curve defined by a set of nodes
    """
    function Bezier(nodes::AbstractArray; δt=0.025)::Bezier
        weights = ones(size(nodes)...)
        return Bezier(nodes, weights; δt=δt)
    end

    """
        Bezier(nodes::AbstractArray, weights::AbstractArray; δt = 0.025)

        Constructor for a degree 1 rational bezier.
    """
    function Bezier(nodes::AbstractArray, weights::AbstractArray; δt = 0.025)::Bezier
        τ, N, η, d = get_curve_parameters(nodes, δt)

        numerator(i) = bernstein(τ; i=i, n=η)' .* nodes[:, i+1] .* weights[i+1]
        denominator(i) = bernstein(τ; i=i, n=η)' .* weights[i+1]
        coordinates = ∑(numerator, 0:η) ./ ∑(denominator, 0:η)

        return Bezier(
            nodes, coordinates, weights, η, d, N, nothing
        )
    end


    # ---------------------------- function evaluation --------------------------- #∂
    """
        Bezier(t)

    Evaluate a `Bezier` curve (N=1) at parameter value `t ∈ [0,1]`
    """
    function (bez::AbstractBezier)(t) 
        numerator(i) = bernstein(t; i=i, n=bez.η)' .* bez.nodes[:, i+1] .* bez.weights[i+1]
        denominator(i) = bernstein(t; i=i, n=bez.η)' .* bez.weights[i+1]

        ∑(numerator, 0:bez.η) ./ ∑(denominator, 0:bez.η)
    end



    # ---------------------------------------------------------------------------- #
    #                                Tensor Product                                #
    # ---------------------------------------------------------------------------- #

    """
        ×(x::AbstractArray, y::AbstractArray)

    Cartesian product between two arrays, without the 3 Euclidean coordinates.
    Given a N x M... array `x` and a P x Q array `y` it returns an array of shape
    N x M... x P x Q... with the product of the two arrays.
    """
    function  ×(x::AbstractArray, y::AbstractArray)::AbstractArray
        if ndims(x) == ndims(y) == 1
            return x .* y'
        else
            # initialize an array of the appropriate shape
            N = zeros(size(x)..., size(y)...)
            
            for Ix in CartesianIndices(x), Iy in CartesianIndices(y)
                N[Ix, Iy] = x[Ix] .+ y[Iy]
            end
        end
        return N
    end

    """
        ⊗(x::AbstractArray, y::AbstractArray)
    
    Tensor product between two arrays.
    Given x as a d x N x M... and y as a d x P x Q...
    it returns an array of shape d x N x M .. x P x Q x ...
    with the product of the two input arrays.

    This is done by iterating over the dimension d using Einstein summatin.
    Then each N x M... and P x Q... array is repeated/permuted to 
    to create two N x M... x P x Q... arrays which are then summed.
    """
    function ⊗(x::AbstractArray, y::AbstractArray)::AbstractArray
        if ndims(x) == ndims(y) == 1
            return x .* y'
        else
            # initialize an array of the appropriate shape
            N = zeros(3, dimensions(x) ..., dimensions(y)...)
            
            for Ix in CartesianIndices(x[1, ..]), Iy in CartesianIndices(y[1, ..])
                N[:, Ix, Iy] = x[:, Ix] .+ y[:, Iy]
            end
        end
        return N
    end

    """
        ⊗(b1::AbstractBezier, b2::AbstractBezier)::AbstractBezier

    Tensor product of two beziers
    """
    function ⊗(b1::T, b2::T) where T <: AbstractBezier
        # get nodes
        nodes = b1.nodes ⊗ b2.nodes
        n, m = nnodes(nodes)

        # get weights
        weights = b1.weights ⊗ b2.weights
        
        # get coordinates
        coordinates = b1.coordinates ⊗ b2.coordinates

        # concatenate inner functions
        β1 = b1.inner_fn == nothing ? [b1] : b1.inner_fn
        β2 = b2.inner_fn == nothing ? [b2] : b2.inner_fn
                
        return T(nodes, coordinates, weights, (n, m), b1.d, b1.N+b2.N, hcat(β1, β2))
    end

    """
        ⊗(b1::AbstractBSpline, b2::AbstractBSpline)::AbstractBSpline

    Tensor product of two beziers
    """
    function ⊗(b1::T, b2::T) where T <: AbstractBSpline
        if b1.degree != b2.degree
            @warn "BSpline ⊗ BSpline, '⊗' is only defined for B Splines with the same degree"
        end

        # get nodes
        nodes = b1.nodes ⊗ b2.nodes
        n, m = nnodes(nodes)

        # get knots
        knots = b1.knots × b2.knots

        # get weights
        weights = b1.weights ⊗ b2.weights
        
        # get coordinates
        coordinates = b1.coordinates ⊗ b2.coordinates

        # concatenate inner functions
        β1 = b1.inner_fn == nothing ? [b1] : b1.inner_fn
        β2 = b2.inner_fn == nothing ? [b2] : b2.inner_fn
                
        return T(nodes, coordinates, weights, knots, b1.degree, (n, m), b1.d, b1.N+b2.N, hcat(β1, β2))
    end


end