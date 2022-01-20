module Beziers
    include("Spline.jl")
    include("Polynomials.jl")
    include("Maths.jl")

    
    import .Spline: AbstractBezier
    import .Polynomials: bernstein
    import .Maths: ∑, nnodes
    import .Maths: ⊗

    export Bezier

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
        nodes::AbstractArray  # coordinates of nodes (d x N (x M ...))
        coordinates::AbstractArray   # coordinates of points along the spline 
        η  # shape of the nodes -1  | degree (linera, quadratic, cubic..)
        d::Int  # dimensionality of the Euclidean space
        N::Int  # spline dimensionality: 1 - curve, 2 - surface...
        inner_fn  # a collection of lower dimensionality bezier curves
    end
    

    """
        Bezier(nodes; δt=0.05)

    Constructor for a degree 1 (`n=1`) `Bezier` curve defined by a set of nodes
    """
    function Bezier(nodes; δt=0.025)::Bezier
        if ndims(nodes) > 2
            @warn "Should only be invoked for deg=1 Bezier curves"
        end

        # parameter interval
        τ = 0:δt:1

        # dimesions
        N = 1
        η = nnodes(nodes)
        d = size(nodes, 1)

        # compute coordinates
        B(i) = bernstein(τ; i=i, n=η)' .* nodes[:, i+1]
        coordinates = ∑(B, 0:η)
    
        return Bezier(nodes, coordinates, η, d, N, nothing)
    end

    """
        bezier(t)

    Turn an instance of `Bezier` into a callable.
    Evaluates a bezier curve (N=1) at parameter value `t ∈ [0,1]`
    """
    function (bez::Bezier)(t)        
        B(i) = bernstein(t; i=i, n=bez.η)' .* bez.nodes[:, i+1]
        ∑(B, 0:bez.η)
    end

    """
    Method for N>1 bezier surfaces. 
    It evaluates each N=1 curve making the surface up at the corresponding parameter value.
    """
    function (bez::Bezier)(params...)        
        if length(params) != bez.N
            error("Number of parameters doesnt match the dimensionality of the Bezier curve")
        end
        

        ∑((n)->bez.inner_fn[n](params[n]), 1:length(params))
    end


    """
        b1 ⊗ b2

    Tensor product of two beziers
    """
    function ⊗(b1::Bezier, b2::Bezier)::Bezier
        nodes = b1.nodes ⊗ b2.nodes
        n, m = nnodes(nodes)

        coordinates = b1.coordinates ⊗ b2.coordinates

        # concatenate inner functions
        β1 = b1.inner_fn == nothing ? [b1] : b1.inner_fn
        β2 = b2.inner_fn == nothing ? [b2] : b2.inner_fn
                
        return Bezier(nodes, coordinates, (n, m), b1.d, b1.N+b2.N, hcat(β1, β2))
    end


end