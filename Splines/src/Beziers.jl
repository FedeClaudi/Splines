module Beziers
    include("Spline.jl")
    include("Polynomials.jl")
    include("Maths.jl")

    import StaticArrays: @SArray

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
    
        return Bezier(nodes, coordinates, η, d, N)
    end

    function Base.show(io::IO, b::Bezier)
        print(io, "Bezier (nodes: $(b.η), N: $(b.N))")
    end


    function ⊗(b1::Bezier, b2::Bezier)::Bezier
        nodes = b1.nodes ⊗ b2.nodes
        n, m = nnodes(nodes)
        coordinates = b1.coordinates ⊗ b2.coordinates

        return Bezier(nodes, coordinates, (n, m), b1.d, b1.N+b2.N)
    end


end