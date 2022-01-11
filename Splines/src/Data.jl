module TestData

    include("Geometry.jl")
    import .Geometry: asPoints

    export circle, sine, circle3D, circle4D

    noise(X; σ=0.1) = σ * randn(size(X))

    function circle(;δ=0.1, σ=.1)
        t =  δ:δ:2*pi
        X = [cos.(t) sin.(t)]
        return asPoints((X + noise(X; σ=σ))')
    end


    function sine(;δ=0.1, σ=.1)
        t =  δ:δ:2*pi
        X = [t sin.(t)]
        return asPoints((X + noise(X; σ=σ))')
    end


    function circle3D(;δ=.1, σ=.02)
        t =  δ:δ:2*pi
        X = [cos.(t) sin.(t) sin.(2t)]
        return asPoints((X + noise(X; σ=σ))')
    end


    function circle4D(;δ=.1, σ=.02)
        t =  δ:δ:2*pi
        X = [cos.(t) sin.(t) sin.(2t) cos.(2t)]
        return asPoints((X + noise(X; σ=σ))')
    end
end
