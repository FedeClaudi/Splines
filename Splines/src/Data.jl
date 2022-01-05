module TestData

    include("Geometry.jl")
    import .Geometry: Points

    export circle, sine

    noise(X; σ=0.1) = σ * randn(size(X))

    function circle()
        t =  .01:.01:2*pi
        X = [cos.(t) sin.(t)]
        return Points((X + noise(X))')
    end


    function sine()
        t =  .01:.01:2*pi
        X = [t sin.(t)]
        return Points((X + noise(X))')
    end



end
