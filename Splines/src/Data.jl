module TestData

    include("Geometry.jl")
    import .Geometry: Points

    export circle, sine, circle3D

    noise(X; σ=0.1) = σ * randn(size(X))

    function circle(;σ=.1)
        t =  .01:.01:2*pi
        X = [cos.(t) sin.(t)]
        return Points((X + noise(X; σ=σ))')
    end


    function sine(;σ=.1)
        t =  .01:.01:2*pi
        X = [t sin.(t)]
        return Points((X + noise(X; σ=σ))')
    end


    function circle3D(;σ=.02)
        t =  .01:.01:2*pi
        X = [cos.(t) sin.(t) sin.(2t)]
        return Points((X + noise(X; σ=σ))')
    end

end
