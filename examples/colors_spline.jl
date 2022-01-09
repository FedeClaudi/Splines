using PlotlyJS


data = range(-8, stop=8, length=40)

X, Y, Z = mgrid(data, data, data)


values = sin.(X .* Y .* Z) ./ (X .* Y .* Z)


display(plot(volume(

    x=X[:],

    y=Y[:],

    z=Z[:],

    value=values[:],

    isomin=0.1,

    isomax=0.8,

    opacity=0.1, # needs to be small to see through all surfaces

    surface_count=17, # needs to be a large number for good volume rendering

)))