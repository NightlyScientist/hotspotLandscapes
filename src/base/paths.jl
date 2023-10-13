function setPath(cli; overwrite::Bool=false)::String
    mainDir = cli.savePath
    density = "D-$(round(cli.density, digits=4))"
    intensity = "I-$(round(cli.intensity, digits=4))"
    radius = "R-$(cli.radius)"
    dimensions = "XY-$(cli.dims)"
    hotspotRegion = "HR-$(cli.hsRegion)"
    trials = "trials-" * string(cli.numberTrials)
    samples = "samples-" * string(cli.numberSamples)

    path = mainDir * "/" * join([density, intensity, radius, dimensions, hotspotRegion, trials, samples], ",")

    # strings = []
    # for (k,v) in zip(keys(a), values(a))
    #     val = v
    #     if typeof(v) <: Float64
    #         val = round(v, digits=4)
    #     end
    #     push!()
    # end
    setPath(path, overwrite, cli)
    return path
end

function setPath(path::String, overwrite::Bool, cli)::Nothing
    if ispath(path)
        overwrite && rm(path, recursive=true)
        overwrite || error(" ! save path not empty ! | $cli ")
    end
    mkpath(path)
    return nothing
end