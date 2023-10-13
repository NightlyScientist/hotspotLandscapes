using ArgParse

function PrintSimInfo(parsedArgs)
    println("Simulation Info:")
    for k in keys(parsedArgs)
        println("  > $k  =>  $(parsedArgs[k])")
    end
end

addArgs!(sts, name; kwargs...) = add_arg_table!(sts, "--$(name)", Dict(kwargs))

function ArgParse.parse_item(::Type{NTuple{2,T}}, x::AbstractString) where T
    return Tuple(convert.(T, parse.(Float64, split(x, ','))))
end

function parseArgs()
    sts = ArgParseSettings()
    sts.prog = "Eden Model of bacterial growth in noisy environments"

    addArgs!(sts, "numberTrials", arg_type=Int, default=1)
    addArgs!(sts, "numberSamples", arg_type=Int, default=50)
    addArgs!(sts, "dims", required=true, arg_type=NTuple{2,Int64})
    addArgs!(sts, "density", required=true, arg_type=Float64)
    addArgs!(sts, "radius", required=true, arg_type=Int64)
    addArgs!(sts, "intensity", required=true, arg_type=Float64)
    addArgs!(sts, "hsRegion", arg_type=NTuple{2,Int64})
    addArgs!(sts, "rngSeed", arg_type=Int, default=1)
    addArgs!(sts, "printInfo", action=:store_true)
    addArgs!(sts, "overwrite", action=:store_true)
    addArgs!(sts, "savePath", required=true, arg_type=String)
    addArgs!(sts, "landscape", arg_type=String)
    addArgs!(sts, "sampleImages", action=:store_true)
    addArgs!(sts, "debug", action=:store_true)
    addArgs!(sts, "roughfront", action=:store_true)
    addArgs!(sts, "animation", action=:store_true)
    addArgs!(sts, "parameter", arg_type=String, default="missing")
    addArgs!(sts, "envType", default="uniform", arg_type=String)
    namedtuple(parse_args(sts))
end
