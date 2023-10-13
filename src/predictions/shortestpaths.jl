using JLD2, FileIO, CodecZlib
using Base.Threads
using ArgParse

include("./Dijkstra.jl")
include("../base/model.jl")
include("../base/environment.jl")
include("../tools//containers.jl")
include("../tools/parseFiles.jl")

import .GraphModel: gNodes, buildMap
import .WeightedGraphs: sPaths, Node

cli = begin 
    sts = ArgParseSettings()
    add_arg_table!(sts, ["--input"], Dict(:required=>true))
    add_arg_table!(sts, ["--output"], Dict(:required=>true))

    namedtuple(parse_args(sts))
end

function shortestpaths(fpath)
    p, env = load(fpath, "cli", "env")

    # dimensions
    width, height = p.dims

    # sinks and sources
    srcs = collect(1:width) .+ width*(height - 1)
    tgts = collect(1:width)

    # rebuild rates
    rates = Float64[1., 1 + p.intensity]

    # build graph and connections
    graph, cntns = buildMap(p.dims, gNodes(), Node, (n,x,y)->Node())

    # build weights from rates
    envT = map(t->1/rates[t], env)

    savePath = let pname = rstrip(fpath, '/')
        tmp = split(pname, '/')
        joinpath(cli.output, join(tmp[end-3:end-1], "/"))
    end

    # perform calculation
    if !isfile(joinpath(savePath, "spaths.jld2"))

        # set path, and NEVER overwrite
        setPath(savePath, false, nothing)

        weights, tracks = sPaths(graph, cntns, srcs, tgts, envT, 25) 

        # save data of shortests paths
        jldsave(joinpath(savePath , "spaths.jld2"); 
            tracks=tracks, weights=weights
        )
    end

    # symbolic link to spaths in predictions/ directory
    if isfile(joinpath(savePath, "spaths.jld2"))

        rp = readlink("./storage")

      
        # hack: check if first path is not already /storage/
        if contains(split(fpath, '/')[2], "jgonzaleznunez")
          rp = "storage" 
        end

        # remove symlink and replace with abs path
        spath = replace(first(splitdir(fpath)), "storage"=>rp) * "/spaths.jld2"

        tgtpath = replace(savePath, "storage"=>rp) * "/spaths.jld2"

        # remove old symlink
        islink(spath) && rm(spath)

        @info tgtpath spath
        symlink(tgtpath, spath)
    end
    return nothing
end

for (root, dirs, files) in walkdir(cli.input)
    # process only "collection-n" directories
    contains(lowercase(root), "collection-") || continue

    @threads for i in eachindex(files)
        if contains(joinpath(root, files[i]), "data.jld2") 
            shortestpaths(joinpath(root, files[i]))
        end
    end
end
