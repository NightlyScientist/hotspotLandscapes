using JLD2, Logging, JuliaSyntax

include("src/base/model.jl")
include("src/base/environment.jl")
include("src/visuals/snapshots.jl")
include("src/tools/trees.jl")
using .GraphModel

JuliaSyntax.enable_in_core!(true)

# command line options (input)
cli = parseArgs()

# set seed for reproducabilty
rngSeed = cli.rngSeed == 0 ? rand() : cli.rngSeed
GraphModel.StatsBase.Random.seed!(rngSeed)

# create containers for observales and associated functions
dataModels = GraphModel.predefinedModels()

# dataModel for animation
if cli.animation 
    dataModels[:animation] = DataModel{Tuple{Float32,Vector{Int32},Vector{Int32}}}(GraphModel.growAnimation)
end

# set landscape from commandLine params, optinally read in env
# bug: is :envType not being updated
env, htspts, cli = createEnvironment!(cli)

# send command line parameters to stdout and mkdirs
savePath = setPath(cli; cli.overwrite)
cli.printInfo && PrintSimInfo(cli)

if cli.debug
    global_logger(ConsoleLogger(stderr, Logging.Debug))
end

# build hex map
graph, cntns = hexGraph(cli.dims)

# save all data in files
jldsave(savePath * "/data.jld2";
    env = env, cli = cli, htspts = htspts, rngSeed = rngSeed
)

function run(trial, dataModels, cli, env, htspts)
    # reset models
    GraphModel.resetModels!(dataModels)

    # reset graph properties
    resetGraph!(graph, cntns)

    # initial population
    active = linearStart!(graph, cntns, cli.dims, env)

    # do simulation
    final = simulate!(graph, cntns, active, env, cli, dataModels)

    # sources
    width, height = cli.dims
    sources = collect(1:width) .+ width*(height - 1)
    
    # save rough-front stopping condition sources
    cli.roughfront && (sources = final.fSnapshot)

    # determine the phylo of the periphery population 
    phylo, branches = GraphModel.genealogy(graph, sources)

    # save trial data
    jldopen(savePath * "/data.jld2", "a+") do file
        ng = JLD2.Group(file, "trial_$trial")
        for (k, v) in GraphModel.unpack(dataModels)
            ng[String(k)] = v
        end
        ng["finalTime"] = final.time
        ng["sources"] = sources
        ng["pIDs"] = getfield.(graph[sources], :ID_2)
        ng["phylo"] = phylo
        ng["tree"] = GraphModel.buildTree(graph, branches, sources)
    end

    # snapshot of sectors, lineages, hotspots
    if trial == 1 && cli.sampleImages
        fig, _ = landVisual(cli.dims, env, getfield.(graph, :ID_2), phylo; drawScatter=true)
        save(savePath * "/sampleImage.png", fig)


        # fig, _ = landVisual(cli.dims, env, getfield.(graph, :ID_2), phylo; drawScatter=true)
        # save(savePath * "/sampleImage_2.png", fig)
    end

    # save full graph for testing
    if trial == 1 && cli.debug
        jldopen(savePath * "/fullData.jld2", "a+") do file
            file["graph"] = graph
            file["sources"] = final.fSnapshot
        end
    end
end

# run all trials
for trial in range(1, cli.numberTrials)
    run(trial, dataModels, cli, env, htspts)
end
