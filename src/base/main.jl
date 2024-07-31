using JLD2
include("model.jl")
include("environment.jl")
include("trees.jl")
include("dataModels.jl")
include("../visuals/snapshots.jl")
using .Model, .DataModels, .GenealogicalTree

function main(trial::Int, dataModels, opts)
  # .reset models
  DataModels.resetModels!(dataModels)

  # .reset graph properties
  resetGraph!(graph, cntns)

  # .initial population
  active = populate!(graph, cntns, opts.dims)

  # .do simulation and generate snapshot of final front
  snapshot = simulate!(graph, cntns, active, opts, dataModels)
  sources = Model.graphSources(opts)

  if opts.detailed_analytics
    phylo, branchPoints, branchTimes = GenealogicalTree.genealogy(graph, sources)
    geneticLabels = getfield.(graph[sources], :ID_2)

    # jldopen(data_path * "/data_trees.jld2", "a+") do file
    #   ng = JLD2.Group(file, "trial_$trial")
    #   ng["tree"] = GenealogicalTree.buildTree(phylo, geneticLabels, branchPoint, branchTimes, sources)
    # end

    jldopen(data_path * "/data_phylo.jld2", "a+") do file
      ng = JLD2.Group(file, "trial_$trial")
      ng["end_time"] = snapshot.time
      ng["source_labels"] = geneticLabels
      ng["phylogeny"] = phylo
      ng["branch_points"] = branchPoints
      ng["branch_times"] = branchTimes
    end

    jldopen(data_path * "/data_sectors.jld2", "a+") do file
      ng = JLD2.Group(file, "trial_$trial")
      ng["sector_sizes"] = DataModels.lateralSectorSize(graph, opts.dims)
    end

    jldopen(data_path * "/data_extras.jld2", "a+") do file
      ng = JLD2.Group(file, "trial_$trial")
      foreach(kv -> ng[String(first(kv))] = last(kv), DataModels.unpack(dataModels))
    end
  end

  # .snapshot of sectors, lineages, hotspots
  if trial == 1
    _phylo = opts.detailed_analytics ? phylo : nothing
    landVisual(opts.dims, graph, _phylo; drawScatter=true, data_path=data_path)
    landscape(graph, opts.dims, data_path)
  end
end

#*******
# .command line options (input)
opts = parseArgs()

# .set seed for reproducabilty
rngSeed = opts.rng_seed == 0 ? rand() : opts.rng_seed
Model.StatsBase.Random.seed!(rngSeed)

# .set landscape from options params, optinally read in env
env, objs, opts = createEnvironment!(opts)

# .create containers for observales and associated functions
dataModels = DataModels.predefinedModels(opts.animate)

# .send command line parameters to stdout and mkdirs
data_path = Model.setPath(opts)

# .build hex map
graph, cntns = hexGraph(opts.dims)
Model.resetGraph!(graph, cntns, env)

# .send env to garbage collection
env = nothing

# .delete item from namedtuple
delete(nt::NamedTuple{names}, keys::Vector{Symbol}) where {names} = NamedTuple{filter(x -> x ∉ keys, names)}(nt)
delete(nt::NamedTuple{names}, keys::Symbol) where {names} = NamedTuple{filter(x -> x != keys, names)}(nt)

# .save all data in files
jldsave(data_path * "/objects.jld2"; objs=objs)
jldsave(data_path * "/Opts.jld2"; delete(merge(opts, (rngSeed=rngSeed,)), :data_path)...)

open(joinpath(data_path, "Opts.csv"), "w") do output
  options = delete(merge(opts, (rngSeed=rngSeed,)), :data_path)
  write(output, join(keys(options), ";") * "\n")
  write(output, join(values(options), ";") * "\n")
end

# .run all trials
for trial in range(1, opts.numberTrials)
  main(trial, dataModels, opts)
end