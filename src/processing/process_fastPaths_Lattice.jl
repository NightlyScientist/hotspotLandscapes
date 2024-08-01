using JLD2, FileIO, ArgParse, Base.Threads
include("../calculations/latticeDijkstra.jl")
include("../../src/base/model.jl")
include("../../src/base/environment.jl")
include("../../src/base/containers.jl")

function execute_main(data_path, radius, width, height, rates, ref_line, nfast; interval=15)
  #task need to update to work with updated lattice dijkstra module
  opts = namedtuple(load(joinpath(input, "Opts.jld2")))

  # .sinks and sources, sampling the front using an interval
  srcs = (collect(1:width) .+ width * (ref_line - 1))[1:interval:end]
  tgts = collect(1:width)

  # .build graph and connections
  graph, cntns = Model.hexGraph((width, height))

  # .get environmental objects (hotspots)
  if ispath(joinpath(data_path, "objects.jld2"))
    objs = load(joinpath(data_path, "objects.jld2"), "objs")
  else
    objs = load(joinpath(data_path, "data.jld2"), "htspts")
  end

  # .build env grid
  env, _ = applyObstacles!(objs, radius, width, height)

  # .build weights from rates
  env = map(t -> 1 / rates[t], env)

  # .perform calculation
  weights, tracks = LatticeDijkstra.optPaths(envObjs, graph, cntns, srcs, tgts, env, nfast)

  # .save data of shortests paths
  data_opt_path = joinpath(data_path, "_opt_paths.jld2")
  jldsave(data_opt_path; tracks=tracks, weights=weights, nfast=nfast)

  #.visited hotspots, MSD, and suriving ancestors
  visitedObjectsIndices = ContinuousOptPaths.objectVisits(voronoiGraph, sourceIndices, optimalPathSets)
  opt_results = ContinuousOptPaths.optimalPathMSD(optimalPathSets, opts.ref_line, opts.width)
  survivingAncestors = ContinuousOptPaths.pathRoots(optimalPathSets)

  # .save data to csv: distance | mean 
  CSV.write(
    joinpath(input, "processed_optimalPathMSD.csv"),
    DataFrame(; distance=collect(1:(opts.ref_line)), msd=opt_results.MSD)
  )
  jldsave(joinpath(input, "processed_visitedObjects_optimalPaths.jld2"); visitedObjectsIndices=visitedObjectsIndices)
  jldsave(joinpath(input, "processed_surivingAncestors_optimalPaths.jld2"); survivingAncestors=survivingAncestors)
end

opts = begin
  sts = ArgParseSettings()
  add_arg_table!(sts, "--parameter_space_table", Dict(:arg_type => String, :required => true))
  add_arg_table!(sts, "--n_fastest", Dict(:arg_type => Int, :required => true))
  add_arg_table!(sts, "--interval", Dict(:arg_type => Int, :required => true))
  add_arg_table!(sts, "--cutoff", Dict(:arg_type => Float64, :required => true))
  parse_args(sts)
end

df = CSV.read(opts["parameter_space_csv"], DataFrame)

sim = copy(df[1, :])
width = sim.width
height = sim.height
dims = (width, height)
radius = sim.radius
nfast = opts["n_fastest"]
interval = opts["interval"]

df = CSV.read(opts["parameter_space_table"], DataFrame)

@threads for single_trial in eachrow(df)
  ref_line = single_trial.ref_line
  num_trials = single_trial.numberTrials
  data_path = single_trial.path

  # rebuild rates
  rates = Float64[1.0, single_trial.intensity + 1]

  execute_main(data_path, radius, width, height, rates, ref_line, nfast; interval=interval)
end
