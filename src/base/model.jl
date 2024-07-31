module Model

using NearestNeighbors, DataStructures, StaticArrays
using StatsBase, Random
using ArgParse

include("containers.jl")
include("dataModels.jl")
include("trees.jl")

import .HashMaps: HashVec, add!, remove!
export hexGraph, simulate!, populate!, resetGraph!, parseArgs

function ArgParse.parse_item(::Type{NTuple{2,T}}, x::AbstractString) where {T}
  return Tuple(convert.(T, parse.(Float64, split(x, ','))))
end

function graphSources(opts)
  ref_line = iszero(opts.ref_line) ? opts.height : opts.ref_line
  sources = collect(1:(opts.width)) .+ opts.width * (ref_line - 1)
  return sources
end

namedtuple(d::Dict) = (; zip(Symbol.(keys(d)), values(d))...)

function parseArgs()
  addArgs!(sts, name; kwargs...) = add_arg_table!(sts, "--$(name)", Dict(kwargs))

  sts = ArgParseSettings()
  addArgs!(sts, "numberTrials"; arg_type=Int, default=1)
  addArgs!(sts, "numberSamples"; arg_type=Int, default=100)
  addArgs!(sts, "width"; required=true, arg_type=Int)
  addArgs!(sts, "height"; required=true, arg_type=Int)
  addArgs!(sts, "density"; required=true, arg_type=Float64)
  addArgs!(sts, "radius"; required=true, arg_type=Int64)
  addArgs!(sts, "intensity"; required=true, arg_type=Float64)
  addArgs!(sts, "gap"; arg_type=Int64, default=0)
  addArgs!(sts, "ref_line"; arg_type=Int64, default=0)
  addArgs!(sts, "rng_seed"; arg_type=Int, default=1)
  addArgs!(sts, "printInfo"; action=:store_true)
  addArgs!(sts, "rewrite"; action=:store_true)
  addArgs!(sts, "data_path"; required=true, arg_type=String)
  addArgs!(sts, "animate"; action=:store_true)
  addArgs!(sts, "parameter"; arg_type=String, default="missing")
  addArgs!(sts, "env_type"; arg_type=String, default="uniform")
  addArgs!(sts, "landscape"; arg_type=String)
  addArgs!(sts, "detailed_analytics"; action=:store_true)
  addArgs!(sts, "separation"; arg_type=Int)
  parsed_args = namedtuple(parse_args(sts))

  parsed_args.printInfo && foreach(k -> println("  > $k  =>  $(parsed_args[k])"), keys(parsed_args))
  return parsed_args
end

mutable struct Node
  const x::Float32
  const y::Float32
  ancestor::UInt32
  filled::Bool
  ID_1::Int32
  ID_2::Int32
  time::Float32
  nbors::UInt32
  Node(n, x, y) = new(x, y, 0, false, 0, 0, 0.0, n)
end

gNodes(typ::Symbol=:hex) = if typ == :hex
    # https://www.redblobgames.com/grids/hexagons/#neighbors
    bottom = [(1, 0), (-1, 0), (-1, -1), (0, -1), (-1, 1), (0, 1)]
    top = [(1, 0), (-1, 0), (0, 1), (1, 1), (0, -1), (1, -1)]
    return Dict{Int,Vector{NTuple{2,Int}}}(1 => bottom, 0 => top)
  end

function buildMap(dimensions, nodes, T::Type, cnstr)
  lx, ly = dimensions
  graph = Vector{T}(undef, lx * ly)
  connections = zeros(UInt32, 6, lx * ly)

  for row in 1:ly, col in 1:lx
    nodeIndx = lx * (row - 1) + col
    nbors::UInt32 = 0

    for (i, (dx, dy)) in enumerate(nodes[row % 2])
      ny = row + dy
      0 < ny <= ly || continue
      nx = mod(col + dx, 1:lx)
      idx = lx * (ny - 1) + nx
      connections[i, nodeIndx] = idx
      nbors += 1
    end
    graph[nodeIndx] = cnstr(nbors, col, row)
  end
  return graph, connections
end

hexGraph(dims) =
  let
    cvr = (x, y) -> (sqrt(3) * (x - 0.5 * (y % 2)), 1 + 1.5 * (y - 1))
    buildMap(dims, gNodes(), Node, (n, x, y) -> Node(n, cvr(x, y)...))
  end

function populate!(graph, cntns, dims; row=1, num=2)
  active = Vector{HashVec{UInt32}}(undef, num)
  foreach(i -> active[i] = HashVec{UInt32}(), collect(1:1:num))

  for col in 1:dims[1]
    # shift group affliation to match environment
    nodeIndx = dims[1] * (row - 1) + col

    graph[nodeIndx].filled = true
    graph[nodeIndx].ID_2 = col

    # don't add site to front if it already has no empty nbors
    graph[nodeIndx].nbors > 0 && add!(active[graph[nodeIndx].ID_1], nodeIndx)

    for nbor in view(cntns, :, nodeIndx)
      nbor == 0 && continue

      # subtract one from all neighbors
      graph[nbor].nbors -= 1

      # if this site, or neighbor, is surrounded, then remove it
      if graph[nbor].nbors == 0 && graph[nbor].filled
        remove!(active[graph[nbor].ID_1], nbor)
      end
    end
  end
  return active
end

@inline function nextIndex(graph, cntns, idx)
  nbors = UInt32[]
  @inbounds for nbor in view(cntns, :, idx)
    (nbor == 0 || graph[nbor].filled) && continue
    push!(nbors, nbor)
  end
  if length(nbors) == 0
    return 0
  elseif length(nbors) == 1
    return first(nbors)
  end
  return rand(nbors)
end

function chooseNext(rates::SVector{2,Float64}, active)::Tuple{Int,Float64}
  r₁::Float64 = rates[1] * length(active[1].data)
  a₀::Float64 = r₁ + rates[2] * length(active[2].data)
  η::Float64 = a₀ * rand()
  return (η <= r₁) ? 1 : 2, a₀
end

function simulate!(graph::Vector{Node}, cntns, active, opts, dataModels)
  rates = SVector{2,Float64}([1, max(1.0 + opts.intensity, 0.0)])
  width, height = opts.dims
  time = 0.0

  populationSizes = Int[length(active[i].data) for i in [1, 2]]

  nrecord::Int64 = cld(width * height, opts.numberSamples)
  itrCntr::Int64 = nrecord
  top_line = width * (height - 1) + 1
  continue_recording::Bool = true

  @inbounds @fastmath while true
    sum(populationSizes) == 0 && break

    if continue_recording && itrCntr == nrecord
      itrCntr = 0
      for (_, measure) in dataModels
        push!(measure.data, measure.func(active, graph, time))
      end
    end
    itrCntr += 1

    # find which cell will grow next
    (label, R) = chooseNext(rates, active)

    # update time
    time += -log(rand()) / R

    # get node indices
    parentIdx = rand(active[label].data)
    childIdx = nextIndex(graph, cntns, parentIdx)

    # fill new node, add to front if viable
    node = graph[childIdx]
    node.filled = true
    node.time = time
    node.ID_2 = graph[parentIdx].ID_2
    node.ancestor = parentIdx

    # don't add site to front if it already has no empty nbors
    if node.nbors > 0 && rates[node.ID_1] > 0.0
      add!(active[node.ID_1], childIdx)
      populationSizes[node.ID_1] += 1
    end

    # iterate through neighbors, updating nbor counts
    for nbor in view(cntns, :, childIdx)
      nbor == 0 && continue

      # subtract one from all neighbors
      graph[nbor].nbors -= 1

      # if this site, or neighbor, is surrounded, then remove it
      if graph[nbor].nbors == 0 && graph[nbor].filled
        remove!(active[graph[nbor].ID_1], nbor)
        populationSizes[graph[nbor].ID_1] -= 1
      end
    end
    childIdx >= top_line && (continue_recording = false)
  end
  return (time=time,)
end

resetGraph!(graph, cntns, env=missing) = @inbounds @simd for i in eachindex(graph)
    graph[i].filled = false
    graph[i].ancestor = 0
    graph[i].ID_2 = 0
    graph[i].time = 0.0
    graph[i].nbors = sum(view(cntns, :, i) .> 0)
    ismissing(env) || (graph[i].ID_1 = env[i])
  end

function setPath(opts)::String
  _opts = ["env_type", "width", "height", "density", "intensity", "radius", "numberTrials", "numberSamples", "rng_seed"]
  _opts = Symbol.(_opts)

  parsedOpts = []
  for opt in _opts
    haskey(opts, Symbol(opt)) || continue
    val = getfield(opts, Symbol(opt))
    if eltype(val) <: Float64
      val = round.(val, digits=3)
    end
    push!(parsedOpts, "$(opt)_$(val)")
  end

  # optional paramters
  opts.gap > 0 && push!(parsedOpts, "gap_$(opts.gap)")
  opts.ref_line > 0 && push!(parsedOpts, "rl_$(opts.ref_line)")
  opts.env_type == "circle" && push!(parsedOpts, "sep_$(opts.separation)")
  opts.detailed_analytics && push!(parsedOpts, "da")
  opts.animate && push!(parsedOpts, "animated")

  path = mkpath(opts.data_path * "/" * join(parsedOpts, ","))

  if ispath(path) && ~isempty(readdir(path))
    if opts.rewrite
      foreach(rm, filter(endswith(".csv"), readdir(path; join=true)))
      foreach(rm, filter(endswith(".png"), readdir(path; join=true)))
      foreach(rm, filter(endswith(".jld2"), readdir(path; join=true)))
      return path
    else
      @info " ! output path not empty"
      exit()
    end
  end
  return path
end

end