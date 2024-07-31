module LatticeDijkstra
include("../common/indexTools.jl")
using FileIO, JLD2, StaticArrays, DataStructures, Parameters

mutable struct Node
  weight::Float32
  ancestor::Vector{UInt32}
  closed::Bool
  ancestors::Int32
  Node(w=Inf, ancst=zeros(UInt32, 6)) = new(w, ancst, false, 0)
end

@with_kw mutable struct Result
  ancestorCounts::Union{Missing,Vector{Float64}} = missing
  lineageMap::Union{Missing,Matrix{Float64}} = missing
  totalSurviving::Union{Missing,Vector{UInt32}} = missing
  MSD::Union{Missing,Vector{Float64}} = missing
end

struct PathSet
  size::Float64
  source::UInt32
  tgts::Vector{UInt32}
  traces::Dict{UInt32,Vector{UInt32}}
end

graphSources(opts::NamedTuple) = graphSources(opts.ref_line, opts.width, opts.height)
function graphSources(ref_line, width, height)
  ref_line = iszero(ref_line) ? height : ref_line
  sources = collect(1:width) .+ width * (ref_line - 1)
  return sources
end

periodicBoundaries(w::Int, x::Int) = (x > cld(w, 2)) ? w - x : x
periodicBoundaries(w::Int, x::Float64) = (x > /(w, 2)) ? w - x : x

_reset!(graph::Vector{Node}) = @inbounds @simd for i in eachindex(graph)
  graph[i].closed = false
  graph[i].weight = Inf
  graph[i].ancestor .= 0
  graph[i].ancestors = 0
end

function dijkstra!(graph, connections, source, env)
  _reset!(graph)
  graph[source].weight = 0

  active = PriorityQueue{UInt32,Float32}()
  active[source] = 0.0

  @fastmath while ~isempty(active)
    #. pick min entry
    nodeIdx = dequeue!(active)

    graph[nodeIdx].closed = true
    weight = graph[nodeIdx].weight

    #. search around node, and update neighbors
    @inbounds for n in range(1, 6)
      nn = connections[n, nodeIdx]
      (nn == 0 || graph[nn].closed) && continue

      oldCost = graph[nn].weight
      newCost = env[nodeIdx] + weight

      if newCost == oldCost
        graph[nn].ancestors += 1
        graph[nn].ancestor[graph[nn].ancestors] = nodeIdx
        # push!(graph[nn].ancestor, nodeIdx)
      elseif newCost < oldCost
        active[nn] = newCost
        graph[nn].weight = newCost
        # graph[nn].ancestor = [nodeIdx]
        graph[nn].ancestor[1] = nodeIdx
        graph[nn].ancestors = 1
      end
    end
  end
  return nothing
end

_addPathNode!(trackSet, position, ancestor) =
  if haskey(trackSet, position)
    if ancestor ∉ trackSet[position]
      push!(trackSet[position], ancestor)
    end
  else
    trackSet[position] = UInt32[ancestor]
  end

function tracePaths!(trackSet::Dict{UInt32,Vector{UInt32}}, graph, targets::Vector{<:Integer})
  for i in eachindex(targets)
    current::UInt32 = targets[i]

    @inbounds while current != 0
      nancestors = graph[current].ancestors
      nancestors == 0 && break

      #. randomly select ancestor (all of equal weight)
      _ancestor = graph[current].ancestor[rand(1:nancestors)]

      _ancestor == 0 || _addPathNode!(trackSet, current, _ancestor)
      current = _ancestor
    end
  end
  return nothing
end

function optPaths(dims, graph, cntns, srcs::Vector{<:Int}, tgts, env, nth; cutoff=1.05)
  trackSets = Dict{UInt32,Vector{UInt32}}()
  # sinks = Vector{Vector{UInt32}}(undef, length(srcs))
  sinks = zeros(Bool, first(dims))
  fastestSinks = zeros(Bool, first(dims))

  for j in eachindex(srcs)
  # begin j = 100
    optPaths!(trackSets, sinks, fastestSinks, graph, cntns, srcs[j], tgts, env, nth; cutoff=cutoff)
    # weightSet, trackSet = optPaths(trackSets, sinks, graph, cntns, srcs[j], tgts, env, nth)
    # push!(trackSets, trackSet)
    # push!(weightSets, weightSet)
  end
  return sinks, fastestSinks, trackSets
end

function optPaths!(trackSet, sinks, fastestSinks, graph, cntns, src::Int, tgts, env, nth; cutoff=1.05)
  dijkstra!(graph, cntns, src, env)

  arrivalTimes = getfield.(graph[tgts], :weight)

  #. use all paths which differ at most by cutoff from the min time, and at most n paths
  sorter = sortperm(arrivalTimes)

  #. group tgts by similar (sorted) weight 
  grouped_tgts, grouped_times = grouparray(arrivalTimes[sorter], tgts[sorter])
  # @info arrivalTimes[sorter][1:10]
  # @info grouped_times[1:2]

  #. randomly select representative element in class
  # representatives = rand.(grouped_tgts)
  representative_times = maximum.(grouped_times)
  _maxIndex = findfirst(representative_times .> (cutoff * first(representative_times)))
  max_index = isnothing(_maxIndex) ? length(representative_times) : _maxIndex

  #. trace paths of nth smallest weights
  tracePaths!(trackSet, graph, collect(Iterators.flatten(grouped_tgts[1:max_index])))
  # tracePaths!(trackSet, graph, representatives[1:max_index])

  #. -> (weight, tgts with this weight)
  _sinks = grouped_tgts[1:max_index]

  #. save the targets reached the "fastest"
  _fastest = first(grouped_tgts)

  @inbounds for i in eachindex(_sinks)
    for j in eachindex(_sinks[i])
      sinks[_sinks[i][j]] = true
    end
  end

  @inbounds for j in eachindex(_fastest)
    fastestSinks[_fastest[j]] = true
  end
  # weightSet = Tuple{Float64,Vector{UInt32}}[]
  # foreach(t -> push!(weightSet, (weights[maximum(t)], t)), bottomSegments[1:max_index])
  # return weightSet, trackSet
  return nothing
end

function grouparray(x, y)
  N = [[t] for t in y]
  groupedX = [[t] for t in x]
  w = fill(true, length(x))
  for i in range(2, length(x))
    if isapprox(x[i - 1], x[i]; atol=0.005) && w[i - 1]
      append!(N[i], N[i - 1])
      append!(groupedX[i], groupedX[i - 1])
      w[i - 1] = false
    end
  end
  return N[w], groupedX[w]
end

function lineageMSD(h, w, rl, srcs, phylo, rev=false; interval=1)
  MSD = zeros(Float64, h)
  indices = zeros(UInt32, 2 * h)

  #. select lineages from intervals of n
  reducedFront = srcs[1:interval:end]

  # .build r_avg(y)
  r_y = zeros(Float64, h)
  _msd = zeros(Float64, h)

  for src in reducedFront
    current = src
    posMarker = 0

    @inbounds while current != 0
      posMarker += 1
      indices[posMarker] = current
      current = phylo[current]
    end

    # .build r_avg(y)
    x, y = reconstructCoordinatesReal(indices[1:posMarker], w, h)

    @inbounds for i in eachindex(x)
      r_y[Int64(y[i])] = x[i]
    end

    # x_ref = rev ? last(r_y) : first(r_y)
    # dx = periodicBoundaries.(w, r_y .- x_ref)

    # .http://www.pnas.org/cgi/doi/10.1073/pnas.0710150104
    @inbounds for L in 1:5:rl
      _counter = 0
      for y₀ in rl:-50:L
        dx = periodicBoundaries.(w, view(r_y, (y₀ - L):y₀) .- r_y[y₀])
        _msd[L] += sum(dx .^ 2) / L
        _counter += 1
      end
      _msd[L] /= _counter
    end

    MSD += _msd
    # .tortuosity calculation
    path_separation = sqrt((x[1] - x[posMarker])^2 + (y[1] - y[posMarker])^2)
    push!(tortuosity, posMarker / path_separation)
  end

  #? what if we apply a sliding window here to get the sliding window values
  return Results(; MSD=MSD ./ (length(reducedFront)), tortuosity=mean_and_std(tortuosity))
  # return Results(; MSD=_msd ./ (_counts .+ (_counts .== 0)), tortuosity=mean_and_std(tortuosity))
end
end
