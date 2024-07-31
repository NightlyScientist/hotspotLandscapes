module ContinuousOptPaths
using StatsBase, SplitApplyCombine, Parameters
using Graphs, SimpleWeightedGraphs, Distances, DelaunayTriangulation
export Results, VoronoiGraph, voronoiTriangulation, weightedGraph, allOptimalPaths, sourceSinkConnections, OptimalPath

@with_kw mutable struct Results
  ancestorCounts::Union{Missing,Vector{Float64}} = missing
  lineageMap::Union{Missing,Matrix{Float64}} = missing
  totalSurviving::Union{Missing,Vector{UInt32}} = missing
  MSD::Union{Missing,Vector{Float64}} = missing
end

@with_kw mutable struct VoronoiGraph
  connections::Union{Missing,Dict{Int,Vector{UInt32}}} = missing
  distance::Union{Missing,Dict{Int,Vector{Float64}}} = missing
  edgeGenerators::Union{Missing,Vector{UInt32}} = missing
  voronoi::Union{Missing,Any} = missing
  generators::Union{Missing,Dict{Int64,NTuple{2,Float64}}} = missing
end

struct OptimalPath
  pathDistance::Float64
  source::UInt32
  sink::UInt32
  pathPositions::Vector{NTuple{2,Float64}}
  pathIndices::Vector{UInt32}
end

periodicBoundaries(w::Int, x::Int) = (x > cld(w, 2)) ? w - x : x
periodicBoundaries(w::Int, x::Float64) = (x > /(w, 2)) ? w - x : x

function graphConnections(_voronoi)
  neighbors = Dict{Int,Vector{Int}}()
  neighbors_distances = Dict{Int,Vector{Float64}}()

  for gen_ref_index in each_generator(_voronoi)
    _neighbors = Int[]
    _neighbors_distances = Float64[]
    _polygon = get_polygon(_voronoi, gen_ref_index)
    pos_1 = get_generator(_voronoi, gen_ref_index)

    for gen_index in each_generator(_voronoi)
      gen_index == gen_ref_index && continue
      poly_indices = get_polygon(_voronoi, gen_index)

      if isneighboring(poly_indices, _polygon)
        pos_2 = get_generator(_voronoi, gen_index)
        dist = sqrt(sum(abs2, pos_1 .- pos_2))

        push!(_neighbors, gen_index)
        push!(_neighbors_distances, dist)
      end
    end
    neighbors[gen_ref_index] = _neighbors
    neighbors_distances[gen_ref_index] = _neighbors_distances
  end
  return VoronoiGraph(; connections=neighbors, distance=neighbors_distances)
end

function graphPosisions(top_edge, bottom_edge, generators, opts)
# function graphPosisions(generators, opts)
  seg_1 = length(generators)
  seg_2 = length(top_edge)
  # seg_2 = opts.width
  positions = Dict{Int,Tuple{Float64,Float64}}()

  _topKeys = collect(keys(top_edge))
  _btmKeys = collect(keys(bottom_edge))
  foreach(t -> positions[t] = generators[t], eachindex(generators))
  foreach(t -> positions[t + seg_1] = (_topKeys[t], opts.ref_line), eachindex(_topKeys))
  foreach(t -> positions[t + seg_1 + seg_2] = (_btmKeys[t], 1), eachindex(_btmKeys))
  # foreach(t -> positions[t + seg_1] = (t, opts.ref_line), 1:(opts.width))
  # foreach(t -> positions[t + seg_1 + seg_2] = (t, 1), 1:(opts.width))
  return positions
end

function edgeNodeIndices(top_edge, btm_edge, generators)
# function edgeNodeIndices(generators, opts)
# seg_2 = opts.width
  seg_1 = length(generators)
  seg_2 = length(top_edge)
  # sources = sort(collect(keys(top_edge))) .+ seg_1
  # sinks = sort(collect(keys(btm_edge))) .+ seg_1 .+ seg_2
  sources = collect((1:length(top_edge)) .+ seg_1)
  sinks = collect((1:length(btm_edge)) .+ seg_1 .+ seg_2)
  # sources = collect((1:(opts.width)) .+ seg_1)
  # sources = collect((1:(opts.width)) .+ seg_1)
  # sinks = collect((1:(opts.width)) .+ seg_1 .+ seg_2)
  return sources, sinks
end

#doc: use flyod warshall algorithm or dijkstra to find shortest paths
function allOptimalPaths(sources, sinks, opts, positions, g; n=25, cutoff=1, floyd=true)
  optimalPathSets = Dict{Int,Vector{OptimalPath}}()
  floyd && (optPaths = floyd_warshall_shortest_paths(g))

  for sourceIndex in sources
    _temp = []
    floyd || (optPaths = dijkstra_shortest_paths(g, sourceIndex))

    for sinkIndex in sinks
      pathIndexSets = floyd ? enumerate_paths(optPaths, sourceIndex, sinkIndex) : enumerate_paths(optPaths, sinkIndex)

      idxSets = typeof(pathIndexSets) <: Vector{Vector{Int}} ? pathIndexSets : [pathIndexSets]

      # .iterate through all paths of the same shortest time
      for pathIndexSet in idxSets
        _pathTrace, idxs = tracePath(positions, opts, pathIndexSet)
        _pathDistance = pathDistance(g, idxs)
        push!(_temp, (_pathDistance, sourceIndex, sinkIndex, _pathTrace, pathIndexSet))
      end
    end
    # .use all paths which differ at most by cutoff from the min time, and at most n paths
    sortedPaths = sort(_temp; by=first)
    _max_traces = findfirst(first.(sortedPaths) .> (cutoff * first(first(sortedPaths))))
    max_traces = (~isnothing(_max_traces) && _max_traces < n) ? _max_traces : n
    optimalPathSets[sourceIndex] = map(t -> OptimalPath(t...), sortedPaths[1:max_traces])
  end
  return optimalPathSets
end

function tracePath(positions::Dict, opts, idxs::Vector{Int})::Tuple{Vector{NTuple{2,Float64}},Vector{Int}}
  segments = NTuple{2,Float64}[]
  width = opts.width

  for j in eachindex(idxs)
    j == firstindex(idxs) && continue
    a = idxs[j - 1]
    b = idxs[j]

    pos_1 = get(positions, a, nothing)
    pos_2 = get(positions, b, nothing)
    x1, y1 = pos_1
    x2, y2 = pos_2

    if abs(x2 - x1) < width / 2
      push!(segments, pos_1, pos_2)
    else
      if x2 - x1 <= 0
        m = (y2 - y1) / (width + x2 - x1)
        mid_y = m * (width - x1) + y1
        push!(segments, pos_1, (width, mid_y))
        push!(segments, (0, mid_y), pos_2)
      else
        m = (y2 - y1) / (-width + x2 - x1)
        mid_y = m * (-x1) + y1
        push!(segments, pos_1, (0, mid_y))
        push!(segments, (width, mid_y), pos_2)
      end
    end
  end
  return segments, idxs
end

function pathDistance(g, idxs)
  dist = 0
  for j in eachindex(idxs)
    j == firstindex(idxs) && continue
    dist += SimpleWeightedGraphs.get_weight(g, idxs[j], idxs[j - 1])
  end
  return dist
end

function regionCentroids(source_sinks)
  _centroids = mean.(collect(group(isequal, (last ∘ Tuple).(source_sinks))))
  return floor.(Int, _centroids)
end

function sourceSinkConnections(vGraph, opts; ref_line=Inf, n_near=2, interval=1)
  # task: handle periodic boundary conditions
  x = collect(1:interval:(opts.width))
  # y1 = ones(Int, opts.width)
  y1 = ones(Int, length(x))

  height = isfinite(ref_line) ? ref_line : opts.height
  y2 = fill(height, length(x))
  # y2 = fill(height, opts.width)

  _keys = collect(keys(vGraph.generators))

  _centers = filter(t -> last(t) <= ref_line, collect(values(vGraph.generators)))
  centerX = getfield.(_centers, 1)
  centerY = getfield.(_centers, 2)

  D = pairwise(Euclidean(), hcat(x, y2)', hcat(centerX, centerY)')
  # (_, top) = findmin(D; dims=2)
  # top = zeros(Int, (length(x), n_near))
  top = Dict{Int,Vector{Integer}}()
  for row in axes(D, 1)
    _mask = sortperm(view(D, row, :))[1:n_near]
    # top[row, :] .= _keys[_mask]
    top[x[row]] = _keys[_mask]
  end

  D = pairwise(Euclidean(), hcat(x, y1)', hcat(centerX, centerY)')
  # (_, btm) = findmin(D; dims=2)
  # btm = zeros(Int, (length(x), n_near))
  btm = Dict{Int,Vector{Integer}}()
  for row in axes(D, 1)
    _mask = sortperm(view(D, row, :))[1:n_near]
    # btm[row, :] .= _keys[_mask]
    btm[x[row]] = _keys[_mask]
  end
  return btm, top
end

function edgeWeight(pos_1, pos_2, radius, strength, Lx)
  dx, dy = abs.(pos_1 .- pos_2)
  dx = dx < (Lx / 2) ? dx : Lx - dx
  λ = sqrt(dx^2 + dy^2)

  λ >= 2 * radius && return 2 * radius / (strength + 1) + (λ - 2 * radius)
  return λ / (strength + 1)
end

function maximumPairDistance(neighbor_list_dist)
  max_gen_pair_dist = 0
  for (_, v) in neighbor_list_dist
    _max = maximum(v)
    max_gen_pair_dist < _max && (max_gen_pair_dist = _max)
  end
  return max_gen_pair_dist
end

function connectBoundary!(vGraph, opts, max_gen_pair_dist)
  edge_generators = vGraph.edgeGenerators
  for i in eachindex(edge_generators)
    ref_gen_index = edge_generators[i]
    pos_1 = get_generator(vGraph.voronoi, ref_gen_index)

    # for gen_index in edge_generators
    for j in (i + 1):lastindex(edge_generators)
      gen_index = edge_generators[j]
      ref_gen_index == gen_index && continue
      pos_2 = get_generator(vGraph.voronoi, gen_index)

      dx, dy = abs.(pos_1 .- pos_2)
      dx = (dx < opts.width / 2) ? dx : opts.width - dx
      dist = sqrt(abs2(dx) + abs2(dy))

      if dist <= max_gen_pair_dist
        push!(vGraph.connections[ref_gen_index], gen_index)
        push!(vGraph.distance[ref_gen_index], dist)

        push!(vGraph.connections[gen_index], ref_gen_index)
        push!(vGraph.distance[gen_index], dist)
      end
    end
  end
end

function isneighboring(a, b)::Bool
  _a = sort(a)
  _b = sort(b)

  first(_a) == first(_b) && return true

  itr_1::Int = firstindex(_a)
  itr_2::Int = firstindex(_b)
  whch = _a[itr_1] > _b[itr_2] ? 1 : 2

  while itr_1 < length(_a) && itr_2 < length(_b)
    if whch == 1
      itr_2 += 1
    else
      itr_1 += 1
    end
    _a[itr_1] == _b[itr_2] && return true
    whch = _a[itr_1] > _b[itr_2] ? 1 : 2
  end
  return false
end

function add_graph_edges!(g, gen_index, index, dist)
  Graphs.add_edge!(g, gen_index, index)
  Graphs.add_edge!(g, index, gen_index)

  g.weights[gen_index, index] = dist
  g.weights[index, gen_index] = g.weights[gen_index, index]
  return nothing
end

#doc construct network of hotspot using totally connected or NN 
function hotspotNetwork(voronoiGraph, g, opts; NN = true)
  for (gen_index, pos_1) in voronoiGraph.generators
    indices = NN ? voronoiGraph.connections[gen_index] : range(1,length(voronoiGraph.generators))
    for neighbor_index in indices
      gen_index == neighbor_index && continue
      pos_2 = voronoiGraph.voronoi.generators[neighbor_index]
      dist = edgeWeight(pos_1, pos_2, opts.radius, opts.intensity, opts.width)
      add_graph_edges!(g, gen_index, neighbor_index, dist)
    end
  end
  return nothing
end

# doc: graph constructed by connecting all generators together
function weightedGraph(voronoiGraph, top_edge, bottom_edge, opts; NN=true)
  nv = length(voronoiGraph.generators) + length(top_edge) + length(bottom_edge)
  # nv = length(voronoiGraph.generators) + 2 * size(top_edge, 1)
  g = SimpleWeightedDiGraph{Int,Float64}(nv)

  # .connect hotspots together
  hotspotNetwork(voronoiGraph, g, opts; NN = NN)

  seg_1 = length(voronoiGraph.voronoi.generators)
  seg_2 = length(top_edge)
  # .add edges to graph by connecting them to near n voronoi generators
  # for i in axes(top_edge, 1)
  _topKeys = collect(keys(top_edge))
  _btmKeys = collect(keys(bottom_edge))

  # for i in keys(top_edge)
  for i in eachindex(_topKeys)
    for nearestGenerator in top_edge[_topKeys[i]]
    # for nearestGenerator in top_edge[i, :]
      pos_1 = voronoiGraph.voronoi.generators[nearestGenerator]
      pos_2 = (convert(Int64, _topKeys[i]), opts.ref_line)
      dist = sqrt(sum(abs2, pos_1 .- pos_2))
      dist = dist >= opts.radius ? dist - opts.radius + opts.radius / (1 + opts.intensity) : dist

      index = i + seg_1
      # index = i + length(voronoiGraph.voronoi.generators
      add_graph_edges!(g, nearestGenerator, index, dist)
    end

    for nearestGenerator in bottom_edge[_btmKeys[i]]
    # for nearestGenerator in bottom_edge[i, :]
      pos_1 = voronoiGraph.voronoi.generators[nearestGenerator]
      pos_2 = (convert(Int64, _btmKeys[i]), 1)
      dist = sqrt(sum(abs2, pos_1 .- pos_2))
      dist = dist >= opts.radius ? dist - opts.radius + opts.radius / (1 + opts.intensity) : dist

      index = i + seg_2 + seg_1
      # index = i + opts.width + length(voronoiGraph.voronoi.generators)
      add_graph_edges!(g, nearestGenerator, index, dist)
    end
  end
  return g
end

function voronoiTriangulation(points, opts)::VoronoiGraph
  _voronoi = voronoi(triangulate(points); clip=true)

  _graph = graphConnections(_voronoi)
  _graph.voronoi = _voronoi

  # .edge voronoi cells and generator positions
  _graph.edgeGenerators = collect(DelaunayTriangulation.get_boundary_polygons(_voronoi))
  _graph.generators = DelaunayTriangulation.get_generators(_voronoi)

  # .connect generators within maximum(pairwise distances)
  maxGenPairDist = ContinuousOptPaths.maximumPairDistance(_graph.distance)
  connectBoundary!(_graph, opts, maxGenPairDist)

  return _graph
end

# doc: the terminal points of the optimal paths should represent the prediced surviving ancestors
function pathRoots(optimalPathSets)
  roots = Float64[]
  for (_, optimalPathSet::Vector{OptimalPath}) in optimalPathSets
    for optimalPath::OptimalPath in optimalPathSet
      x_position, _ = last(optimalPath.pathPositions)
      push!(roots, x_position)
    end
  end
  return unique(Int64.(roots))
end

# doc: msd using an interpolation between "knots" of the path trace 
function optimalPathMSD(optimalPathSets, ref_line, width)
  #task measure MSd from the bottom instead of the top
  _msd = zeros(Float64, ref_line)
  _counts = zeros(UInt32, ref_line)

  for (_, optimalPathSet::Vector{OptimalPath}) in optimalPathSets
    for optimalPath::OptimalPath in optimalPathSet
      pathPositions = optimalPath.pathPositions

      x_reference = first(first(pathPositions))
      y_reference = last(first(pathPositions))

      for j in 2:2:length(pathPositions)
        (x1, y1) = pathPositions[j - 1]
        (x2, y2) = pathPositions[j]
        if y1 != y2
          y_range = min(y1, y2):1:max(y1, y2)

          dx = @. x_reference - x1 - (x2 - x1) / (y2 - y1) * (y_range - y1)
          dx = periodicBoundaries.(width, dx)
          y = @. round(Int, abs(y_range - y_reference)) .+ 1
          _msd[y] .+= abs.(dx)
          _counts[y] .+= 1
        else
          y = round(Int64, y1)
          _msd[y] += abs.(x2 - x1)
          _counts[y] += 1
        end
      end
    end
  end
  return Results(; MSD=_msd ./ (_counts .+ (_counts .== 0)))
end

# doc: msd using on the "knots" of the path trace (no interpolation)
function optimalPathMSDAlt(optimalPathSets, ref_line, width)
  _msd = zeros(Float64, ref_line)
  _counts = zeros(UInt32, ref_line)

  for (_, optimalPathSet::Vector{OptimalPath}) in optimalPathSets
    for optimalPath::OptimalPath in optimalPathSet
      pathPositions = optimalPath.pathPositions

      x_reference = first(first(pathPositions))
      y_reference = last(first(pathPositions))

      for (_, optimalPathSet::Vector{OptimalPath}) in optimalPathSets
        for optimalPath::OptimalPath in optimalPathSet
          pathPositions = optimalPath.pathPositions
          _pathTrace = unique(pathPositions)

          x_reference = first(first(pathPositions))
          y_reference = last(first(pathPositions))

          dx = abs.(first.(_pathTrace) .- x_reference)
          # .bin in unit steps along the vertical 
          dy = round.(Int64, y_reference .- last.(_pathTrace)) .+ 1
          _msd[dy] .+= abs.(periodicBoundaries.(width, dx))
          _counts[dy] .+= 1
        end
      end
    end
  end
  return Results(; MSD=_msd ./ (_counts .+ (_counts .== 0)))
end

# doc: get a list of visited generator indices
function objectVisits(VoronoiGraph, sourceIndices, optimalPathSets)
  visitedGenerators = Dict{UInt32,Bool}()
  for source_index in sourceIndices
    visitedGeneratorIndices = unique(vcat(getfield.(optimalPathSets[source_index], :pathIndices)...))
    filtered = filter(t -> t ∉ sourceIndices, visitedGeneratorIndices)
    foreach(t -> visitedGenerators[t] = true, filtered)
  end
  return collect(keys(visitedGenerators))
end

end