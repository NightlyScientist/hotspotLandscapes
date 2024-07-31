module DataModels
using StatsBase
include("trees.jl")
using .GenealogicalTree

mutable struct DataModel{T}
  data::Vector{T}
  func::Function
  DataModel{T}(func) where {T} = new{T}(Vector{T}(), func)
end

resetModels!(dm) =
  for symb in keys(dm)
    empty!(dm[symb].data)
  end

function predefinedModels(animation=false)::Dict{Symbol,DataModel}
  dataModels = Dict{Symbol,DataModel}()
  dataModels[:sectors] = DataModel{NTuple{3,Float32}}(sectorCoarsening)
  dataModels[:front] = DataModel{NTuple{5,Float32}}(frontMotion)

  if animation
    _type = Tuple{Float64,Vector{Int32},Vector{Int32}}
    dataModels[:animation] = DataModel{_type}(growAnimation)
  end
  return dataModels
end

# doc: generate list of pairs -> (symbol, data)
function unpack(dataModels)
  p = []
  foreach(x -> push!(p, Pair(x[1], x[2].data)), (dataModels))
  return p
end

# doc: measures the mean and variance of the front height
function frontMotion(active, graph, time)
  front = vcat(active[1].data, active[2].data)
  front_y = getfield.(graph[front], :y)
  m, v = mean_and_var(front_y)
  return (time, m, v, extrema(front_y)...)
end

# doc: measures the mean and variance of the sector size at the front, and number of surviving roots
function sectorCoarsening(active, graph, time)
  front = vcat(active[1].data, active[2].data)
  uniqueElems = countmap(getfield.(graph[front], :ID_2))
  # nroots = length(keys(uniqueElems))
  m, v = mean_and_var(values(uniqueElems))
  return (time, m, v)
  # return (time, m, v, nroots)
end

# doc: take snapshots of the front to later make an animation
function growAnimation(active, graph, time)
  lineages = Int32[]
  if time > 0.0
    front = vcat(active[1].data, active[2].data)
    phylogeny, _ = GenealogicalTree.genealogy(graph, front)
    lineages = collect(keys(phylogeny))
  end
  growth = getfield.(graph, :ID_2)
  return (time, growth, lineages)
end

function lateralSectorSize(graph, dims)
  sectorSizes = zeros(Float32, last(dims), 2)
  grid_2 = reshape(getfield.(graph, :ID_2), dims)

  for y in axes(grid_2, 2)
    uniqueElems = countmap(view(grid_2, :, y))
    # nroots = length(keys(uniqueElems))
    m, v = mean_and_var(values(uniqueElems))
    sectorSizes[y, 1] = m
    sectorSizes[y, 2] = v
    # sectorSizes[y, 3] = nroots
  end
  return sectorSizes
end

end