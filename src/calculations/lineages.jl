module LineageTracing
include("../common/indexTools.jl")
using JLD2, NearestNeighbors, StaticArrays
import StatsBase: mean_and_std

@kwdef mutable struct Results
  ancestorCounts::Union{Missing,Vector{Float64}} = missing
  lineageMap::Union{Missing,Matrix{Float64}} = missing
  MSD::Union{Missing,Vector{Float64}} = missing
  tortuosity::Union{Missing,NTuple{2,Float64}} = missing
end

periodicBoundaries(w::Int, x::Int) = (abs(x) > cld(w, 2)) ? sign(x) * (w - abs(x)) : x
periodicBoundaries(w::Int, x::Float64) = (abs(x) > /(w, 2)) ? sign(x) * (w - abs(x)) : x

mask(haystack::Vector{String}, needle::String) = filter(t -> contains(t, needle), haystack)
mask(haystack::Base.KeySet, needle::String) = filter(t -> contains(t, needle), collect(haystack))

# doc: trace lineages individually, removing overlaps, and map onto plane
mapLineageTrace!(target::Vector{<:Real}, phylo) = target[collect(keys(phylo))] .+= 1

# doc: trace lineages individually and map onto plane
mapLineageTrace!(target::Vector{<:Real}, srcs, phylo) =
  for src in srcs
    current = src
    @inbounds while current != 0
      target[current] += 1
      current = phylo[current]
    end
  end

graphSources(opts::NamedTuple) = graphSources(opts.ref_line, opts.width, opts.height)
function graphSources(ref_line, width, height)
  ref_line = iszero(ref_line) ? height : ref_line
  sources = collect(1:width) .+ width * (ref_line - 1)
  return sources
end

function lineageTraces(path, opts; full=true, name="phylogeny")::Results
  lineageMap = zeros(Float64, opts.width * opts.height)
  ancestorCount = zeros(Float64, opts.width)

  sources = graphSources(opts)
  jldopen(joinpath(path, "data_phylo.jld2"), "r") do file
    for trial in mask(keys(file), "trial")
      phylo = file[trial][name]
      geneticLabels = file[trial]["source_labels"]
      full ? mapLineageTrace!(lineageMap, sources, phylo) : mapLineageTrace!(lineageMap, phylo)

      _uniqueLabels = collect(unique(geneticLabels))
      ancestorCount[_uniqueLabels] .+= 1
    end
  end
  return Results(; ancestorCounts=ancestorCount, lineageMap=reshape(lineageMap, (opts.width, opts.height)))
end

# doc: remove lateral moving lineages at the front
function reduceBoundary(srcs::Vector{<:Int}, ref_line::Int, phylo::Dict{UInt32,UInt32}, width::Int, height)
  _reducedBoundary = UInt32[]

  for src in srcs
    current = src
    _index = src

    # task: construct array of lineage trace then find last instance of crossing
    @inbounds while phylo[current] != 0
      y = last(reconstructCoordinates(current, width, height))

      y >= ref_line && (_index = current)
      (y <= ref_line - 3) && break

      current = phylo[current]
    end

    push!(_reducedBoundary, _index)
  end
  sort!(_reducedBoundary)
  return unique!(_reducedBoundary)
end

function lineageMSD(h, w, rl, srcs, phylo, rev=false; interval=1)
  MSD = zeros(Float64, h)
  # _counts = zeros(UInt32, h)
  tortuosity = Float64[]

  indices = zeros(UInt32, 2 * h)

  #. select lineages from intervals of n
  reducedFront = sort(reduceBoundary(srcs, rl, phylo, w, h))[1:interval:end]

  # .build r_avg(y)
  r_y = zeros(Float64, h)
  _msd = zeros(Float64, h)

  for src in reducedFront
    # if ~rev
    #   x_ref, y_ref = reconstructCoordinatesReal(current, w, h)
    # else
    #   @inbounds while phylo[current] != 0
    #     current = phylo[current]
    #   end
    #   x_ref, y_ref = reconstructCoordinatesReal(current, w, h)
    #   current = src
    # end

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

    # .http://www.pnas.org/cgi/doi/10.1073/pnas.0710150104
    # @inbounds for L in 1:5:rl
    #   _counter = 0
    #   for y₀ in 1:50:(rl - L)
    #     dx = periodicBoundaries.(w, view(r_y, y₀:(y₀ + L)) .- r_y[y₀])
    #     _msd[L] += sum(dx .^ 2) / L
    #     _counter += 1
    #   end
    #   _msd[L] /= _counter
    # end
    # MSD .+= _msd

    #. set reference point to be the bottom
    dx = periodicBoundaries.(w, view(r_y, 1:h) .- r_y[1])
    MSD += (dx .^2)

    # .tortuosity calculation
    path_separation = sqrt((x[1] - x[posMarker])^2 + (y[1] - y[posMarker])^2)
    push!(tortuosity, posMarker / path_separation)
  end

  return Results(; MSD=MSD ./ (length(reducedFront)), tortuosity=mean_and_std(tortuosity))
end

function lineageMSD_window(h, w, rl, srcs, phylo, rev=false; y_max=750, interval=1)
  MSD = zeros(Float64, y_max)
  tortuosity = Float64[]

  indices = zeros(UInt32, 3 * h)

  _reducedBoundary = reduceBoundary(srcs, rl, phylo, w, h)
  for src in _reducedBoundary
    current = src
    posMarker = 0

    @inbounds while current != 0
      posMarker += 1
      indices[posMarker] = current
      current = phylo[current]
    end

    x, y = reconstructCoordinatesReal(indices[1:posMarker], w, h)

    # .build r_avg(y)
    r_y = zeros(Float64, h)

    @inbounds for i in eachindex(x)
      # ?: how to handle mutable values with periodic boundaries
      r_y[Int64(y[i])] = x[i]
    end

    # .http://www.pnas.org/cgi/doi/10.1073/pnas.0710150104
    _msd = zeros(Float64, y_max)
    @inbounds for L in 1:y_max
      for y₀ in 1:(h - L)
        dx = periodicBoundaries.(w, abs.(view(r_y, y₀:(y₀ + L)) .- r_y[y₀]))
        _msd[L] += sum(dx .^ 2) / L
      end
      _msd[L] /= (h - L)
    end

    MSD .+= _msd
    push!(tortuosity, posMarker / sqrt((x[1] - x[posMarker])^2 + (y[1] - y[posMarker])^2))
  end
  return Results(; MSD=MSD ./ length(_reducedBoundary), tortuosity=mean_and_std(tortuosity))
end

# doc: surviving ancestors as the roots of the phylogeny trees
function phylogenyTreeRoots(phylo, width)
  # .find all keys with indices less than 'width', corresponding to bottom edge
  return sort(collect(filter(t -> t <= width, keys(phylo))))
end

# doc: use spatial partition to find visited environmental objects
function objectVisits(generators::Vector{NTuple{2,Float64}}, input, nTrials, width, height, radius)::Matrix{Int64}
  reached = zeros(Int64, (length(generators), nTrials))
  kdtree = KDTree(SVector.(generators))

  jldopen(joinpath(input, "data_phylo.jld2"), "r") do file
    for (i, trial) in enumerate(mask(keys(file), "trial"))
      phylo = file[trial]["phylogeny"]
      x, y = reconstructCoordinatesReal(collect(keys(phylo)), width, height)
      idxs, dist = nn(kdtree, SVector.(zip(x, y)))
      @inbounds foreach(j -> reached[j, i] = 1, idxs[dist .<= radius])
    end
  end
  return reached
end

end