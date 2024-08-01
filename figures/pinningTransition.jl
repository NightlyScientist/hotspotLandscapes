#> this file generates images of the pinning transition, the overlap between optimal paths and lineages
using JLD2, FileIO, Revise, Parameters, CSV, DataFrames, Base.Threads
include("common/theme.jl")
include("common/binning.jl")
include("common/modeling.jl")
includet("../src/calculations/lineages.jl")

# import .Model: gNodes, buildMap
# typemap = Dict("Main.GraphModel.TreeNode" => TreeNode)
module PinningResults
using Measurements
@kwdef mutable struct Results
  lineageVisits::Union{Missing,Matrix{Float64}} = missing
  optimalVisits::Union{Missing,Vector{Integer}} = missing
  lineageAncestors::Union{Missing,Vector{Integer}} = missing
  optimalAncestors::Union{Missing,Vector{Integer}} = missing
  generators::Union{Missing,Vector{NTuple{2,Float64}}} = missing
  color::Union{Missing,Any} = missing
  parameters::Union{Missing,NamedTuple} = missing
end
@kwdef mutable struct Averaging
  ancestorOverlap::Union{Missing,Measurement{Float64}} = missing
  objectOverlap::Union{Missing,Measurement{Float64}} = missing
  color::Union{Missing,Any} = missing
  parameters::Union{Missing,NamedTuple} = missing
end
end

function backgroundEstimate(spinput, linput, xy, dims, density, radius; nl=1, ntests=5)
  bckgrndEstmt = Vector{Measurement}(undef, max(1, ntests))

  for i in range(1, max(1, ntests))
    objs = _uniformEnv(radius, density, dims...)

    # build kdtree from hotspot locations
    kdtree = buildHSKDT(xy, objs, dims)
    sphs = last(spTracks(spinput, xy, dims, kdtree, radius; nl=nl))

    acc, _ = subsetHotspots(kdtree, xy, radius, sphs, linput)
    # bckgrndEstmt[i] = (acc_hotspots = acc,freq_hotspots = frq)
    bckgrndEstmt[i] = acc
  end
  return mean_and_std(bckgrndEstmt)
end

function subsetAncestors(data, dims, roots)
  accuracy = Float64[]
  cnts = 0
  survivalFreq = zeros(Float64, dims[1])

  # interpolate values between surviving ancestors
  interpol = deepcopy(roots)
  interpol[(circshift(roots, -1) .== 1) .& (circshift(roots, 1) .== 1)] .= 1

  for trial in mask(keys(data), "trial")
    pids = unique(data[trial]["pIDs"])
    cnts += 1
    matches = 0

    for id in pids
      interpol[id] == 1 && (matches += 1)
      survivalFreq[id] += 1
    end
    push!(accuracy, matches / length(pids))
  end

  # fraction captured
  auc = sum(interpol .* survivalFreq) / sum(survivalFreq)
  acc = measurement(mean_and_std(accuracy)...)
  return (acc_ancestors=acc, freq_ancestors=survivalFreq ./ cnts, auc=auc)
end

function accImg!(ax, results)
  for (i, dms) in dataCasks
    # hotspots -> agreement with geoemtric paths
    x, xerr = unzip(getfield.(dms, :parameter))
    y, yerr = unzip(getfield.(dms, :acc_hotspots))

    # draw accuracy graph for hotspots
    scatter!(ax, x, y; markersize=35, color=get(ColorSchemes.roma, 1), label="Hotspots (M)")
    errorbars!(ax, x, y, xerr; color=get(ColorSchemes.roma, 1), direction=:x)
    errorbars!(ax, x, y, yerr; color=get(ColorSchemes.roma, 1))

    # estimates of background for hotspot measure
    bg, bgerr = unzip(getfield.(dms, :acc_hotspots_bg))

    lines!(ax, x, fill(value(first(bg)), length(x)); color=:red)

    if i == 1
      kwargs = (fontsize=23, rotation=0.0)
      text!(ax, "Random Hotspot Selection"; position=(2.5, 0.34), align=(:left, :center), kwargs...)
    end

    # ancestors -> agreement with geoemtric paths
    y, yerr = unzip(getfield.(dms, :acc_ancestors))

    # draw accuracy graph for ancestors
    _x, _y = mask(x, y; xlbound=first(x))
    scatter!(ax, _x, _y; markersize=35, color=get(ColorSchemes.roma, 0), label="Ancestors (K)")
    errorbars!(ax, _x, _y, xerr[2:end]; color=:black, direction=:x)
    errorbars!(ax, _x, _y, yerr[2:end]; color=get(ColorSchemes.roma, 0))
  end
  # axislegend(ax, position=:rb, labelsize=25)
  return axislegend(ax; position=:rb, labelsize=25, merge=true)
end

# doc: env objects vists by lineages and optimal paths
function objectVisits(path, opts)
  _objs = load(joinpath(path, "objects.jld2"), "objs")

  generators = filter(t -> last(t) <= opts.ref_line, unique!(_objs))
  generators = convert.(NTuple{2,Float64}, generators)

  # .visited objects by lineages 
  lineageVisits = LineageTracing.objectVisits(generators, path, opts.numberTrials, opts.width, opts.height, opts.radius)

  # .visited objects by optimalPaths
  _visited = Int64.(load(joinpath(path, "processed_visitedObjects_optimalPaths.jld2"), "visitedObjectsIndices"))
  optimalVisits = filter(<=(length(generators)), _visited)
  # visitedRoots =filter(>(length(generators)), visited) .- length(generators) .- opts.width

  return PinningResults.Results(; lineageVisits=lineageVisits, optimalVisits=optimalVisits, generators=generators)
end

# doc: surviving ancestors and optimal path roots
function ancestorVisits(path, width)
  optimalAncestors = load(joinpath(path, "processed_surivingAncestors_optimalPaths.jld2"), "survivingAncestors")

  ancestorCount = zeros(Float64, width)
  overlaps = Float64[]

  jldopen(joinpath(path, "data_phylo.jld2"), "r") do file
    for trial in LineageTracing.mask(keys(file), "trial")
      geneticLabels = file[trial]["source_labels"]
      _uniqueLabels = collect(unique(geneticLabels))
      ancestorCount[_uniqueLabels] .+= 1

      _survivingAncestors = zeros(Int64, width)
      _survivingAncestors[_uniqueLabels] .= 1
      push!(overlaps, sum(_survivingAncestors[optimalAncestors]) / sum(_survivingAncestors))
    end
  end

  # s = sum(ancestors.lineageAncestors[ancestors.optimalAncestors])
  # s / sum(ancestors.lineageAncestors)
  return (lineageAncestors=ancestorCount, optimalAncestors=optimalAncestors, overlap=mean_and_std(overlaps))
end

# doc: lineage map with an averaging kernel
function lineageSnapshot(path, opts)
  # .smoothing lineage map with a 7x7 averaging kernel
  results::LineageTracing.Results = LineageTracing.lineageTraces(path, opts; full=false)
  smoothed_lineageMap = copy(results.lineageMap)
  for x in axes(smoothed_lineageMap, 1), y in axes(smoothed_lineageMap, 2)
    _xrange = clamp.((x - 3):(x + 3), 1, opts.width)
    _yrange = clamp.((y - 3):(y + 3), 1, opts.ref_line)
    smoothed_lineageMap[x, y] = sum(results.lineageMap[_xrange, _yrange]) / 49
  end
  return (smoothed_lineageMap ./= mean(smoothed_lineageMap[smoothed_lineageMap .> 0]))
end

# doc: draw the environment with lineages and visited hotspots
function environment(generators, lineageMap, lineageVisits, optimalVisits; showimg=false)
  fig, ax = heatmap(lineageMap; lowclip=:transparent, colorrange=(0.01, 0.7), colormap=ColorSchemes.hot)

  # .add env objects
  scatter!(ax, generators; color=(:green, 0.5), markersize=15)
  # .add env objects highlighted by optimal path visits
  scatter!(ax, generators[optimalVisits]; color=(:blue, 0.5), markersize=15)

  limits!(ax, (1, opts.width), (1, opts.ref_line))
  hidedecorations!(ax)
  display(fig)

  lineageVisited = findall(>=(0.05), mean(lineageVisits; dims=2)[:, 1])
  scatter!(ax, generators[lineageVisited]; color=(:pink, 0.5), markersize=10)
  showimg ? display(fig) : return fig
end

# doc: metric -> |l ∩ s| ÷ |l ∪ s| == lis / lus
function overlapFraction(lineageVisits, optimalVisits)
  overlaps = Float64[]
  for j in axes(lineageVisits, 2)
    slice = lineageVisits[:, j]
    lis = sum(slice[optimalVisits])
    slice[optimalVisits] .= 1
    lus = sum(slice)
    push!(overlaps, lis / lus)
  end
  return mean_and_std(overlaps)
end

# doc: main process
function process(gdf, cmap, vrange, colorby)
  results = PinningResults.Averaging[]

  # .process subframes grouped by, e.g., intensity
  for subframe in DataFrameTools.convert(gdf, colorby)
    _paths = subframe.path
    _value = mean(subframe[!, colorby])
    ref_line = first(subframe.ref_line)
    width = first(subframe.width)
    height = first(subframe.height)
    radius = first(subframe.radius)
    numberTrials = first(subframe.numberTrials)

    _result = PinningResults.Averaging()

    _tmp_1 = zeros(Measurement{Float64}, length(_paths))
    _tmp_2 = zeros(Measurement{Float64}, length(_paths))

    @threads for j in eachindex(_paths)
      opts = (radius=radius, numberTrials=numberTrials, ref_line=ref_line, width=width, height=height)
      _visits = objectVisits(_paths[j], opts)
      _overlap = measurement(overlapFraction(_visits.lineageVisits, _visits.optimalVisits)...)
      _tmp_1[j] = _overlap
      _overlap = measurement(ancestorVisits(_paths[j], width).overlap...)
      _tmp_2[j] = _overlap
    end

    _result.objectOverlap = mean(_tmp_1)
    _result.ancestorOverlap = mean(_tmp_2)

    _result.color = get(cmap, _value, vrange)
    _result.parameters = (
      density=mean(subframe.density), intensity=mean(subframe.intensity), radius=mean(subframe.radius)
    )
    push!(results, _result)
  end
  return results
end

function cacheSync(cachePath, toplevel, df, _index, _values, opts; rewrite=false)
  mkpath(joinpath(cachePath, toplevel))
  _path = joinpath(cachePath, toplevel, "cache_pinningTransition_$(_index)_$(opts...).jld2")
  rewrite || (ispath(_path) && return namedtuple(load(_path)))

  gdf = DataFrameTools.groupframe(df, _values, _index; opts=opts)
  _bounds = DataFrameTools.bounds(gdf, last(opts))
  results = process(gdf, cmap_1, _bounds, last(opts))
  jldsave(_path; opts=opts, values=_values, index=_index, bounds=_bounds, results=results)
  return (results=results, bounds=_bounds, opts=opts, index=_index, values=_values)
end

function saveImage(fig, path, imgPath, name)
  imgFolderPath = last(splitdir(path))
  img_path = joinpath(imgPath, imgFolderPath)
  figPath = joinpath(img_path, name)
  mkpath(img_path)
  save(figPath, fig)
  return figPath
end

#. set path to be the csv file generated from src/processing/generate_parameter_space_table.py for the data set that you would like to analyze, for example
path = "/home/Projects/disorderedLandscapes/workspace/experiments/2024_25_02/logs/PS_rf:1000_H:1100_W:2000_R:10_G:0.csv"
alt_path = "/storage/disorderedLandscapes/simulations/pathOptimization/"
img_path = "/home/Images/disorderedLandscapes/pinningTransition"

toplevel = first(filter(t -> contains(t, "2024"), splitpath(path)))
cachePath = "/home/Projects/disorderedLandscapes/workspace/experiments/cache/pinningTransition/$(toplevel)_alt"

# >compare hotspots visited by lineages and optimal paths
opts = namedtuple(load(joinpath(path, "Opts.jld2")))
results = objectVisits(path, opts)
lineageMap = lineageSnapshot(path, opts)
environment(results.generators, lineageMap, results.lineageVisits, results.optimalVisits; showimg=true)

objOverlaps = overlapFraction(results.lineageVisits, results.optimalVisits)
ancestors = ancestorVisits(path, opts)

# .defined dataframe and parameters used in search
df = CSV.read(path, DataFrame)
opts = (:density, :intensity)
# opts = (:intensity, :density)
_values = sort(unique(df[!, DataFrameTools._replace(first(opts), :density, :group)]))

cmap_1 = alphaColor(ColorSchemes.brg, 0.5)
customTheme!(28)
relabel = Dict(:objectOverlap => "M (Hotspots)", :ancestorOverlap => "K (Ancestors)")

# .cache results
caches = []
let _index = 3
  # .iterate through parameter space slice, load from cache
  cache = cacheSync(cachePath, toplevel, df, _index, _values, opts)
  push!(caches, cache)
end

# for _index in eachindex(_values)
let
  # >draw transition in visited hotspots
  fig = Figure(; size=(800, 650))

  cache = last(caches)
  _val = round(getfield(first(cache.results).parameters, first(opts)); digits=2)
  # titleKwargs = (title="$(first(opts)): $(round(_val, digits=2))", titlealign=:left)

  kwargs = (xlabelsize=45, ylabelsize=35, yticklabelsize=25, xticklabelsize=25)
  ax = Axis(fig[1, 1]; ylabel="Overlap Fraction", xlabel=L"I", kwargs...)

  # .iterate through results to generate plot
  for result in Iterators.flatten(getfield.(caches, :results))
    for (color, measure) in [(:firebrick, :ancestorOverlap), (:steelblue3, :objectOverlap)]
      config = (color=color, markersize=20, label=relabel[measure])
      x = getfield(result.parameters, last(opts))
      y = value(getfield(result, measure))
      y_err = uncertainty(getfield(result, measure))

      # .skip undefined ancestor measure value at v = 0
      if (measure == :ancestorOverlap && result.parameters.intensity == 0)
        # hlines!(ax, value(result.objectOverlap), color=:red)
      end
      pointPlot!(ax, x, y, 0, y_err, config)
    end
  end

  axislegend(ax; merge=true, position=:lt)
  ylims!(ax, 0, 1.0)
  display(fig)
  fig_name = "$(first(opts)):$(_val).png"
  figPath = saveImage(fig, path, img_path, fig_name)
end
