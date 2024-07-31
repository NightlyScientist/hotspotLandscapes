using JLD2, FileIO, Revise
using CSV, DataFrames, StatsBase, Base.Threads
includet("common/theme.jl")
includet("common/binning.jl")
includet("common/modeling.jl")
includet("../src/calculations/lineages.jl")
include("../src/base/environment.jl")
includet("../src/calculations/OptimalPaths.jl")
include("common/theme.jl")
include("../src/calculations/lineages.jl")

function saveImage(fig, path, imgPath, name)
  imgFolderPath = last(splitdir(path))
  img_path = joinpath(imgPath, imgFolderPath)
  figPath = joinpath(img_path, name)
  mkpath(img_path)
  save(figPath, fig)
  return figPath
end

# doc: add heat map of surviving ancestors below lineage map
function addSurvivalFreqImg!(fig, result, cmap, cli)
  sub_axis = Axis(fig[3, 1]; backgroundcolor=:white, xlabelsize=24, ylabelsize=24)

  freq_ancestors = result.ancestorCounts ./ cli.numberTrials
  # freq_ancestors = result.ancestorCounts ./ sum(result.ancestorCounts)
  img = reshape(repeat(freq_ancestors, 2), (length(freq_ancestors), 2))
  crange = (0.005, maximum(img))
  heatmap!(sub_axis, img; lowclip=:white, colormap=cmap, colorrange=crange)

  hidedecorations!(sub_axis; grid=true)

  kwargs = (ticksize=5, ticklabelsize=20, vertical=false, flipaxis=false, labelsize=24, label="Survival Probability")
  Colorbar(fig[4, 1]; colormap=cmap, colorrange=crange, kwargs..., lowclip=:white)

  xlims!(sub_axis, (1, length(freq_ancestors)))
  rowsize!(fig.layout, 3, Fixed(25))
  rowsize!(fig.layout, 4, Fixed(5))
  rowgap!(fig.layout, 1)
  return nothing
end

# doc: add main heat map of lineage positions
function addLineagemapImg!(fig, result, cmap, opts, data_path)
  ax = Axis(fig[2, 1]; backgroundcolor=:white)

  # task: normalize to something better
  hm = result.lineageMap ./ maximum(result.lineageMap)
  colorrange = (0.001, 0.5 * maximum(hm))
  h = heatmap!(ax, hm; lowclip=:white, colormap=cmap, colorrange=colorrange)

  kwargs = (ticksize=10, ticklabelsize=20, tickalign=1, labelsize=24, flipaxis=true, vertical=false)
  Colorbar(fig[1, 1], h; label="Lineage Visitation Frequency", kwargs...)

  colors = alphaColor(ColorSchemes.grays, 0.3; ncolors=3)
  if opts.intensity > 0
    objs = load(joinpath(data_path, "objects.jld2"), "objs")
    env = first(applyObstacles!(objs, opts.radius, opts.width, opts.height))
    env = reshape(env, opts.dims)
    heatmap!(ax, env; colorrange=(2, 3), lowclip=:transparent, colormap=colors)
  end

  rowgap!(fig.layout, 1)
  colgap!(fig.layout, 1)
  rowsize!(fig.layout, 2, Fixed(300))
  rowsize!(fig.layout, 1, Fixed(10))
  hidedecorations!(ax)
  scalebar!(ax, (150, 50); width=200, height=3, offset=50, color=:black)
  top_line = opts.ref_line == 0 ? opts.height : opts.ref_line
  limits!(ax, 1, opts.width, 1, top_line)
  return ax
end

# >collect data files
# task: convert to dataframes
data_path = ""

_dataFileName = "data_phylo.jld2"
_fieldNames = ["end_time", "source_labels", "phylogeny", "branch_points", "branch_times"]

opts = namedtuple(FileIO.load(joinpath(data_path, "Opts.jld2")))
results::LineageTracing.Results = LineageTracing.lineageTraces(data_path, opts)

customTheme!(22)
cmap = alphaColor(ColorSchemes.dense, 1.0)
cmap_2 = alphaColor(ColorSchemes.tempo, 1.0)

# >pinning heatmap and ancestor Survival Frequencies
fig = Figure(; backgroundcolor=:white, size=(600, 600 * 0.87))

# .lineage heat map of Visitation Frequency
addLineagemapImg!(fig, results, cmap, opts)

# .surviving ancestors heatmap
addSurvivalFreqImg!(fig, results, cmap_2)
display(fig)

# >survival maps sequenced by intensity
imgPath = "/home/Projects/disorderedLandscapes/workspace/images/survivalSequence"

paths = [""]

opts = (:density, :intensity)

df = CSV.read(first(paths), DataFrame)

intensities = df.intensity[sortperm(df.intensity)]
sortedPaths = df.path[sortperm(df.intensity)]
uniqueValuesIndices = Int64[]

for i in eachindex(intensities)
  if i == lastindex(intensities)
    push!(uniqueValuesIndices, i)
  elseif intensities[i] != intensities[i + 1]
    push!(uniqueValuesIndices, i)
  end
end

uniqueValuesIndices
intensities[uniqueValuesIndices]
paths = sortedPaths[uniqueValuesIndices]

colorGrid = zeros(Float64, (2000, length(uniqueValuesIndices)))

@threads for i in eachindex(paths)
  cli = namedtuple(FileIO.load(joinpath(paths[i], "Opts.jld2")))
  result::LineageTracing.Results = LineageTracing.lineageTraces(paths[i], cli)

  freq_ancestors = result.ancestorCounts ./ cli.numberTrials
  # freq_ancestors = result.ancestorCounts ./ sum(result.ancestorCounts)
  colorGrid[:, i] .= freq_ancestors
end

# .heatmap sequence of survival probabilities
let
  crange = (0, 0.2 * maximum(colorGrid))
  fig = Figure(; size=(800, 600))
  ax = Axis(
    fig[1, 1];
    xlabel=L"x",
    ylabel=L"I",
    xlabelsize=25,
    ylabelsize=25,
    xticklabelsize=16,
    yticklabelsize=16,
    xaxisposition=:top,
    ylabelrotation=pi / 2
  )

  heatmap!(ax, colorGrid; colorrange=crange, colormap=ColorSchemes.hot)
  # hidedecorations!(ax; grid=true)

  # .set y-ticks to reflect intensity
  ax.yticks = (collect(1:2:length(uniqueValuesIndices)), string.(intensities[uniqueValuesIndices])[1:2:end])

  # .color bar to show survival Probability
  kwargs = (ticksize=5, ticklabelsize=20, vertical=false, flipaxis=false, labelsize=21, label="Survival Probability")

  ticks = ([0, crange[2]], ["0", string(last(crange))])
  Colorbar(fig[2, 1]; colormap=ColorSchemes.hot, kwargs..., colorrange=crange)
  # Colorbar(fig[2, 1]; colormap=cmap, ticks=([0, 0.5, 1], ["0", "0.5", "1"]), kwargs...)
  rowgap!(fig.layout, 1)

  save(imgPath * "/survialProbability.png", fig)
  display(fig)
end

# .sample landscape images
cmap = alphaColor(ColorSchemes.dense, 1.0)
cmap_2 = alphaColor(ColorSchemes.amp, 1.0)

# let
level = 20
# level = 29
data_path = paths[level]

opts = namedtuple(FileIO.load(joinpath(data_path, "Opts.jld2")))
result::LineageTracing.Results = LineageTracing.lineageTraces(data_path, opts)

if ispath(joinpath(data_path, "objects.jld2"))
  objs = load(joinpath(data_path, "objects.jld2"), "objs")
else
  objs = load(joinpath(data_path, "data.jld2"), "htspts")
end

htspts = filter(t -> last(t) <= opts.ref_line, unique!(objs))
points = transpose(hcat(Float64.(first.(htspts)), Float64.(last.(htspts))))

voronoiGraph = ContinuousOptPaths.voronoiTriangulation(points, opts)

# .nearest hotspot to each source (sink) node
sinks_to_htpts, sources_to_htspts = ContinuousOptPaths.sourceSinkConnections(
  voronoiGraph, opts; ref_line=opts.ref_line, n_near=30, interval=2
)

# .combine sources, sinks, and hotspots as nodes in the graph
g = ContinuousOptPaths.weightedGraph(voronoiGraph, sources_to_htspts, sinks_to_htpts, opts; NN=false)

sourceIndices, sinkIndices = ContinuousOptPaths.edgeNodeIndices(
  sources_to_htspts, sinks_to_htpts, voronoiGraph.generators
)

gPositions = ContinuousOptPaths.graphPosisions(sources_to_htspts, sinks_to_htpts, voronoiGraph.generators, opts)

optimalPathSets = @time ContinuousOptPaths.allOptimalPaths(
  sourceIndices, sinkIndices, opts, gPositions, g; n=200, cutoff=1.01, floyd=false
)

fig = Figure(; backgroundcolor=:white, size=(600, 600 * 0.87))

# .lineage heat map of Visitation Frequency
ax = addLineagemapImg!(fig, result, cmap, opts, data_path)

# .surviving ancestors heatmap
addSurvivalFreqImg!(fig, result, cmap_2, opts)
# display(fig)

#. NN network paths
# task: hash table of NTuples of positions to fix horizontal lines
segments = Dict{NTuple{4,Float64},Int64}()

for source_index in sourceIndices
  for optimalPath in optimalPathSets[source_index]
    for j in 2:2:length(optimalPath.pathPositions)
      (x1, y1) = optimalPath.pathPositions[j - 1]
      (x2, y2) = optimalPath.pathPositions[j]
      if ~haskey(segments, (x1, y1, x2, y2))
        segments[(x1, y1, x2, y2)] = 1
      end
    end
  end
end

for (key, _) in segments
  x1, y1, x2, y2 = key
  lines!(ax, [x1, x2], [y1, y2]; color=(:brown, 0.5), linewidth=2)
  # lines!(ax, [], yvecs[interval]; color=(:brown, 0.5), linewidth=1)
end

display(fig)
save(imgPath * "/lineagemap_$(opts.intensity).png", fig)
# end
