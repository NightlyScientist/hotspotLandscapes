using JLD2, FileIO, Base.Threads, Parameters, Revise, CSV, DataFrames
include("../src/calculations/lineages.jl")
include("../src/base/trees.jl")
include("common/theme.jl")
include("common/binning.jl")
include("common/modeling.jl")
include("routines/ancestry.jl")

module Tools
using Measurements
@kwdef mutable struct Result
  x::Union{Missing,Vector{Float64}} = missing
  y::Union{Missing,Vector{Float64}} = missing
  T2::Union{Missing,Measurement{Float64}} = missing
  color::Union{Missing,Any} = missing
  parameters::Union{Missing,NamedTuple} = missing
  slope::Union{Missing,Measurement} = missing
  intercept::Union{Missing,Measurement} = missing
end
end

function CommonAncestryTime(coalRatetx, times, max_dist)
  t2x = NTuple{2,Float64}[]
  for i in eachindex(coalRatetx)
    i <= max_dist || continue
    sum(coalRatetx[i, :]) > 1 || continue
    push!(t2x, (i, sum((coalRatetx[i, :] .* times)) / sum(coalRatetx[i, :])))
  end
  return getfield.(t2x, 1), getfield.(t2x, 2)
end

function tauHist(results, max_dist)
  _x = Float64[]

  for dx in axes(results.jtx, 1)
    (Float64(dx) > maxdist || sum(results.jtx[dx, :]) == 0) && continue
    _x = vcat(_x, (results.jtx[dx, :] .* times)[results.jtx[dx, :] .> 0])
  end

  p = sortperm(vcat(_x, results.bins))
  bins = findall(>(length(_x)), p)
  bin_count = Int.(diff(bins) .- 1)
  return bin_count, times
end

function saveImage(fig, path, imgPath, name)
  imgFolderPath = last(splitdir(path))
  img_path = joinpath(imgPath, imgFolderPath)
  figPath = joinpath(img_path, name)
  mkpath(img_path)
  save(figPath, fig)
  return figPath
end

# doc: main process
function process(gdf, cmap, vrange, colorby, max_distance, scale_factor)
  results = Tools.Result[]

  # .process subframes grouped by, e.g., intensity
  for subframe in DataFrameTools.convert(gdf, colorby)
    _paths = subframe.path
    _value = mean(subframe[!, colorby])

    fits = Vector{FitModels.LinearFitModel}(undef, length(_paths))
    rescaledXs = Vector{Vector{Float64}}(undef, length(_paths))
    rescaledYs = Vector{Vector{Float64}}(undef, length(_paths))
    Ns = Vector{Float64}(undef, length(_paths))

    @threads for j in eachindex(_paths)
      file = jldopen(joinpath(_paths[j], "processed_ancestry.jld2"), "r")
      result::Ancestry.AncestryTools.Results = file["result"]
      nbins = result.nbins

      # .rescale coalescence rate
      rescaledX, rescaledY = Ancestry.rescaleCoalescenceRate(result, max_distance; nbins=nbins, scale=scale_factor)
      rescaledXs[j] = rescaledX
      rescaledYs[j] = rescaledY

      # .find power law for τ/x >> 1 
      fits[j] = FitModels.linearfit(rescaledX, rescaledY; xlbound=1.5 * 10^-1, xhbound=8.0 * 10^0, scale=log10)
      # fits[j] = FitModels.linearfit(rescaledX, rescaledY; xlbound=2.5 * 10^0, xhbound=9.0 * 10^1, scale=log10)

      # .mean common ancestry time
      distances, t2xt = CommonAncestryTime(result.jtx, result.bins, 55)
      mft = last(result.nbins * result.binSize)

      # .set normalization to null case
      Ns[j] = mean(t2xt) / mft
      # y = @. measurement(mean_and_std(t2xt)...) / mft / N₀
    end

    _x = vcat(rescaledXs...)
    _y = vcat(rescaledYs...)
    binEdges = collect(logspace(log.(2, extrema(_x))..., 1000; base=2.0))
    x_binned, y_binned, _ = binning(_x, _y, binEdges)

    _result::Tools.Result = Tools.Result(; x=x_binned, y=y_binned[:, 1])

    _result.T2 = measurement(mean_and_std(Ns)...)
    _result.slope = mean(getfield.(fits, :a) .± getfield.(fits, :a_err))
    _result.intercept = mean(getfield.(fits, :b) .± getfield.(fits, :b_err))
    _result.color = get(cmap, _value, vrange)
    _result.parameters = (
      density=mean(subframe.density), intensity=mean(subframe.intensity), radius=mean(subframe.radius)
    )
    push!(results, _result)
  end
  return results
end

# ******************************
path = ""
img_path = "/home/Images/disorderedLandscapes/ancestry/processed/"

# >global constructs
max_distances = (250, 10)
scale_factor = 1.5

df = CSV.read(path, DataFrame)
# opts = (:density, :intensity)
opts = (:intensity, :density)

values = sort(unique(df[!, DataFrameTools._replace(first(opts), :density, :group)]))
cmap_1 = alphaColor(ColorSchemes.roma, 0.5)
labels = Dict(:intensity => L"I", :density => L"\phi")

# .constant appearing in front of J for null case
constant = let
  constants = Float64[]
  gdf = DataFrameTools.groupframe(df, values, 1; opts=opts)
  _bounds = DataFrameTools.bounds(gdf, last(opts))
  results = process(gdf, cmap_1, _bounds, last(opts), max_distances[1], scale_factor)

  foreach(t -> push!(constants, (abs ∘ value)(t.intercept)), results)
  mean(constants)
end

all_results = Dict{Int,Vector{Tools.Result}}()
second = copy(all_results)

# for _index in 10:11
# for _index in eachindex(values)
let _index = 3
  gdf = DataFrameTools.groupframe(df, values, _index; opts=opts)
  _bounds = DataFrameTools.bounds(gdf, last(opts))

  #. cache results in dict
  if ~haskey(all_results, _index)
    results = process(gdf, cmap_1, _bounds, last(opts), max_distances[1], scale_factor)
    all_results[_index] = results
  end

  #> main figure of coalescence rate
  customTheme!(22)
  xlabel_1 = L"\tau / \Delta x_0^{3/2}"
  xlabel_2 = L"\beta \tau / \Delta x_0^{3/2}"
  ylabel = L"\Delta x_0 \cdot J(\tau \; | \; \Delta x_0)"
  kwargs = (xlabelsize=35, ylabelsize=35, xscale=log10, yscale=log10)

  fig = Figure(; size=(1000, 400))
  ax = Axis(fig[1, 1]; xlabel=xlabel_1, ylabel=ylabel, kwargs...)
  twin = Axis(fig[1, 2]; xlabel=xlabel_2, ylabel=ylabel, kwargs...)

  for (j, result) in enumerate(all_results[_index])
    maskedX, maskedY = FilterTools.mask(result.x, result.y; xlbound=1e-4, xhbound=1e2, scale=log10)
    scatter!(ax, maskedX, maskedY; markersize=10, color=result.color)

    # .rescale factor, b, appearing in (b τ) / Δx
    rescaling = exp(-result.intercept)

    # .scaled coalescence for τ / x >> 1
    scatter!(twin, maskedX .* value(rescaling), maskedY; markersize=10, color=result.color)

    # maskedX, maskedY = FilterTools.mask(result.x, result.y; xlbound=2.0 * 10^-1, xhbound=5.0 * 10^0, scale=log10)
    # scatter!(ax, maskedX, maskedY; markersize=10, color=:black)
  end

  #. colorbar of all samples in the figure
  colran = extrema(map(t -> getfield(t.parameters, last(opts)), all_results[_index]))
  config = (colorrange=colran,)
  kwargs = (vertical=false, ticksize=10, tickalign=1, flipaxis=true, labelsize=35, label=L"\phi")
  Colorbar(fig; colormap=ColorSchemes.roma, bbox=rect(700, 20, 200, 1), kwargs..., config...)

  colgap!(fig.layout, 1)
  # figPath = saveImage(fig, path, img_path, fig_name)
  # println(figPath)
  display(fig)
end

# >slope values graph (inset)
let _index = 2
  gdf = DataFrameTools.groupframe(df, values, _index; opts=opts)
  _bounds = DataFrameTools.bounds(gdf, last(opts))

  #. cache results in dict
  if ~haskey(all_results, _index)
    results = process(gdf, cmap_1, _bounds, last(opts), max_distances[1], scale_factor)
    all_results[_index] = results
  end

  fig = Figure(; backgroundcolor=:white)

  kwargs = (xlabelsize=50, ylabelsize=100, xticks=WilkinsonTicks(4))
  ax = Axis(fig[1, 1]; yticklabelcolor=:grey, rightspinecolor=:grey, ytickcolor=:grey, kwargs...)

  twin = Axis(fig[1, 1]; yaxisposition=:right, yticklabelcolor=:red, rightspinecolor=:red, ytickcolor=:red, kwargs...)

  betas = Measurement{Float64}[]

  for (j, result) in enumerate(all_results[_index])
    x = getfield(result.parameters, last(opts))
    slope = result.slope
    intercept = result.intercept

    # .plot alpha and uncertainty
    kwargs = (color=:gray, marker=:circle, markersize=20, label=L"\alpha")
    pointPlot!(ax, x, (abs ∘ value)(slope), 0.0, uncertainty(slope), kwargs)

    # .rescale factor, c, appearing in (c τ) / Δx
    _scale = exp(-result.intercept / result.slope) * (constant^(1 / result.slope))
    push!(betas, _scale)

    # .point for β in inset plot
    kwargs = (color=:red, marker=:dtriangle, markersize=20, label=L"\beta")
    pointPlot!(twin, x, value(_scale), 0, uncertainty(_scale), kwargs)
  end

  # .calculate λ(ϕ)
  x = collect(LinRange(_bounds..., 100))
  y = (lambda.(x, 10, lambda_sqr) .+ 20) .^ (4 / 3)
  y = y .* (maximum(betas) / maximum(y))
  lines!(twin, x .+ 0.02, value.(y); color=:blue, linewidth=5, label=L"λ(\phi, R)")

  # axislegend(ax; position=:rt, labelsize=45, linewidth=80, markersize=80)
  ylims!(ax, (1.3, 1.9))
  # ylims!(twin, (0.25, 0.45))
  hidespines!(ax, :t, :r)
  hidespines!(twin, :t, :r)
  display(fig)
end

#> common Ancestry time (T_2)
#. prepare T2 values
allresults = Dict{Int,Vector{Tools.Result}}()
second = copy(allresults)
third = copy(allresults)

for _index in eachindex(values)
  gdf = DataFrameTools.groupframe(df, values, _index; opts=opts)
  _bounds = DataFrameTools.bounds(gdf, last(opts))

  #. cache results in dict
  results = process(gdf, cmap_1, _bounds, last(opts), max_distances[1], scale_factor)
  allresults[_index] = results
end

_p = "/home/jgonzaleznunez/Projects/disorderedLandscapes/workspace/experiments/cache/ancestry/ancestry.jld2"
jldsave(_p; allresults=allresults, second=second, third=third, ts_results=ts_results)

ts_results = Dict{Float64,Vector{Tools.Result}}()
for _set in [allresults, second, third]
  for (k, v) in _set
    intensity = getfield(first(v), :parameters).intensity
    if haskey(ts_results, intensity)
      push!(ts_results[intensity], v...)
    else
      ts_results[intensity] = v
    end
  end
end

intensities = Float64[]
Ns = Float64[]
for (k, v) in ts_results
  is = getfield.(getfield.(v, :parameters), :intensity)
  push!(intensities, is...)
  if mean(is) == 0
    append!(Ns, value.(getfield.(v, :T2)))
  end
end
N₀ = mean(Ns)

# .computer precision for intensities give different values when they should be the same
proper_container = Dict{Float64,Vector{Tools.Result}}()
for (k,v) in ts_results
  _k = round(k, digits=2)
  if haskey(proper_container, _k)
    push!(proper_container[_k], v...)
  else
    proper_container[_k] = v
  end
end

customTheme!(28)
# theme(:white)
set_theme!(theme_latexfonts(); figure_padding=16, fontsize=20)

i_index = 0
let cmap = ColorSchemes.rainbow1
  global i_index += 1
  fig = Figure(; size=(830, 600))
  ax = Axis(fig[1, 1]; xlabel=L"\phi", ylabel=L"T_2/T_{2,0}", xlabelsize=40, ylabelsize=35,
    xticklabelsize=30,
    yticklabelsize=30,
    xticksize=10,
    yticksize=10,
    xticklabelsvisible=true,
    yticklabelsvisible=true,
    topspinevisible = false,
    rightspinevisible = false,
  )
  hidedecorations!(ax, ticklabels = false, ticks=false, minorticks = false, label=false)

  _keys = sort(collect(keys(proper_container)))
  sets = []
  for k in _keys
    k < 1.0 && (round(Int64, k * 10) % 2 == 0 || continue)
    # 1 < k < 2.0 && continue
    v = proper_container[k]
    _set = []

    for (j, result) in enumerate(v)
      (; density, intensity, radius) = result.parameters
      color = get(cmap, intensity, extrema(intensities))
      # if j == 1
      # else
      T2 = result.T2 / N₀
      kwargs = (color=color, marker=:circle, markersize=22, label="I = $(round(intensity, digits=2))")
      # λ = √(-π * radius^2 / log(1 - density))
      pointPlot!(ax, density, value(T2), 0, uncertainty(T2) / √(10), kwargs)
      push!(_set, (color, density, value(T2)))
      # end
    end
    push!(sets, _set)
  end

  # .plot all points
  for _set in sets
    color = getfield.(_set, 1)
    density = getfield.(_set, 2)
    sorter = sortperm(density)

    y = getfield.(_set, 3)
    lines!(ax, density[sorter], y[sorter], color=color[sorter], linewidth=4)
  end

  c_r = extrema(intensities)

  Colorbar(
    fig;
    bbox=BBox(600, 805, 100, 120),
    colormap=cmap,
    colorrange=c_r,
    flipaxis=true,
    label=L"I",
    labelsize=35,
    ticklabelsize=25,
    ticksvisible=true,
    vertical=false,
    xticks=collect(0.0:0.2:1.0),
  )

  # vlines!(ax, 0.54; color=:black)
  # text!(ax, L"\phi = 0.55"; position=(0.55, 0.71), align=(:left, :center), rotation=0, fontsize=40)

  # vlines!(ax, 0.3; color=:black)
  # text!(ax, L"\phi = 0.3"; position=(0.30, 0.86), align=(:left, :center), rotation=0, fontsize=40)
  # ylims!(ax, 0.4, 0.8)
  # axislegend(ax; position=:rt, merge=true, bbox=BBox(850,900,500,600))
  display(fig)
  figPath = saveImage(fig, path, img_path, "T2_$i_index.pdf")
end