using JLD2, FileIO, Revise, Parameters, CSV, DataFrames, Base.Threads
include("common/theme.jl")
include("common/binning.jl")
include("common/modeling.jl")
includet("../src/calculations/lineages.jl")

module SurvivalResults
using Measurements
@kwdef mutable struct Result
  survivingAncestors::Union{Missing,Measurement{Float64}} = missing
  color::Union{Missing,Any} = missing
  parameters::Union{Missing,NamedTuple} = missing
end
end

function saveImage(fig, path, imgPath, name)
  imgFolderPath = last(splitdir(path))
  img_path = joinpath(imgPath, imgFolderPath)
  figPath = joinpath(img_path, name)
  mkpath(img_path)
  save(figPath, fig)
  return figPath
end

function cacheSync(cachePath, df, _index, _values, opts; rewrite=false)
  _path = joinpath(cachePath, "cache_pinningTransition_$(_index)_$(opts...).jld2")
  rewrite || (ispath(_path) && return namedtuple(load(_path)))

  gdf = DataFrameTools.groupframe(df, _values, _index; opts=opts)
  _bounds = DataFrameTools.bounds(gdf, last(opts))
  results = process(gdf, cmap_1, _bounds, last(opts))

  mkpath(first(splitdir(_path)))
  jldsave(_path; opts=opts, values=_values, index=_index, bounds=_bounds, results=results)
  return (results=results, bounds=_bounds, opts=opts, index=_index, values=_values)
end

# doc: main process
function process(gdf, cmap, vrange, colorby)
  results = SurvivalResults.Result[]

  # .process subframes grouped by, e.g., intensity
  for subframe in DataFrameTools.convert(gdf, colorby)
    _paths = subframe.path
    _value = mean(subframe[!, colorby])
    opts = copy(subframe[1, :])
    (; ref_line, width, height, radius, numberTrials) = copy(subframe[1, :])

    result = SurvivalResults.Result()

    survival = zeros(Measurement{Float64}, length(_paths))

    @threads for j in eachindex(_paths)
      survival[j] = survivingAncestors(_paths[j])
    end

    result.survivingAncestors = mean(survival)
    result.color = get(cmap, _value, vrange)
    result.parameters = (
      density=mean(subframe.density), intensity=mean(subframe.intensity), radius=mean(subframe.radius)
    )
    push!(results, result)
  end
  return results
end

function survivingAncestors(path)::Measurement{Float64}
  ancestorCount = Float64[]

  jldopen(joinpath(path, "data_phylo.jld2"), "r") do file
    for trial in FilterTools.mask(keys(file), "trial")
      geneticLabels = file[trial]["source_labels"]
      push!(ancestorCount, length(unique(geneticLabels)))
    end
  end
  return measurement(mean_and_std(ancestorCount)...)
end

paths = [
  "/home/Projects/disorderedLandscapes/workspace/experiments/2024_25_02/logs/PS_rf:1000_H:1100_W:2000_R:10_G:0.csv",
  "/home/Projects/disorderedLandscapes/workspace/experiments/2024_26_02/logs/PS_rf:1000_H:1100_W:2000_R:10_G:0.csv"
]
alt_path = "/storage/disorderedLandscapes/simulations/pathOptimization/"
toplevel = first(filter(t -> contains(t, "2024"), splitpath(first(paths))))
img_path = "/home/Images/disorderedLandscapes/survival/$(toplevel)"
cachePath = "/home/Projects/disorderedLandscapes/workspace/experiments/cache/survival/$(toplevel)"

# .theme options
customTheme!(28)
cmap_1 = alphaColor(ColorSchemes.brg, 0.95)
getColor(v, cRange, cmap) = get(cmap, v, cRange)

# .storage for cached results
cached = Dict()

# .iterate through all paths to get the msd values and cache them
opts = (:density, :intensity)
# opts = (:intensity, :density)

for path in paths
  df = CSV.read(path, DataFrame)
  values = sort(unique(df[!, DataFrameTools._replace(first(opts), :density, :group)]))

  # .process through parameter space slice
  for _index in eachindex(values)
    gdf = DataFrameTools.groupframe(df, values, _index; opts=opts)
    _bounds = DataFrameTools.bounds(gdf, last(opts))
    results = process(gdf, cmap_1, _bounds, last(opts))
    _val = round(getfield(first(results).parameters, first(opts)); digits=2)
    print(_index, " ")

    # .add result to cache and find colorrange
    if haskey(cached, _val)
      append!(cached[_val], results)
    else
      cached[_val] = results
    end
  end
end

colorrange = extrema(collect(keys(cached)))

# >draw transition in visited hotspots
let
  fig = Figure(; size=(800, 750))

  # xscale=log10, yscale=log10)
  kwargs = (xlabelsize=38, ylabelsize=33, yticklabelsize=25, xticklabelsize=25)
  ax = Axis(fig[1, 1]; ylabel=L"N / N^{\text{uniform}}", xlabel=L"I", kwargs...)
  inset = Axis(
    fig;
    ylabel=L"N / N^{\text{uniform}}",
    xlabel=L"\Lambda^{-\frac{1}{3}} \frac{I}{I + 1}",
    # xlabel=L"\left(\frac{1}{\Lambda} \right) ^\frac{-1}{3} \frac{I}{I + 1}",
    # xlabel=L"\sqrt{\phi} \frac{I}{I + 1}",
    xtrimspine=false,
    ytrimspine=(false, false),
    xlabelsize=30,
    ylabelsize=35,
    xticklabelsize=28,
    yticklabelsize=28,
    bbox=BBox(rect(420, 450, 370, 0.8)...)
  )

  for opt_val in sort(collect(keys(cached)))
    # .iterate through results to generate plot
    cache = cached[opt_val]
    sorted = sort(cache; by=t -> getfield(t.parameters, :intensity))
    density = getfield.(getfield.(sorted, :parameters), :density)
    colors = getColor(density, colorrange, cmap_1)

    density[1] <= 0.45 || continue
    # 0.06 <= density[1] < 0.45 || continue

    _measurements = getfield.(sorted, :survivingAncestors)
    intensities = getfield.(getfield.(sorted, :parameters), :intensity)

    # .normalize by first element (uniform case)
    y = value.(_measurements ./ first(_measurements))
    y_err = uncertainty.(_measurements ./ first(_measurements))

    #? should this be in terms of the density or hotspot separation?
    r = getfield.(getfield.(sorted, :parameters), :radius)
    λ = @. √(-π * r^2 / log(1 - density))

    # _y = @. y 
    # _x = @. (intensities  / (1 + intensities))  / λ ^ (1/4)
    _y = @. y
    _x = @. 100 * (1 / (2 * r + λ))^(1 / 3) * (intensities / (1 + intensities))  / (10^1.5)
    config(i, d) = (color=colors[i], markersize=15, label="$(round(d[i], digits=2))")

    for i in eachindex(y)
      # i > 2 || continue
      pointPlot!(ax, intensities[i], y[i], 0, y_err[i], config(i, density))
      pointPlot!(inset, _x[i], _y[i], 0, y_err[i], config(i, density))
    end
  end

  axislegend(ax, ax, L"\phi"; merge=true, position=:lb, titlesize=40, orientation=:horizontal, nbanks=3, fontsize=20)

  ylims!(ax, 10^(-0.6), 1)
  # ax.xscale = log10
  ax.yscale = log10

  # inset.xscale = log10
  inset.yscale = log10
  display(fig)

  fig_name = "fig_5_C.pdf"
  figPath = saveImage(fig, first(paths), img_path, fig_name)
end

let 
  errors = []
  for c in 1:0.01:10
    X, Y = Float64[], Float64[]
    for opt_val in sort(collect(keys(cached)))
      cache = cached[opt_val]
      sorted = sort(cache; by=t -> getfield(t.parameters, :intensity))
      density = getfield.(getfield.(sorted, :parameters), :density)
      colors = getColor(density, colorrange, cmap_1)

      density[1] <= 0.45 || continue

      _measurements = getfield.(sorted, :survivingAncestors)
      intensities = getfield.(getfield.(sorted, :parameters), :intensity)

      # .normalize by first element (uniform case)
      y = value.(_measurements ./ first(_measurements))
      y_err = uncertainty.(_measurements ./ first(_measurements))
      r = getfield.(getfield.(sorted, :parameters), :radius)
      λ = @. √(-π * r^2 / log(1 - density))
      _y = @. y
      _x = @. 100 * (2000 / (2 * r + λ))^(1 / c) * (intensities / (1 + intensities))  / (10^1.5)
      # _x = @. (intensities  / (1 + intensities))  * density ^(1/c)
      push!(X, _x[2:end]...)
      push!(Y, _y[2:end]...)
    end

    fit = FitModels.linearfit(X, log.(Y))
    push!(errors, fit.sumSqRes)
  end

  
  fig = Figure(; size=(800, 750))
  # println(collect(1:0.01:10)[argmin(errors)])
  x_min = minimum(errors)
  println(errors[argmin(errors)-1], errors[argmin(errors)], errors[argmin(errors)+1])
  right = findlast(errors .<= x_min*1.02)
  left = findfirst(errors .<= x_min*1.02)
  

  kwargs = (xlabelsize=38, ylabelsize=33, yticklabelsize=25, xticklabelsize=25)
  ax = Axis(fig[1, 1]; ylabel="Sum Square Residuals", xlabel=L"\eta", kwargs...)
  scatter!(ax, collect(1:0.01:10), Float64.(errors))
  display(fig)
  figPath = saveImage(fig, first(paths), img_path, "smsqres.pdf")
end
