#> this file generates images for the sector size scaling in the disordered landscape model
using JLD2, ArgParse, FileIO, Base.Threads, CSV, DataFrames
include("common/theme.jl")
include("common/modeling.jl")

module BoundaryResults
using Measurements
@kwdef mutable struct Results
  x::Union{Missing,Vector{Float64}} = missing
  sectorSizes::Union{Missing,Vector{Float64}} = missing
  color::Union{Missing,Any} = missing
  parameters::Union{Missing,NamedTuple} = missing
end
end

function process(gdf, cmap, vrange, colorby)
  results = BoundaryResults.Results[]

  # .iterate subframes grouped by, e.g., intensity
  for subframe in DataFrameTools.convert(gdf, colorby)
    _paths = subframe.path
    _value = mean(subframe[!, colorby])
    height = first(subframe.height)
    num_trials = first(subframe.numberTrials)

    _result = BoundaryResults.Results()
    _result.sectorSizes = zeros(Float64, height)
    _result.x = collect(1:height)

    for j in eachindex(_paths)
      _result.sectorSizes .+= meanSectorSizes(_paths[j], height, num_trials)
    end

    _result.sectorSizes ./= length(_paths)
    _result.color = get(cmap, _value, vrange)
    _result.parameters = (
      density=mean(subframe.density), intensity=mean(subframe.intensity), radius=mean(subframe.radius)
    )
    push!(results, _result)
  end
  return results
end

function meanSectorSizes(data_path, height, num_trials)
  sectorSize_mean = zeros(Float64, height)

  jldopen(joinpath(data_path, "data_sectors.jld2"), "r") do file
    trials = FilterTools.mask(keys(file), "trial")
    @inbounds for j in eachindex(trials)
      sector_sizes = file[trials[j]]["sector_sizes"]
      sectorSize_mean .+= sector_sizes[:, 1]
    end
  end
  return sectorSize_mean ./ num_trials
end

function saveImage(fig, path, imgPath, name)
  imgFolderPath = last(splitdir(path))
  img_path = joinpath(imgPath, imgFolderPath)
  figPath = joinpath(img_path, name)
  mkpath(img_path)
  save(figPath, fig)
  return figPath
end

# ******************
#. path needs to be set to the csv file generated from src/processing/generate_parameter_space_table.py, for example
path = "workspace/experiments/2024_25_02/logs/PS_rf:1000_H:1100_W:2000_R:10_G:0.csv"
alt_path = "/storage/disorderedLandscapes/simulations/pathOptimization/"

toplevel = first(filter(t -> contains(t, "2024"), splitpath(path)))
img_path = "/home/Images/disorderedLandscapes/boundaries/$(toplevel)"
cachePath = "/home/Projects/disorderedLandscapes/workspace/experiments/cache/boundaries/$(toplevel)"

df = CSV.read(path, DataFrame)
# opts = (:density, :intensity)
opts = (:intensity, :density)
values = sort(unique(df[!, DataFrameTools._replace(first(opts), :density, :group)]))

customTheme!(28)
# cmap_1 = alphaColor(ColorSchemes.brg, 0.5)
cmap_1 = alphaColor(ColorSchemes.roma, 0.9; ncolors=20)
# cmap_2 = alphaColor(ColorSchemes.roma, 0.6; ncolors=20)
labels = Dict(:intensity => L"I", :density => L"\phi")

null_case = []

begin
  _index = 5
  # .process through parameter space slice
  gdf = DataFrameTools.groupframe(df, values, _index; opts=opts)
  _bounds = DataFrameTools.bounds(gdf, last(opts))
  results = process(gdf, cmap_1, _bounds, last(opts))
  _val = round(getfield(first(results).parameters, first(opts)); digits=2)

  fig = Figure(; size=(700, 600), figure_padding=5, backgroundcolor=:transparent)

  if first(opts) == :intensity
    inset_1 = Axis(
      fig;
      xlabel=labels[last(opts)],
      ylabel=L"\Lambda",
      xtrimspine=false,
      ytrimspine=(false, false),
      xlabelsize=32,
      ylabelsize=32,
      xticklabelsize=20,
      yticklabelsize=20,
      bbox=BBox(rect(495, 145, 200, 0.8)...)
    )

    inset_2 = Axis(
      fig;
      xlabel=L"h / h_c",
      ylabel=L"l_s / \Lambda",
      xtrimspine=false,
      ytrimspine=(false, false),
      xlabelsize=32,
      ylabelsize=32,
      xticklabelsize=20,
      yticklabelsize=20,
      xscale=log10,
      yscale=log10,
      bbox=BBox(rect(178, 400, 200, 1.0)...)
    )
  end

  ax = Axis(
    fig[1, 1];
    yscale=log10,
    xscale=log10,
    # xlabel=L"h",
    # ylabel=L"l_s ",
    xlabel=L"h",
    ylabel=L"l_s",
    xlabelsize=40,
    ylabelsize=35,
    xticklabelsize=25,
    yticklabelsize=25,
    # xticks=([0, 10, 10^2, 10^3, 10^4], ["0", "10¹", "10²", "10³", "10⁴"]),
    # yticks=([0, 10, 10^2, 10^3, 10^4], ["0", "10¹", "10²", "10³", "10⁴"])
  )

  _null_case = []

  fits = FitModels.LinearFitModel[]
  for (i, result) in enumerate(results)
    _value = getfield(result.parameters, last(opts))
    println("Processing: ", _value)

    (; density, radius, intensity) = result.parameters

    λ = lambda(density, radius, lambda_sqr) 
    # λ = lambda(density, radius, lambda_sqr) + 2 * radius

    # .draw separation length scale
    :intensity == first(opts) && scatter!(inset_1, _value, λ; color=result.color, markersize=19)

    x, y = FitModels.mask(result.x, result.sectorSizes; scale=log10, xlbound=1, xhbound=500)
    fit = FitModels.linearfit(result.x, result.sectorSizes; scale=log10, xlbound=1, xhbound=500)

    push!(fits, fit)

    # .scaled curves
    _x_ref = result.x[findfirst(result.sectorSizes .> λ)]
    x, y = FitModels.mask(result.x ./ _x_ref, result.sectorSizes ./ λ; scale=log10, xlbound=0)
    scatter!(inset_2, x, y; color=result.color, markersize=5)

    # .unscaled curves
    x, y = FitModels.mask(result.x, result.sectorSizes; scale=log10, xlbound=0)
    scatter!(ax, x, y; color=result.color)
    isapprox(_val, 0.00) && push!(_null_case, (x, y))
  end

  # .add null_case to main ax
  if isempty(null_case)
    null_case = _null_case
  else
    for (x, y) in null_case
      scatter!(ax, x, y; color=:black, markersize=5)
    end
  end

  if :intensity == first(opts)
    if _index == 1
      addtext!(inset_2, "α = 2/3", (0.3, 1.0), 0.7)
      addline!(inset_2, 0.1, 2.3; pow=2 / 3, shift=0.2)
    else
      addtext!(inset_2, "α = 2/3", (0.2, 1.5), 0.66)
      addline!(inset_2, 0.1, 15.5; pow=2 / 3, shift=0.3)
      limits!(ax, 5, 900, 3, 400)
    end
  else
    addtext!(inset_2, "α = 2/3", (120, 80), 0.6)
    addline!(inset_2, 50, 700.0; pow=2 / 3, shift=0.4)
    xlims!(ax, 10, 1200)
    ylims!(ax, 5, 200)
  end

  # .display scene
  display(fig)
  fig_name = "$(first(opts)):$(_val).pdf"
  figPath = saveImage(fig, path, img_path, fig_name)
end