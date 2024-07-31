using JLD2, FileIO, Revise, CSV, DataFrames, DataFramesMeta, Parameters, Base.Threads
includet("common/theme.jl")
includet("common/binning.jl")
includet("common/modeling.jl")
includet("../src/calculations/lineages.jl")

module MSDResults
using Measurements
@kwdef mutable struct Results
  x::Union{Missing,Vector{Float64}} = missing
  lineageMSD::Union{Missing,Vector{Measurement{Float64}}} = missing
  optimalPathMSD::Union{Missing,Vector{Float64}} = missing
  color::Union{Missing,Any} = missing
  parameters::Union{Missing,NamedTuple} = missing
  tortuosity::Union{Missing,NTuple{2,Float64}} = missing
end
end

function process(gdf, cmap, vrange, colorby)
  results = MSDResults.Results[]

  # .process subframes grouped by, e.g., intensity
  for subframe in DataFrameTools.convert(gdf, colorby)
    _paths = subframe.path
    _value = mean(subframe[!, colorby])
    ref_line = first(subframe.ref_line)

    _result = MSDResults.Results()
    _result.lineageMSD = zeros(Measurement{Float64}, ref_line)
    _result.x = collect(1:ref_line)

    for j in eachindex(_paths)
      _dataframe = CSV.read(joinpath(_paths[j], "processed_lineageMSD.csv"), DataFrame)
      # _result.lineageMSD .+= accumulate(+, (_dataframe.msd .± _dataframe.std_error))[1:ref_line]
      _result.lineageMSD .+= (_dataframe.msd .± _dataframe.std_error)[1:ref_line]
    end

    _result.optimalPathMSD = zeros(Float64, ref_line)
    for j in eachindex(_paths)
      if ispath(joinpath(_paths[j], "processed_optimalPathMSD.csv"))
        _dataframe = CSV.read(joinpath(_paths[j], "processed_optimalPathMSD.csv"), DataFrame)
        _result.optimalPathMSD .+= _dataframe.msd
      end
    end

    _result.lineageMSD ./= length(_paths)
    _result.optimalPathMSD ./= length(_paths)
    _result.color = get(cmap, _value, vrange)
    _result.parameters = (
      density=mean(subframe.density), intensity=mean(subframe.intensity), radius=mean(subframe.radius)
    )
    push!(results, _result)
  end
  return results
end

function saveImage(fig, path, imgPath, name)
  imgFolderPath = last(splitdir(path))
  img_path = joinpath(imgPath, imgFolderPath)
  figPath = joinpath(img_path, name)
  mkpath(img_path)
  save(figPath, fig)
  return figPath
end

# ******************************
paths = [
  "/home/Projects/disorderedLandscapes/workspace/experiments/2024_25_02/logs/PS_rf:1000_H:1100_W:2000_R:10_G:0.csv",
  "/home/Projects/disorderedLandscapes/workspace/experiments/2024_26_02/logs/PS_rf:1000_H:1100_W:2000_R:10_G:0.csv"
]

toplevel = first(filter(t -> contains(t, "2024"), splitpath(first(paths))))
img_path = "/home/jgonzaleznunez/Images/disorderedLandscapes/msd/$(toplevel)"

# .theme options
customTheme!(28)
cmap_1 = alphaColor(ColorSchemes.brg, 0.7)
cmap_2 = alphaColor(ColorSchemes.brg, 0.9)
labels = Dict(:intensity => L"I", :density => L"\phi")


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

    if haskey(cached, _val)
      append!(cached[_val], results)
    else
      cached[_val] = results
    end
  end
end

# .uniform landscape msd
nullMSD = let 
  x = collect(1:1000)
  y = zeros(Float64, 1000)
  for result in cached[0.0]
    y += result.lineageMSD 
  end
  y /= length(cached[0.0])
end

# .generate figure for specific opt value
let _val = 8.0
  fig = Figure(; size=(700, 600), backgroundcolor=:white)

  kwargs = (xlabelsize=38, ylabelsize=28, backgroundcolor=:white, xscale=log10, yscale=log10)
  titleKwargs = (title="$(first(opts)): $(round(_val, digits=2))", titlealign=:left)
  ax = Axis(fig[1, 1]; xlabel=L"l_y", ylabel="Lateral MSD", kwargs...)
  # titleKwargs...)

  kwargs = (xtrimspine=false, ytrimspine=(false, true), xlabelsize=30, ylabelsize=30)
  # inset_1 = Axis(fig; xlabel=labels[last(opt0)], ylabel="slope", bbox=BBox(rect(200, 380, 180, 0.8)...), kwargs...)

  inset = Axis(
    fig;
    xlabel=L"\phi",
    ylabel=L"\text{MSD}(h=10^3)",
    xtrimspine=false,
    ytrimspine=(false, false),
    xlabelsize=32,
    ylabelsize=25,
    xticklabelsize=28,
    yticklabelsize=28,
    yscale=log10,
    xticks=WilkinsonTicks(4),
    bbox=BBox(rect(210, 430, 190, 0.8)...)
  )

  fits = Tuple{Int64,FitModels.LinearFitModel}[]

  all_results = cached[_val]

  sort!(all_results; by=t -> getfield(t.parameters, last(opts)), rev=true)

  value_range = getfield.(getfield.(all_results, :parameters), last(opts))
  colors = get(cgrad(cmap_1), value_range, :extrema)

  MSD_at_y = NTuple{3,Float64}[]
  parameterVals = Float64[]

  # .add each msd to figure
  for (i, result) in enumerate(all_results)
    label = "$(round(getfield(result.parameters, last(opts)), digits=2))"

    _val_2 = round(getfield(result.parameters, last(opts)); digits=2)
    push!(parameterVals, _val_2)

    #. set color according to intensity value
    color = colors[i]

    # .lineage msd
    x, y = FitModels.mask(result.x, value.(result.lineageMSD); scale=log10, xlbound=1, xhbound=1000)
    lines!(ax, x, y; color=color, linewidth=2)

    # .add msd value as a specific height
    # push!(MSD_at_y, (_val_2, mean_and_std(y[840:842])...))
    push!(MSD_at_y, (_val_2, mean_and_std(y[(end - 3):end])...))
    # push!(MSD_at_y, (_val_2, mean_and_std(y[840:842])...))

    # lines!(ax, x, y; color=color, label=label)
    # lines!(ax, x, y; color=result.color, label=label)
    # lines!(ax, x, x .^ (4/3); color=:black, label = "KPZ (4/3)")

    # .add uniform landscape msd
    if i == 1
      lines!(ax, collect(2:1000), value.(nullMSD)[2:end]; linewidth=5, color=:black, label="Uniform Landscape")
      push!(MSD_at_y, (0, mean_and_std(value.(nullMSD)[(end - 3):end])...))

      # .add vertical line to indicate 1 - 0.54 and 1 - 0.68
      vlines!(inset, 1 - 0.68,  color=:brown)
      addtext!(inset, L"\phi = 0.3", (0.33, 10^3.7), 0.0, 26)
    end

    # .lineage msd masked to a fitting region
    # task: fit regimes to better estimate slope, consider sliding window
    x, y = FitModels.mask(result.x, value.(result.lineageMSD); scale=log10, xlbound=25, xhbound=300)

    # xlbound = i == 1 ? 5 : 20
    # xhbound = i == 1 ? 400 : 400

    fitlimits = (xlbound=10, xhbound=700)

    fit = FitModels.linearfit(result.x, value.(result.lineageMSD); scale=log10, fitlimits...)
    x, y = FitModels.mask(result.x, value.(result.lineageMSD); scale=log10, fitlimits...)
    # lines!(ax, x, y; color=:black, label=label)

    # push!(fits, (i, fit))

    # .optimal path msd 
    # x, y = FitModels.mask(result.x, value.(result.optimalPathMSD) .^ 2; scale=log10, xlbound=1)
    # lines!(ax, x, y; color=result.color, linestyle=:dash)
  end

  # for (i, fit) in fits
  # x = getfield(results[i].parameters, last(opts))
  # pointPlot!(inset_1, x, fit.a, 0, 0, (color=results[i].color,))
  # end

  # .fill inset with msd(h_max) values
  for i in eachindex(MSD_at_y)
    (x, y, y_err) = MSD_at_y[i]
    pointPlot!(inset, x, y, 0, sqrt(y_err), (color=:black,))
  end

  # axislegend(
  # ax, ax; merge=true, unique=true, position=:lt, titlesize=33, orientation=:vertical, nbanks=2
  # ax, ax, labels[last(opts)]; merge=true, unique=true, position=:lt, titlesize=33, orientation=:vertical, nbanks=2
  # )

  kwargs = (flipaxis=true, vertical=false, labelsize=40, colormap=cmap_1, label=L"\phi")
  crange = extrema(parameterVals)
  ticks = round.(LinRange(crange[1], crange[2], 5), digits=2)
  allequal(crange) || @info crange
  Colorbar(fig; bbox=rect(310, 100, 300, 0.1), colorrange=crange, ticks=(ticks, string.(ticks)), kwargs...)

  # if _val > 0 && first(opts) == :intensity
  # if _val > 0 && first(opts) == :intensity
  addline!(ax, 50, 800; pow=4 / 3, shift=0.7, color=:black)
  addtext!(ax, "2α = 4/3", (160, 8000), 0.5)

  addline!(ax, 20, 500; pow=4 / 3, shift=-1, color=:black)
  addtext!(ax, "2α = 4/3", (50, 10), 0.50)
  # elseif first(opts) == :intensity
  #   addline!(ax, 50, 800; pow=4 / 3, shift=-1.2)
  #   addtext!(ax, "α = 4/3", (150, 30), 0.67)
  # end

  # ylims!(ax, 10^-1, 10^5)

  fig_name = "$(first(opts)):$(_val).png"
  figPath = saveImage(fig, paths |> first, img_path, fig_name)
  # ylims!(inset_1, 1.0, 1.5)
  display(fig)
end

_path = df[4, :path]
_nt = df[4, :numberTrials]
_h = df[4, :height]
_w = df[4, :width]
_rf = df[4, :ref_line]

sources = LineageTracing.graphSources(_rf, _w, _h)
MSD = zeros(Float64, (_h, _nt))
_processed = zeros(Bool, _nt)

dataFile = joinpath(_path, "data_phylo.jld2")
jldopen(dataFile, "r") do file
  trials = LineageTracing.mask(keys(file), "trial")

  for j in eachindex(trials)
    phylo = file[trials[j]]["phylogeny"]
    result = LineageTracing.lineageMSD(_h, _w, _rf, sources, phylo; interval=20)
    MSD[:, j] .= result.MSD
    _processed[j] = true

    # .check why we get a huge spike in the msd
    # if result.MSD[6] > 10
    #   x = findall(result.MSD .> 0)
    #   y = result.MSD[x]
    #   f, a = scatter(x, y)

    #   ax.xscale = log10
    #   ax.yscale = log10
    #   display(f)
    #   @info j
    #   break
    # end

    # j == 10 && break
  end
end

index_mask = findall(_processed)
y = mean(MSD[:, index_mask]; dims=2)[1:_rf]

# _dataframe = CSV.read(joinpath(_path, "processed_lineageMSD.csv"), DataFrame)
x = findall(!iszero, y)
f, ax = scatter(x[1:(end - 5)], y[x][1:(end - 5)])
ax.xscale = log10
ax.yscale = log10
display(f)

fit = FitModels.linearfit(x, y[x]; scale=log10, xlbound=10, xhbound=800)
# kwargs = (flipaxis=true, vertical=false, labelsize=40, colormap=cmap_1)
# ticks = round.(LinRange(_bounds[1], _bounds[2], 4), digits=2)
# allequal(_bounds) ||
# Colorbar(fig; bbox=rect(310, 100, 300, 0.1), colorrange=_bounds, ticks=(ticks, string.(ticks)), kwargs...)