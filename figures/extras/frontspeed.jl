using JLD2, FileIO, Revise, CSV, DataFrames, DataFramesMeta, Parameters, Base.Threads
includet("../common/theme.jl")
includet("../common/binning.jl")
includet("../common/modeling.jl")

module SpeedResults
using Measurements
@kwdef mutable struct Results
  frontspeed::Union{Missing,Vector{Measurement{Float64}}} = missing
  color::Union{Missing,Any} = missing
  parameters::Union{Missing,NamedTuple} = missing
end
end

function process(gdf, cmap, vrange, colorby)
  results = SpeedResults.Results[]

  # .process subframes grouped by, e.g., intensity
  for subframe in DataFrameTools.convert(gdf, colorby)
    _paths = subframe.path
    _value = mean(subframe[!, colorby])
    ref_line = first(subframe.ref_line)
    numberSamples = first(subframe.numberSamples)

    _result = SpeedResults.Results()
    _result.frontspeed = zeros(Measurement{Float64}, numberSamples)

    for j in eachindex(_paths)
      jldopen(joinpath(_paths[j], "data_extras.jld2"), "r") do file
        data = file["trial_1"]["front"]
        time = getfield.(data, 1)
        h_t = getfield.(data, 2)
        speed = FitModels.linearfit(time, h_t).a
        _result.frontspeed[j] = speed
      end
      # break
    end

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
paths = String[]

dfs = [CSV.read(path, DataFrame) for path in paths]
df = outerjoin(dfs[1], dfs[2], on=Symbol.(names(dfs[1])), makeunique=false)

opts = (:density, :intensity)
opts = (:intensity, :density)
values = sort(unique(df[!, DataFrameTools._replace(first(opts), :density, :group)]))

customTheme!(28)
# set_theme!(theme_latexfonts(); figure_padding=16, fontsize=20)
cmap_1 = alphaColor(ColorSchemes.brg, 0.7)
cmap_2 = alphaColor(ColorSchemes.brg, 0.9)
labels = Dict(:intensity => L"I", :density => L"\phi")

for _index in eachindex(values)
# let
  # _index = 10
  # _index = 15
  gdf = DataFrameTools.groupframe(df, values, _index; opts=opts)
  _bounds = DataFrameTools.bounds(gdf, last(opts))
  results = process(gdf, cmap_1, _bounds, last(opts))
  _val = round(getfield(first(results).parameters, first(opts)); digits=2)
  println(_val)

  fig = Figure(; size=(700, 600), figure_padding=5, backgroundcolor=:white)
  ax = Axis(
    fig[1, 1]; xlabel=labels[last(opts)], ylabel="Front Speed", xlabelsize=45, ylabelsize=35, xticklabelsize=25, yticklabelsize=25
  )

  for result in results
    _value = getfield(result.parameters, last(opts))
    _frontspeed = first(result.frontspeed)
    _mean = mean(_frontspeed)
    _error = std(_frontspeed)
    scatter!(ax, _value, value(_mean); color="black", markersize=15)
  end

  display(fig)

  p = "Projects/disorderedLandscapes/workspace/images/frontspeed/frontspeed_$(first(opts))_$(_val).svg"
  save(p, fig)
end