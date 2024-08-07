#> this file generates an animation of the pinning effect
#. input needs to be directory where all subdirectories contain the same landscape at verying intensity values
include("../../src/base/model.jl")
include("../../src/base/containers.jl")
include("../common/modeling.jl")

using JLD2, CairoMakie, ColorSchemes, Colors
using ArgParse
import .Model: gNodes, buildMap

function lineageTracks(dims, data)
  lhmap = zeros(Int64, reduce(*, dims))
  for trial in FilterTools.mask(keys(data), "trial")
    lhmap[collect(keys(data[trial]["phylo"]))] .+= 1
  end
  return reshape(lhmap, dims), maximum(lhmap)
end

function colorscheme_alpha(cscheme::ColorScheme, alpha::T=0.5; ncolors=12) where {T<:Real}
  return ColorScheme([RGBA(get(cscheme, k), alpha) for k in range(0, 1; length=ncolors)])
end

cli = let sts = ArgParseSettings()
  addargs!(sts, name; kwargs...) = add_arg_table!(sts, "--$(name)", Dict(kwargs))

  addargs!(sts, "input"; required=true)
  addargs!(sts, "output"; required=true)
  addargs!(sts, "overwrite"; action=:store_true)
  addargs!(sts, "hotspots"; action=:store_true)
  addargs!(sts, "fps"; default=14, arg_type=Int)
  namedtuple(parse_args(sts))
end
@info "cli: $cli \nsaving recording to: $(cli.output)"

fn = let
  s = splitdir(cli.input)
  last(s) == "" ? last(splitdir(first(s))) : last(s)
end

@info fn
setPath(cli.output * "/" * fn, cli.overwrite, cli)

#. cli.input needs to be a specific env_# directory
#. additionally, the intensity need to be varied in this directory
paths =  readdir(cli.input, join=true)

begin
  #. get env and cli from file, all should be the same so pick one
  env, prms = load(joinpath(paths[1], "Opts.jld2"), "env", "cli")

  # theme and color
  set_theme!(; figure_padding=00)

  # colors for the hotspots
  colors = colorscheme_alpha(ColorSchemes.grays, 0.5; ncolors=3)

  # create figure
  fig = Figure(; resolution=prms.dims .* (√(3), 1.5), backgroundcolor=:transparent)
  ax = Axis(fig[1, 1]; backgroundcolor=:white)

  # set image limits
  hidedecorations!(ax)

  limits!(ax, 1, prms.dims[1], 1, prms.dims[2])

  # total number of trials
  uprbnd = prms.numberTrials

  record(fig, joinpath(cli.output, fn, "pinning.mp4"); framerate=cli.fps) do io

    # iterate through paths
    for path in paths

      # data from simulation containing n trials
      data = jldopen(joinpath(path, "data_phylo.jld2"), "r")

      # get lineage heat map
      lhmap, maxtrials = lineageTracks(sample.dims, data)

      close(data)

      # empty previous content from axis
      empty!(ax)

      # hotspots
      if cli.hotspots
        heatmap!(ax, reshape(env, sample.dims); colorrange=(2, 3), lowclip=:transparent, colormap=colors)
      end

      # lineages
      heatmap!(ax, lhmap; colorrange=(1, uprbnd), lowclip=:transparent, colormap=ColorSchemes.amp)

      recordframe!(io)  # record a new frame
    end
  end
end