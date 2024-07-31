include("../base/model.jl")
include("../tools/containers.jl")
include("../tools/parseFiles.jl")

using JLD2, CairoMakie, ColorSchemes, Colors
using ArgParse
import .Model: gNodes, buildMap
import ColorSchemes: tab10, ColorScheme
import Colors: RGBA


function colorscheme_alpha(cscheme::ColorScheme, alpha::T = 0.5; 
        ncolors=12) where T<:Real
    return ColorScheme([RGBA(get(cscheme, k), alpha) for k in range(0, 1, length=ncolors)])
end

addArgs!(sts, name; kwargs...) = add_arg_table!(sts, "--$(name)", Dict(kwargs))

function ParseArgs()
    sts = ArgParseSettings()
    addArgs!(sts, "input", required=true)
    addArgs!(sts, "output", required=true)
    addArgs!(sts, "overwrite", action=:store_true)
    addArgs!(sts, "fps", default=14, arg_type=Int)
    addArgs!(sts, "lineages", action=:store_true)
    addArgs!(sts, "hotspots", action=:store_true)
    addArgs!(sts, "colors", action=:store_true)
    return parse_args(sts) |> namedtuple
end

cli = ParseArgs()
@info "cli:" cli

if ~cli.lineages && ~cli.colors 
    @warn "need at least one of: colors, lineages"
    exit()
end

@info "saving recording to:" cli.output 

fn = let
    s = splitdir(cli.input)
    last(s) == "" ? first(s) |> splitdir |> last : last(s)
end

@info "using data from: " fn
setPath(cli.output * "/" * fn, cli.overwrite, cli, true)

let cmdline = cli
    file = jldopen(joinpath(cmdline.input, "data.jld2"), "r")
    env = file["env"]
    cli = file["cli"]

    cli.animation == false && return println("no animation found")

    # theme and color
    set_theme!(figure_padding = 00)
    colors = colorscheme_alpha(ColorSchemes.grays, 0.3, ncolors=3)
    greens = colorscheme_alpha(ColorSchemes.Greens, 0.5, ncolors=2)
    
    # create figure
    fig = Figure(resolution=cli.dims .* (√(3), 1.5), backgroundcolor=:transparent)
    ax = Axis(fig[1,1], backgroundcolor=:white)

    # set image limits
    hidedecorations!(ax)
    limits!(ax, 1, cli.dims[1], 1, cli.dims[2])

    # used for updating lineages
    tmp = zeros(Int64, size(env))

    shiftValues(x) = x == 0 ? 0 : 5

    record(fig, joinpath(cmdline.output, fn, "growth.mp4"); framerate=cmdline.fps) do io
        for trial in fContains(keys(file), "trial")
            # create object for snapshots: list->tuple(growth, lineages)
            snapshots = file[trial]["animation"]

            for i in eachindex(snapshots)
                # empty previous content from axis
                empty!(ax)
                _, ids, lineages = snapshots[i]

                # id colors
                if cmdline.colors
                    heatmap!(ax, reshape(ids, cli.dims), colormap=:hsv, colorrange=(1,cli.dims[1]), lowclip=:transparent)
                else
                    heatmap!(ax, shiftValues.(reshape(ids, cli.dims)), colorrange=(1,2), lowclip=:transparent, highclip=(:green, 0.2))
                end

                # lineages
                if cmdline.lineages && length(lineages) > 0
                    tmp .= 0; tmp[lineages] .= 3
                    color = cmdline.colors ? :black : :black
                    heatmap!(ax, reshape(tmp, cli.dims), colorrange=(1,2), lowclip=:transparent, highclip=color)
                end
          
                # hotspots
                if cmdline.hotspots
                    heatmap!(ax, reshape(env, cli.dims), colorrange=(2,3), lowclip=:transparent, colormap=colors)
                end

                # scale bar
                lines!(ax, [10,110], [490, 490], linewidth=14, color=:black)
                text!(ax, L"100", position = (60, 470), align = (:center, :center), fontsize = 40)

                recordframe!(io)  # record a new frame
            end
        end
    end

    close(file)
end
