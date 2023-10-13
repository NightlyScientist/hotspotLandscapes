import Base.Threads: @spawn, fetch, nthreads
import JLD2: jldopen, jldsave
import FileIO: load
import StatsBase: mean, kurtosis, std, median, mode
import Colors: RGBA
import ColorSchemes: get, amp, roma, hot, ColorScheme
import Measurements: value, uncertainty, ±, measurement
using EasyFit
using CairoMakie
using DataStructures
using LaTeXStrings

include("../../src/tools/parseFiles.jl");
include("../../src/tools/binUtils.jl");
include("../../src/base/model.jl");
include("../../src/tools/containers.jl")

import .GraphModel: gNodes, buildMap

# using P.B.C get the arc distance between two points
function periodicBoundaries(halfwidth::Int, width::Int)
    f = let h=halfwidth, w=width
        x -> (x > h) ? w - x : x
    end
    return f
end;

function passFltr(y , x = collect(1:length(y)); uplim = 1000, offset = 2)
    
    ps = (x .!== NaN) .& (y .!== NaN) .& (x .> 0) .& (y .> 0) .& (x .< uplim)
    return x[ps] .+ offset, y[ps]
end

"""
    measure the msd from genetic lineages
"""
function lineageMSD!(msd, counts, gx, gy, srcs, phylo, pbc; offset = 2)
    for src in srcs
        yref = gy[src]
        xref = gx[src]
        current = src; check = true
        @inbounds while current != 0
            if check && abs(gy[current] - gy[src]) < offset
                xref = gx[current]
                yref = gy[current]
                current = phylo[current]
                continue
            else
                check = false
            end
            x = pbc(abs(xref - gx[current])) ^2
            y = abs(yref-gy[current]) + 1
            msd[y] += x
            counts[y] += 1
            current = phylo[current]
        end
    end
end;

"""
    measure the msd from the shortest paths
"""
function spMSD!(msd, counts, gx, gy, tracks,  pbc; nl=1, offset = 2)
    for (src, trcks) in tracks
        yref = gy[src]
        xref = gx[src]
        for i in eachindex(trcks)
            i <= nl || continue
            trck = trcks[i]
            check = true
            for j in eachindex(trck)
                if check && abs(gy[trck[j]] - gy[src]) < offset
                    xref = gx[trck[j]]
                    yref = gy[trck[j]]
                    continue
                else
                    check = false
                end
                x = pbc(abs(xref - gx[trck[j]])) ^2
                y = abs(yref-gy[trck[j]]) + 1
                msd[y] += x
                counts[y] += 1
            end
        end
    end 
end

function processSamples(partition; smpl=true, nl = 1, offset = 2)

    params = rand(partition)
    width, height = params.dims
    pvals = []

    # build map geometry
    gxy, _ = buildMap(params.dims, gNodes(), NamedTuple, (n,x,y)->(x=x, y=y))
    gx = getfield.(gxy, :x); gy = getfield.(gxy, :y)

    # periodic boundaries
    pbc = periodicBoundaries(div(width, 2), width)

    # lineage msd
    msd = zeros(Float64, height + 1)
    counts = zeros(Float64, height + 1)

    # SP msd
    spmsd = zeros(Float64, height + 1)
    spcounts = zeros(Float64, height + 1)

    for ithp in partition

        push!(pvals, getfield(ithp, Symbol(ithp.parameter)))

        # measure tracks msd
        if params.density > 0.0 && islink(ithp.path * "/spaths.jld2")
            tracks = load(readlink(ithp.path * "/spaths.jld2"), "tracks")
            spMSD!(spmsd, spcounts, gx, gy, tracks, pbc; 
                nl = nl, offset = offset
            )
        end

        # measure lineage msd
        jldopen(joinpath(ithp.path, "data.jld2"), "r") do file
            for trial in fContains(keys(file), "trial")
                srcs = file[trial]["sources"]
                phylo = file[trial]["phylo"]

                lineageMSD!(msd, counts, gx, gy, srcs, phylo, pbc;
                    offset = offset 
                )
            end
        end
    end
    return (msd ./ counts), (spmsd ./ spcounts), pvals
end;    addargs!(sts, name; kwargs...) = add_arg_table!(sts, "--$(name)", Dict(kwargs))

inputs = begin 
    sts = ArgParseSettings()
    addargs!(sts, "path", arg_type=String, required=true)
    addargs!(sts, "offset", arg_type=Int, default=2) 
    addargs!(sts, "smpl", action=:store_true)
    addargs!(sts, "nl", default=5)

    namedtuple(parse_args(sts))
end

fSets = append!([getCollections(pth) for pth in pths]...)

function psub(basedir, ptn, smpl, offset)
    fname = splitdir(first(ptn).path) |> last

    isfile(basedir * "/$(fname).jld2") && return nothing

    msd, spmsd, pvals = processSamples(ptn; 
        smpl=smpl, offset = offset
    )

    jldsave(basedir * "/$(fname).jld2"; 
        msd = msd,
        spmsd = spmsd,
        config = ptn |> first,
        configs = ptn,
        pval = measurement(mean_and_std(pvals)...),
    )
    return mean(pvals)
end

pths = [inputs.path]

pvals = fetch.([@spawn psub($basedir, $ptn, $false, $(inputs.offset)) for ptn in fSets])

@info pvals