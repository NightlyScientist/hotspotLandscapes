using JLD2
import DataFrames: DataFrame
import StatsBase: mean, std, mean_and_std
include("binUtils.jl")

function filterRead(path, needle, join=true)
    haystack = readdir(path, join=join)
    isempty(haystack) && return haystack
    filter(t->contains(t,needle), haystack)
end

function setPath(pth::String, overwrite::Bool, cli, recursive=false)::String
    path = islink(pth) ? readlink(pth) : pth
    if ispath(path) && ~isempty(readdir(path))
        overwrite && rm(path, recursive=recursive)
        overwrite || error(" ! $path not empty ! | $cli")
    end
    mkpath(path)
    return path
end

function parse_input(pth)
    d = Dict{Symbol,Any}()
    for line in readlines(pth)
        sp = split(line, ":")
        length(sp) > 2 && error("unable to parse to dictionary")
        d[first(sp) |> Symbol] = strip(String(last(sp)), " ")
    end
    return d
end

function fContains(haystack, needle)
    filter(t->contains(t,needle), haystack)
end

function collectionTable(path, radius::Float64, gb = :density)
    data = []
    for root in readdir(path, join=true)
        prtns = getCollections(root, radius, gb)
        foreach(prtn->append!(data, prtn), prtns)
    end
    DataFrame(data)
end

"""
    group elements in y using data from first axis of x within dist radius
"""
function group(xi, yi, radius::Float64)
    str = sortperm(xi)
    x = xi[str]
    y = yi[str]
    # bug: some elements are different, e.g. don't have :envType field
    N = [Any[t] for t in y]
    w = fill(true, size(x,1))
    for i in range(2, size(x,1))
        (abs(x[i-1] - x[i]) < radius && w[i-1]) || continue
        append!(N[i], N[i-1])
        w[i-1] = false
    end
    return N[w]
end

function getCollections(path, radius::Float64 = 0.005, gb = :density)
    in(gb, [:density, :intensity, :radius]) || return

    data = []
    parameters = []
    for grp in filterRead(path, "collection")
        for rnf in readdir(grp, join=true)
            rprms = merge(load(rnf * "/data.jld2", "cli"), (:path => rnf,))
            push!(data, rprms)
            push!(parameters, haskey(rprms, :parameter) ? rprms.parameter : gb)
        end
    end

    @assert allequal(parameters)
    prms = getfield.(data, parameters |> first |> Symbol)

    partitions = group(prms, data, radius)
    prtns = Vector{Vector{Any}}(undef, length(partitions))

    for i in eachindex(partitions)
        prtn = partitions[i]
        prm = first(prtn).parameter |> Symbol
        mn, sd = getfield.(prtn, prm) |> mean_and_std
        prtns[i] = Any[]
        for j in eachindex(prtn)
            push!(prtns[i], merge(prtn[j], (:meanParam => mn, :stdParam => sd,)))
        end
    end
    return prtns
end