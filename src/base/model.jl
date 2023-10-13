module GraphModel

using NearestNeighbors
using DataStructures
using StaticArrays
using StatsBase
using Random

include("../tools/containers.jl")
include("./commandLine.jl")
include("./dataModels.jl")
include("./paths.jl")
include("../tools/trees.jl")

import .HashMaps: HashVec, add!, remove!

export hexGraph, simulate!, linearStart!, resetGraph!, PrintSimInfo, parseArgs, DataModel, setPath

mutable struct Node
    const x::Float32
    const y::Float32
    ancestor::UInt32
    filled::Bool
    ID_1::Int32
    ID_2::Int32
    time::Float32
    nbors::UInt32
    Node(n, x, y) = new(x, y, 0, false, 0, 0, 0., n)
end

function gNodes(typ::Symbol=:hex)
    if typ == :hex 
        # https://www.redblobgames.com/grids/hexagons/#neighbors
        bottom = [(1, 0), (-1, 0), (-1, -1), (0, -1), (-1, 1), (0, 1)]
        top = [(1, 0), (-1, 0), (0, 1), (1, 1), (0, -1), (1, -1)]
        return Dict{Int,Vector{NTuple{2,Int}}}(1 => bottom, 0 => top)
    end
end

function buildMap(dimensions, nodes, T::Type, cnstr)
    lx, ly = dimensions
    graph = Vector{T}(undef, lx*ly)
    connections = zeros(UInt32, 6, lx * ly)

    for row = 1:ly, col in 1:lx
        nodeIndx = lx*(row - 1) + col
        nbors::UInt32 = 0

        for (i, (dx, dy)) in enumerate(nodes[row % 2])
            ny = row + dy
            0 < ny <= ly || continue
            nx = mod(col + dx, 1:lx)
            idx = lx*(ny - 1) + nx
            connections[i, nodeIndx] = idx
            nbors += 1
        end
        graph[nodeIndx] = cnstr(nbors, col, row)
    end
    return graph, connections
end

hexGraph(dims) = begin 
    cvr = (x,y) -> (sqrt(3)*(x - 0.5*(y%2) ), 1 + 1.5*(y-1))
    buildMap(dims, gNodes(), Node, (n,x,y) -> Node(n,cvr(x,y)...))
end

function linearStart!(graph, cntns, dimensions, env; row = 1)
    active = Vector{HashVec{UInt32}}(undef, 2)
    foreach(i->active[i] = HashVec{UInt32}(), [1,2])

    for col in 1:dimensions[1]
        nodeIndx = dimensions[1]*(row - 1) + col
        addNode!(graph, cntns, active, 0, nodeIndx, 0., env[col], col)
    end
    return active
end

@inline function nextIndex(graph, cntns, idx)
    nbors = UInt32[]
    @inbounds for nbor in view(cntns, :, idx)
        (nbor == 0 || graph[nbor].filled) && continue
        push!(nbors, nbor)
    end
    if length(nbors) == 0
        return 0
    elseif length(nbors) == 1
        return first(nbors)
    end
    return rand(nbors)
end

function addNode!(graph::Vector{Node}, cntns, active, prevIdx, newIdx, time, envID, ID_2, viable=true)

    node = graph[newIdx]

    node.filled && error("selecting an already filled node")

    node.filled = true
    node.time = time
    node.ID_1 = envID
    node.ID_2 = ID_2
    node.ancestor = prevIdx

    # don't add site to front if it already has no empty nbors
    node.nbors > 0 && viable && add!(active[envID], newIdx)

    for nbor in view(cntns, :, newIdx)
        nbor == 0 && continue

        # subtract one from all neighbors
        graph[nbor].nbors -= 1

        # if this site, or neighbor, is surrounded, then remove it
        if graph[nbor].nbors == 0 && graph[nbor].filled
            remove!(active[graph[nbor].ID_1], nbor)
        end
    end
end

function chooseNext(rates::SVector{2,Float64}, active)::Tuple{Int,Float64}
    r₁::Float64 = rates[1]*length(active[1].data)
    a₀::Float64 = r₁ + rates[2]*length(active[2].data)
    η::Float64 = a₀ * rand()
    return (η <= r₁) ? 1 : 2, a₀
end

function simulate!(graph::Vector{Node}, cntns, active, env, cli, dataModels; stopline=0)
    rates = SVector{2,Float64}([1, max(1. + cli.intensity, 0.)])
    width, height = cli.dims; time = 0.

    nrecord::Int64 = cld(width * height, cli.numberSamples)
    itrCntr::Int64 = nrecord; continueRecording::Bool = true
    
    # save snapshot of the front as it touches the top
    fSnapshot = Int32[]

    stopCondition = iszero(stopline) ? width * (height - 1) + 1 : width * (stopline - 1) + 1

    @inbounds @fastmath while true
        length(active[1].data) + length(active[2].data) == 0 && break

        if continueRecording && itrCntr == nrecord
            itrCntr = 0
            for (_, measure) in dataModels
                push!(measure.data, measure.func(active, graph, time))
            end
        end
        itrCntr += 1
        
        # find which cell will grow next
        (label, R) = chooseNext(rates, active)

        # update time
        time += -log(rand()) / R

        # get node indices
        nodeIndx = rand(active[label].data)
        newIdx = nextIndex(graph, cntns, nodeIndx)

        # fill new node, add to front if viable
        addNode!(graph, cntns, active, nodeIndx, newIdx, time, env[newIdx], graph[nodeIndx].ID_2, !(rates[env[newIdx]] == 0.))

        # stop recording data for later analysis
        if newIdx >= stopCondition && isempty(fSnapshot) 
            fSnapshot = vcat(active[1].data, active[2].data)
            continueRecording = false
            iszero(stopline) || break
        end
    end
    return (fSnapshot=fSnapshot, time=time)
end

function resetGraph!(graph, cntns)
    @inbounds @simd for i in eachindex(graph)
        graph[i].filled = false
        graph[i].ancestor = 0
        graph[i].ID_1 = 0
        graph[i].ID_2 = 0
        graph[i].time = 0.
        graph[i].nbors = sum(view(cntns, :, i) .> 0)
    end
end

end