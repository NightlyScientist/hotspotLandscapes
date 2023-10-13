module WeightedGraphs
    using FileIO, JLD2
    using StaticArrays, DataStructures
    export sPaths

    mutable struct Node
        weight::Float32
        ancestor::Vector{UInt32}
        closed::Bool
        Node(w = Inf, ancst = [0]) = new(w, ancst, false)
    end

    function dijkstra!(graph, connections, source, env)
        # reset fields (is this faster?)
        @inbounds @simd for i in eachindex(graph)
            graph[i].closed = false
            graph[i].weight = Inf
            graph[i].ancestor = [0]
        end

        graph[source].weight = 0

        active = PriorityQueue{UInt32,Float32}()
        active[source] = 0.

        while ~isempty(active)
            # pick min entry
            nodeIdx = dequeue!(active)

            # investigate: nodeIdx already in SP list,then stop
            # haskey(sp, nodeIdx) && break

            graph[nodeIdx].closed = true
            weight = graph[nodeIdx].weight

            # search around node, and update neighbors
            @inbounds for n in range(1,6)
                nn = connections[n, nodeIdx]
                (nn == 0 || graph[nn].closed) && continue
                oldCost = graph[nn].weight
                newCost = env[nodeIdx] + weight
                if newCost == oldCost
                    push!(graph[nn].ancestor, nodeIdx)
                elseif newCost < oldCost
                    active[nn] = newCost
                    graph[nn].weight = newCost
                    graph[nn].ancestor = [nodeIdx]
                end
            end
        end
    end

    function tracePaths(graph, targets::Vector{Int})
        tracks = Vector{Vector{UInt32}}(undef, length(targets))
        for i in eachindex(targets)
            track = UInt32[]
            current = targets[i]
            @inbounds while current != 0
                push!(track, current)
                # randomly select ancestor (all of equal weight)
                current = rand(graph[current].ancestor)
            end
            tracks[i] = reverse(track)
        end
        return tracks
    end

    function sPaths(graph, cntns, srcs::Vector{<:Int}, tgts, env, nth)
        tracks = Dict{Int,Vector{Vector{Int}}}()
        weights = Dict{Int,Vector{Tuple{Float64,Vector{Int}}}}()
        # sps = fill(false, reduce(*, dims))

        for src in srcs
            wp, tks = sPaths(graph, cntns, src, tgts, env, nth)
            tracks[src] = tks
            weights[src] = wp
        end
        return weights, tracks
    end

    function sPaths(graph, cntns, src::Int, tgts, env, nth)
        dijkstra!(graph, cntns, src, env)

        weights = getfield.(graph[tgts], :weight)
        srtr = sortperm(weights)

        # group tgts by similar (sorted) weight 
        partitions = groupSame(weights[srtr], tgts[srtr])

        # randomly select representative element in class
        subset = rand.(partitions)

        # trace paths of nth smallest weights
        tracks = tracePaths(graph, subset[1:min(length(subset), nth)])

        @debug (@assert allequal(map(p->allequal(weights[p]), partitions)))

        # -> (weight, tgts with this weight)
        wp = []; foreach(p->push!(wp, (weights[first(p)], p)), partitions)
        return wp, tracks
    end

    function groupSame(x, y)
        N = [[t] for t in y]
        w = fill(true, length(x))
        for i in range(2,length(x))
            if x[i-1] == x[i] && w[i-1]
                append!(N[i], N[i-1])
                w[i-1] = false
            end
        end
        return N[w]
    end
end