import Measurements: ±, value, uncertainty, measurement, Measurement

function buildHSKDT(grid, hotspots, dims)
    hs = grid[map(xy->dims[1] * (xy[2] - 1) + xy[1], hotspots)]
    return KDTree(SVector.(hs))
end

function backgroundEstimate(spinput, linput, xy, dims, density, radius; nl = 1, ntests = 5)
    bckgrndEstmt = Vector{Measurement}(undef, max(1, ntests))

    for i in range(1, max(1,ntests))
        hotspots = uniformHotspots(dims, density, radius, (1,dims[2]))[2]

        # build kdtree from hotspot locations
        kdtree = buildHSKDT(xy, hotspots, dims)
        sphs = spTracks(spinput, xy, dims, kdtree, radius; nl=nl) |> last

        acc, _ = subsetHotspots(kdtree, xy, radius, sphs, linput)
        # bckgrndEstmt[i] = (acc_hotspots = acc,freq_hotspots = frq)
        bckgrndEstmt[i] = acc
    end

    return mean_and_std(bckgrndEstmt)
end

""" returns: shortest path tree root locations, shortests paths indices, visited hotspots"""
function spTracks(data, xy, dims, kdtree, radius; nl=1)

    # find shortest-paths sink locations
    roots = zeros(Int, dims[1])
    for (_, ws) in data["weights"]
        sorted = sort(ws, by = first)
        for i = 1:min(nl, length(sorted))
            @inbounds for idx in last(sorted[i])
                roots[idx] = 1
            end
        end
    end

    # get the site visits by shortest paths
    sps = zeros(Int64, reduce(*, dims))
    for (_, trcks) in data["tracks"]
        for i in range(1, min(nl, length(trcks)))
            @inbounds for idx in trcks[i]
                sps[idx] = 1
            end
        end
    end

    # # use visited sites to get visited hotspots
    r = ceil(Int, √(3) * radius)
    spHs = overlap(r, xy, findall(sps .== 1), kdtree)
    return (sps=sps, roots=roots, spHs=spHs)
end

function subsetHotspots(kdtree, xy, radius, spHS, data)
    isempty(kdtree.data) && return

    # compare hotspots visited by SP and lineages
    cnts = 0; r = ceil(Int, √(3) * radius)
    accuracy = Float64[]
    freq = zeros(Float64, size(kdtree.data))

    # calculate visited hotspots
    for trial in fContains(keys(data), "trial")

        cnts += 1
        phylo = data[trial]["phylo"]
        lHS = overlap(r, xy, collect(keys(phylo)), kdtree)

        # count how many times a HS is visited
        freq .+= lHS

        # metric -> |l ∩ s| ÷ |l ∪ s| == lis / lus
        lis = sum((spHS .== 1) .& (lHS .== 1))
        lus = sum((spHS .== 1) .| (lHS .== 1))
        push!(accuracy, lis / lus)
    end

    return (
        acc_hotspots = measurement(mean_and_std(accuracy)...), 
        freq_hotspots = freq ./ cnts
    )
end

""" get the predicted ancestors from shortest paths"""
function subsetAncestors(data, dims, roots)
    accuracy = Float64[]; cnts = 0
    survivalFreq = zeros(Float64, dims[1])

    # interpolate values between surviving ancestors
    interpol = deepcopy(roots)
    interpol[(circshift(roots, -1) .== 1) .& (circshift(roots, 1) .== 1)] .= 1

    for trial in fContains(keys(data), "trial")
        pids = unique(data[trial]["pIDs"])
        cnts += 1; matches = 0 
        
        for id in pids
            interpol[id] == 1 && (matches += 1)
            survivalFreq[id] += 1
        end
        push!(accuracy, matches / length(pids))
    end

    # figure for survival frequency
    f, _ = lines(survivalFreq ./ cnts, color=interpol,
        axis = (ylabel="1/N", xlabel="Ancestor",)
    )
    display(f)

    # fraction captured
    auc = sum(interpol .* survivalFreq ) / sum(survivalFreq)

    return (
        acc_ancestors=measurement(mean_and_std(accuracy)...),
        freq_ancestors=survivalFreq ./ cnts, 
        auc = auc
    )
end

function lineageTracks(dims, data)
    lhmap = zeros(Int64, reduce(*, dims))
    survivalnumbers = Int64[]

    for trial in fContains(keys(data), "trial")
        pids = length(unique(data[trial]["pIDs"]))
        push!(survivalnumbers, pids)

        lhmap[collect(keys(data[trial]["phylo"]))] .+= 1
    end

    return  (
        lhmap = reshape(lhmap, dims),
        surviving_ancestors = measurement(mean_and_std(survivalnumbers)...)
    )
end

function pinning(lhmap, y₀, updwn=10)
    bucket = Float64[]
    for y in range(y₀ - fld(updwn, 2), y₀ + fld(updwn, 2))
        append!(bucket, view(lhmap, :, y))
    end
    return bucket ./ maximum(lhmap)
end