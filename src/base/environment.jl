# using Accessors
using FileIO
using JLD2

function createEnvironment!(cli)
    upnt(c, k, v) = merge(c, (Symbol(k) => v,))

    tgtDens = min(cli.density, 1.0)

    # set region to place hotspots
    if isnothing(cli.hsRegion) 
        # cli = @set cli.hsRegion = (1,cli.dims[2])
        cli = upnt(cli, "hsRegion", (1, cli.dims[2]))
    end

    if !isnothing(cli.landscape)
        input = load(cli.env)
        inputCli = input["cli"]
        cli = upnt(cli, "density", inputCli.density)
        cli = upnt(cli, "radius", inputCli.radius)
        cli = upnt(cli, "dims", inputCli.dims)

        # cli = @set cli.density = inputCli.density
        # cli = @set cli.radius = inputCli.radius
        # cli = @set cli.dims = inputCli.dims
        hotspots = input["hotspots"]
    else
        if cli.envType == "circle" uniformHotspots
            # env, hotspots, density = singleCircle(cli.dims, cli.radius, fld(cli.dims[1], 2), 2*cli.radius)
            env, hotspots, density = singleCircle(cli.dims, cli.radius, fld.(cli.dims, 2)...)
        else
            env, hotspots, density = uniformHotspots(cli.dims, tgtDens, cli.radius, cli.hsRegion)
        end

        # cli = @set cli.density = density
        cli = upnt(cli, "density", density)
    end
    return env, hotspots, cli
end

function uniformHotspots(dims, density, radius, regionBounds)
    env = ones(UInt8, *(dims...))
    hotspots = Tuple{Int64,Int64}[]

    lx, ly = dims
    low, hi = regionBounds

    radii = ceil(Int, sqrt(3) * radius) 
    xpool = collect(radii+1:lx-radii)
    ypool = collect(radii+low:hi-radii)

    cvr = (x,y) -> (sqrt(3)*(x - 0.5*(y%2) ), 1 + 1.5*(y-1))

    nfill = 1; cntr = 0
    # task: change to account for density of just a small region
    while nfill <= floor(Int, *(dims...) * density)
        cx = rand(xpool)
        cy = rand(ypool)
        in((cx,cy), hotspots) && continue
        nfill += fillCircle!(lx, ly, radii, env, cx, cy, cvr)
        push!(hotspots, (cx, cy))
        cntr += 1
        cntr == 100_000 && break
    end
    return env, hotspots, sum(env .== 2) / *(dims...)
end

function singleCircle(dims, radii, cx, cy)
    # bug: r -> sqrt(3) * radius
    radius = ceil(Int, sqrt(3) * radii) 
    cvr = (x,y) -> (sqrt(3)*(x - 0.5*(y%2) ), 1 + 1.5*(y-1))
    env = ones(UInt8, *(dims...))
    hotspots = NTuple{2,Int64}[]
    fillCircle!(dims..., radius, env, cx, cy, cvr)
    push!(hotspots, (cx, cy))
    return env, hotspots, sum(env .== 2) / *(dims...)
end

function fillCircle!(lx, ly, radii, env, cx, cy, cvr)
    x0, y0 = cvr(cx, cy)
    mx = radii ^2; nfilled = 0
    for dx in (-radii-9):(radii+9), dy in (-radii-9):(radii+9)
        x, y = cvr(cx + dx, cy + dy)

        if (x - x0) ^ 2 + (y - y0) ^2 <= mx
            ny = mod(dy + cy, 1:ly)
            nx = mod(dx + cx, 1:lx)
            nidx = lx*(ny - 1) + nx
            env[nidx] == 1 && (nfilled += 1)
            env[nidx] = 2
        end
    end
    return nfilled
end