using FileIO, JLD2

function createEnvironment!(opts)
  upnt(c, k, v) = merge(c, (Symbol(k) => v,))

  opts = upnt(opts, "dims", (opts.width, opts.height))

  if !isnothing(opts.landscape)
    inputopts, hotspots = load(opts.landscape, "opts", "htspts")
    opts = upnt(opts, "density", inputopts.density)
    opts = upnt(opts, "radius", inputopts.radius)
    opts = upnt(opts, "dims", inputopts.dims)
    return env, hotspots, opts
  end

  if opts.env_type == "circle"
    isnothing(opts.separation) && error("must provide a separation between circles")
    obs = obstacles_vertical(opts.radius, opts.separation, opts.width, opts.height)
  else
    _density = min(opts.density, 1.0)
    obs = obstacles_uniform(opts.radius, _density, opts.width, opts.height, opts.gap)
  end

  env, num = applyObstacles!(obs, opts.radius, opts.width, opts.height, min(opts.density, 1.0))
  density = sum(env .== 2) ./ length(env)
  opts = upnt(opts, "density", density)
  return env, obs[1:num], opts
end

function obstacles_hexgrid(Radius, LatticeSurfSep, width, height, gap=0)
  ObsCenterList = []

  # The separation between obstacle centers is more useful
  LatticeCentreSep = LatticeSurfSep + 2.0 * Radius

  #Shortest dist is between horizontal neighbours
  xObstaclesNum = floor(Int, float(width) / LatticeCentreSep)
  yObstaclesNum = floor(Int, float(height) / LatticeCentreSep)

  for obsrow in range(0, yObstaclesNum)
    ypos = ceil(Int, (2 / sqrt(3)) * Radius + obsrow * LatticeCentreSep)
    for obscol in range(0, xObstaclesNum)
      xpos = ceil(Int64, Radius + (obscol + 0.5 * (obsrow % 2)) * LatticeCentreSep)
      if 1 <= xpos <= width && 1 <= ypos + gap <= height
        push!(ObsCenterList, (xpos, ypos + gap))
      end
    end
  end
  return ObsCenterList, width
end

function obstacles_uniform(radius, phi, DomainWidth, DomainHeight, gap=0, safety=true)
  # r = 2 * radius / sqrt(3)
  r = radius

  #Units: obstacles per site. hence the 2/root(3)
  # numberdensity = log(1 - phi) / (-π * r^2 * 2 / sqrt(3))
  numberdensity = log(1 - phi) / (-π * r^2)

  ObsNum = numberdensity * (DomainHeight - r) * DomainWidth
  ObsNum = safety ? 2 * ObsNum : ObsNum

  ObsCenterList = []

  #Create the list of obstacle center positions, using uniform random values
  for _ in range(1, round(Int, ObsNum))
    x = rand(1:DomainWidth)
    y = rand((gap + round(Int, r)):DomainHeight)
    push!(ObsCenterList, (x, y))
  end
  return ObsCenterList
end

function obstacles_vertical(radius::Int, sep::Int, lx::Int, ly::Int)
  obstacles = NTuple{2, Int64}[]
  center = floor(Int, lx / 2)
  r = ceil(2 * radius / sqrt(3))
  foreach(xy -> push!(obstacles, xy), [(center, y) for y in (2 * r):sep:(ly - r)])
  return obstacles
end

function applyObstacles!(obsCenters, radius, lx, ly, area_fraction=Inf; envObj=false)
  envObjs = Vector{Vector{Int}}(undef, lx * ly)
  env = ones(UInt8, (lx * ly))
  r = sqrt(3) * radius
  R = ceil(Int, r)

  offset(x, y) = (sqrt(3) * (x - 0.5 * (y % 2)), 1 + 1.5 * (y - 1))

  _sites_covered = 0
  _obs_number = 0
  for (j, (ox, oy)) in enumerate(obsCenters)
    for dx in (-R):R, dy in (-R):R
      sum(abs2, offset(ox + dx, oy + dy) .- offset(ox, oy)) <= r^2 || continue
      ny = dy + oy
      1 <= ny <= ly || continue
      nx = mod(dx + ox, 1:lx)

      node_index = lx * (ny - 1) + nx
      if env[node_index] == 1
        _sites_covered += 1
        #issue i dont' think that the indices are properly tracked
        if isassigned(envObjs, node_index)
          push!(envObjs[node_index], j)
        else
          envObjs[node_index] = Int[j]
        end
      end
      env[node_index] = 2
    end

    _obs_number += 1
    _area_fraction = _sites_covered / length(env)
    if _area_fraction >= area_fraction ||
       isapprox(_area_fraction, area_fraction, rtol=0.005)
      break
    end
  end

  # update obs to satisfy the area_fraction condition
  envObj && return env, _obs_number, envObjs
  return env, _obs_number
end

function uniform(lx, ly, density, radius, regionBounds)
  env = ones(UInt8, (lx * ly))
  hotspots = Tuple{Int64, Int64}[]

  low, hi = regionBounds

  radii = ceil(Int, sqrt(3) * radius)
  xpool = collect((radii + 1):(lx - radii))
  ypool = collect((radii + low):(hi - radii))

  cvr = (x, y) -> (sqrt(3) * (x - 0.5 * (y % 2)), 1 + 1.5 * (y - 1))

  nfill = 1
  cntr = 0
  while nfill <= floor(Int, lx * ly * density)
    cx = rand(xpool)
    cy = rand(ypool)
    in((cx, cy), hotspots) && continue
    nfill += fillCircle!(lx, ly, radii, env, cx, cy, cvr)
    push!(hotspots, (cx, cy))
    cntr += 1
    cntr == 100_000 && break
  end
  return env, hotspots, sum(env .== 2) / (lx * ly)
end

function circles(lx, ly, radii, sep::Int)
  cvr = (x, y) -> (sqrt(3) * (x - 0.5 * (y % 2)), 1 + 1.5 * (y - 1))
  radius = ceil(Int, sqrt(3) * radii)
  env = ones(UInt8, lx * ly)

  hotspots = NTuple{2, Int64}[]

  for (cx, cy) in [(floor(Int, lx / 2), y) for y in (2 * radius):sep:(ly - radius)]
    fillCircle!(lx, ly, radius, env, cx, cy, cvr)
    push!(hotspots, (cx, cy))
  end
  return env, hotspots, sum(env .== 2) / (lx * ly)
end

function fillCircle!(lx, ly, radii, env, cx, cy, cvr)
  x0, y0 = cvr(cx, cy)
  mx = radii^2
  nfilled = 0
  for dx in (-radii - 9):(radii + 9), dy in (-radii - 9):(radii + 9)
    x, y = cvr(cx + dx, cy + dy)

    if (x - x0)^2 + (y - y0)^2 <= mx
      ny = mod(dy + cy, 1:ly)
      nx = mod(dx + cx, 1:lx)
      nidx = lx * (ny - 1) + nx
      env[nidx] == 1 && (nfilled += 1)
      env[nidx] = 2
    end
  end
  return nfilled
end
