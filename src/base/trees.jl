module GenealogicalTree
include("../common/indexTools.jl")
export TreeNode

mutable struct TreeNode
  parent::UInt32
  children::Set{UInt32}
  const time::Float32
  const position::UInt32
end

graphSources(opts::NamedTuple) = graphSources(opts.ref_line, opts.width, opts.height)
function graphSources(ref_line, width, height)
  ref_line = iszero(ref_line) ? height : ref_line
  return convert.(UInt32, collect(1:width) .+ width * (ref_line - 1))
end

# doc: construct ancestry trees as a dictionary (julia implemented as hash table)"
function genealogy(graph, front)
  phylogeny = Dict{UInt32,UInt32}()

  # .associate originID to list of branch splits
  branchPoints = Dict{UInt32,Vector{UInt32}}()
  branchTimes = Dict{UInt32,Float32}()

  @inbounds for nodeIdx in front
    ID::UInt32 = graph[nodeIdx].ID_2
    current::UInt32 = nodeIdx

    # refactor: added source times to branch points
    branchTimes[current] = graph[current].time

    while true
      if haskey(phylogeny, current)
        haskey(branchPoints, ID) || (branchPoints[ID] = UInt32[])
        push!(branchPoints[ID], current)
        branchTimes[current] = graph[current].time
        break
      end

      phylogeny[current] = graph[current].ancestor
      current = graph[current].ancestor
      current == 0 && break
    end
  end
  # @debug foreach(k->println(Pair(k, phylogeny[k])), front)
  foreach(unique!, values(branchPoints))
  return phylogeny, branchPoints, branchTimes
end

# doc: construct ancestry trees as a dictionary of common ancestors
function buildTree(genealogy, geneticLabels, branchPoints, branchTimes, sources)
  tree = Dict{UInt32,TreeNode}()

  for i in eachindex(sources)
    source::UInt32 = sources[i]
    id::UInt32 = geneticLabels[i]
    # .lineage can have just one line, skip it
    haskey(branchPoints, id) || continue
    divs = branchPoints[id]

    # .do not reset source nodes if previous track passed through here
    haskey(tree, source) || (tree[source] = TreeNode(0, Set([]), branchTimes[source], source))

    current::UInt32 = source
    child::UInt32 = source

    while current != 0
      if in(current, divs)
        # .do not push current into children list
        if current != source
          if haskey(tree, current)
            push!(tree[current].children, child)
          else
            # refactor: previously initialized parent with -1, now typemax(UInt32)
            tree[current] = TreeNode(typemax(UInt32), Set(child), branchTimes[current], current)
          end
          tree[child].parent = current
        end
        child = current
      end
      @inbounds current = genealogy[current]
      # @inbounds current = graph[current].ancestor
    end
  end
  return tree
end

function traverse(start::UInt32, tree::Dict{UInt32,TreeNode})
  parent = start
  track = UInt32[]
  nodes = UInt32[]
  while parent != 0
    push!(track, parent)
    push!(nodes, tree[parent])
    parent = tree[parent].parent
  end
  # @debug join(reverse(track), " --> ")
  return track, nodes
end

function treeMatrix(tree::Dict{UInt32,TreeNode})
  maxLength = 1
  treeMatrix = zeros(UInt32, 10, 100)
  for i in eachindex(sources)
    tmp = traverse(sources[i], tree)
    cntr = length(tmp)
    if cntr > maxLength
      maxLength = cntr
    end
    treeMatrix[i, 1:cntr] = reverse(tmp)
  end
end

periodicBoundaries(w::Int, x::Integer) = (x > cld(w, 2)) ? w - x : x
periodicBoundaries(w::Int, x::Float64) = (x > /(w, 2)) ? w - x : x

# doc: find the least common ancestor of two lineages
function lca_times(tree::Dict{UInt32,TreeNode}, a::Int32, b::Int32, ai::UInt32, bi::UInt32)::UInt32
  ai != bi && return 0

  STOP = typemax(UInt32)
  parentOne::UInt32 = a
  parentTwo::UInt32 = b
  next = tree[a].time > tree[b].time ? 1 : 2

  # .check that one of the sources is not already a node of the other
  parentOne == parentTwo && return 0

  while parentOne != STOP && parentTwo != STOP
    if next == 1
      parentOne = tree[parentOne].parent
    else
      parentTwo = tree[parentTwo].parent
    end
    parentOne == parentTwo && break
    next = tree[parentOne].time > tree[parentTwo].time ? 1 : 2
  end
  @debug @assert isequal(parentOne, parentTwo) " falsely found a common ancestor"
  # .skip comparisons when common ancestor is one of the sample points
  (parentTwo == a || parentOne == b) && return 0
  parentOne == parentTwo && return parentOne
  return 0
end
end