mutable struct TreeNode
    parent::Int64
    children::Set{Int64}
    const time::Float32
    const position::Int64
end

function splitLocations(graph, branching)
    # select one of the trees 
    for (originID, divs) in branching
        println("id is ", originID)
        bn = graph[filter(x->x!=0, divs)]
        return getfield.(bn, :x), getfield.(bn, :y)
    end
end

function buildTree(graph, branching, sources)
    tree = Dict{Int,TreeNode}()
    for i in eachindex(sources)
        source = sources[i]
        id = graph[source].ID_2
        divs = branching[id]

        # do not reset source nodes if previous track passed through here
        haskey(tree, source) || (tree[source] = TreeNode(0, Set([]), graph[source].time, source))

        current = source
        child = source

        while current != 0
            if in(current, divs)
                # do not push current into children list
                if current != source
                    if haskey(tree, current)
                        push!(tree[current].children, child)
                    else
                        tree[current] = TreeNode(-1, Set(child), graph[current].time, current)
                    end
                    tree[child].parent = current
                end
                child = current
            end
            @inbounds current = graph[current].ancestor
        end
    end
    # @debug foreach(println, collect(sort(tree)));
    return tree
end

function traverse(start, tree)
    parent = start
    track = []
    nodes = []
    while parent != 0
        push!(track, parent)
        push!(nodes, tree[parent])
        parent = tree[parent].parent
    end
    # @debug join(reverse(track), " --> ")
    return track, nodes
end

function treeMatrix(tree)
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

function lca_times(tree, a, b, ai, bi)
    # continue if two are the diff ids (diff trees)
    ai != bi && return 0

    # parents = Int64[a, b]
    parentOne::Int = a
    parentTwo::Int = b
    whch = tree[a].time > tree[b].time ? 1 : 2

    # check that one of the sources is not already a node of the other
    parentOne == parentTwo && return parentOne

    # while parents[1] != 0 || parents[2] != 0
    while parentOne != 0 && parentTwo != 0
        if whch == 1
            parentOne = tree[parentOne].parent
        else
            parentTwo = tree[parentTwo].parent
        end
        parentOne == parentTwo && return break
        whch = tree[parentOne].time > tree[parentTwo].time ? 1 : 2
    end
    @debug @assert isequal(parentOne, parentTwo) " falsely found a common ancestor"
    return parentOne
end