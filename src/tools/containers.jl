module HashMaps
    export HashVec, add!, remove!

    abstract type HashContainer end

    mutable struct HashVec{Typ} <: HashContainer
        data::Vector{Typ}
        hashd::Dict{Typ,Int64}

        function HashVec{Typ}() where {Typ<:Any}
            self = new{Typ}()
            self.data = Vector{Typ}(undef, 0)
            self.hashd = Dict{Typ,Int64}()
            return self
        end
    end

    # function add!(c::HashArr{T}, x::T) where {T} 
    function add!(c, x) 
        haskey(c.hashd, x) && return false
        push!(c.data, x)
        c.hashd[x] = length(c.data)
        return
    end

    # function remove!(c::HashArr{T}, x::T) where {T}
    function remove!(c, x) 
        i = get(c.hashd, x, Nothing)
        i === Nothing && return
        laste = last(c.data)
        lastI = lastindex(c.data)
        c.data[i], c.data[lastI] = c.data[lastI], c.data[i]
        c.hashd[laste] = i
        delete!(c.hashd, x)
        pop!(c.data)
        return
    end
end

# swap last element and pop
function swapRemove!(c::Vector{<:Real}, idx::Int)
    lastI = lastindex(c)
    c[idx], c[lastI] = c[lastI], c[idx]
    pop!(c)
    return
end

dictkeys(d::Dict) = (collect(keys(d))...,)
dictvalues(d::Dict) = (collect(values(d))...,)

function namedtuple(d::Dict{Symbol,T}) where {T}
    NamedTuple{dictkeys(d)}(dictvalues(d))
end

function namedtuple(d::Dict{String,T}) where {T}
    NamedTuple{Symbol.(dictkeys(d))}(dictvalues(d))
end