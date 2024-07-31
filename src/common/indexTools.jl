# doc: convert linear indices to cartesian (x,y)
reconstructCoordinates(index::Integer, w, h)::NTuple{2,Integer} = Tuple(CartesianIndices((w, h))[index])
function reconstructCoordinates(index::Vector{<:Integer}, w, h)::NTuple{2,Vector{Integer}}
  i2s = CartesianIndices((w, h))
  xy = Tuple.(i2s[index])
  return first.(xy), last.(xy)
end

# doc: convert linear indices to cartesian (x,y) in real space
function reconstructCoordinatesReal(index::Integer, w, h)::NTuple{2,Float64}
  x, y = reconstructCoordinates(index, w, h)
  return x - 0.5 * isodd(y), y
end

function reconstructCoordinatesReal(index::Vector{<:Integer}, w, h)::NTuple{2,Vector{Float64}}
  x, y = reconstructCoordinates(index, w, h)
  return x .- 0.5 * isodd.(y), y
end