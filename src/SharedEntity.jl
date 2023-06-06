"""
"""
abstract type SharedEntity end

"""
"""
remote_entities(::SharedEntity)

# TODO the following three structs could be merged to one struct with type parameter.
# We might want to think about the design a bit.
"""
"""
struct SharedVertex <: SharedEntity
    local_idx::VertexIndex
    remote_vertices::Dict{Int,Vector{VertexIndex}}
end

@inline remote_entities(sv::SharedVertex) = sv.remote_vertices

"""
"""
struct SharedFace <: SharedEntity
    local_idx::FaceIndex
    remote_faces::Dict{Int,Vector{FaceIndex}}
end

@inline remote_entities(sf::SharedFace) = sf.remote_faces

"""
"""
struct SharedEdge <: SharedEntity
    local_idx::EdgeIndex
    remote_edges::Dict{Int,Vector{EdgeIndex}}
end

@inline remote_entities(se::SharedEdge) = se.remote_edges
