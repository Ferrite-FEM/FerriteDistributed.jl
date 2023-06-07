"""
    SharedEntity

Supertype for shared entities.
"""
abstract type SharedEntity end

"""
    remote_entities(::SharedEntity)

Get an iterable of pairs, containing the 1-based rank and a collection of EntityIndices on the rank.
"""
remote_entities(::SharedEntity)

# TODO the following three structs could be merged to one struct with type parameter.
# We might want to think about the design a bit.
"""
    SharedVertex <: SharedEntity

A shared vertex induced by a local vertex index and all remote vertex indices on all remote ranks.
"""
struct SharedVertex <: SharedEntity
    local_idx::VertexIndex
    remote_vertices::Dict{Int,Vector{VertexIndex}}
end

@inline remote_entities(sv::SharedVertex) = sv.remote_vertices

"""
    SharedFace <: SharedEntity

A shared face induced by a local face index and all remote face indices on all remote ranks.
"""
struct SharedFace <: SharedEntity
    local_idx::FaceIndex
    remote_faces::Dict{Int,Vector{FaceIndex}}
end

@inline remote_entities(sf::SharedFace) = sf.remote_faces

"""
    SharedEdge <: SharedEntity

A shared edge induced by a local edge index and all remote edge indices on all remote ranks.
"""
struct SharedEdge <: SharedEntity
    local_idx::EdgeIndex
    remote_edges::Dict{Int,Vector{EdgeIndex}}
end

@inline remote_entities(se::SharedEdge) = se.remote_edges
