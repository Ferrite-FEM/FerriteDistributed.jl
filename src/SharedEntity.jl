"""
    SharedEntity

Supertype for shared entities.
"""
abstract type SharedEntity <: Entity end

"""
    remote_entities(::SharedEntity)

Get an iterable of pairs, containing the 1-based rank and a collection of EntityIndices on the rank.
"""
remote_entities(::SharedEntity)


"""
    SharedVertex <: SharedEntity

A shared vertex induced by a local vertex index and all remote vertex indices on all remote ranks.
"""
struct SharedVertex <: SharedEntity
    unique_local_representation::VertexRepresentation # Identify via node in grid
    local_vertices::Vector{VertexIndex}
    remote_vertices::Dict{Int,Vector{VertexIndex}}
end

@inline local_entities(sv::SharedVertex) = sv.local_vertices
@inline remote_entities(sv::SharedVertex) = sv.remote_vertices


"""
    SharedFace <: SharedEntity

A shared face induced by a local face index and all remote face indices on all remote ranks.
"""
struct SharedFace <: SharedEntity
    unique_local_representation::FaceRepresentation
    local_face::FaceIndex
    remote_face::Dict{Int,FaceIndex}
end

@inline local_entities(sf::SharedFace) = (sf.local_face,)
@inline remote_entities(sf::SharedFace) = sf.remote_face


"""
    SharedEdge <: SharedEntity

A shared edge induced by a local edge index and all remote edge indices on all remote ranks.
"""
struct SharedEdge <: SharedEntity
    unique_local_representation::EdgeRepresentation # Identify via node in grid
    local_edges::Vector{EdgeIndex}
    remote_edges::Dict{Int,Vector{EdgeIndex}}
end

@inline local_entities(se::SharedEdge) = se.local_edges
@inline remote_entities(se::SharedEdge) = se.remote_edges
