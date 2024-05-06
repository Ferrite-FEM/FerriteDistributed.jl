"""
    VertexRepresentation

We can identify an edge uniquely by the sorted node numbers associated with the end points.
"""
struct VertexRepresentation
    node::Int
end

"""
    EdgeRepresentation

We can identify an edge uniquely by the sorted node numbers associated with the end points.
"""
struct EdgeRepresentation
    a::Int
    b::Int
end

function EdgeRepresentation(ab::Tuple{Int,Int}) 
    se = Ferrite.sortedge_fast(ab)
    return EdgeRepresentation(se[1], se[2])
end

"""
    FaceRepresentation

We can identify a face uniquely by 3 sorted node numbers associated with the vertices.
3 points are sufficient, because 3 (non-aligned) points can uniquely describe a surface in 3D.
"""
struct FaceRepresentation
    a::Int
    b::Int
    c::Int
end


function FaceRepresentation(ab::Tuple{Int,Int})
    @warn "Fixme after https://github.com/Ferrite-FEM/Ferrite.jl/pull/789" maxlog=1
    sf = Ferrite.sortface_fast(ab)
    return FaceRepresentation(-1, sf[1], sf[2])
end

function FaceRepresentation(abc::Union{Tuple{Int,Int,Int}, Tuple{Int,Int,Int,Int}}) 
    sf = Ferrite.sortface_fast(abc)
    return FaceRepresentation(sf[1], sf[2], sf[3])
end


"""
    Entity

Supertype for geometric entities.
"""
abstract type Entity end

# !!!THE STUFF BELOW IS CURRENTLY UNUSED!!!

"""
    Vertex <: Entity

A shared vertex induced by a local vertex index and all remote vertex indices on all remote ranks.
"""
struct Vertex <: Entity
    unique_local_representation::VertexRepresentation # Identify via node in grid
    local_vertices::Vector{VertexIndex}
end

"""
    Face <: Entity

A face induced by a local face index and all remote face indices on all remote ranks.
"""
struct Face <: Entity
    unique_representation::FaceRepresentation
    local_faces::Pair{FaceIndex,FaceIndex}
end

"""
    Edge <: Entity

An edge induced by a local edge index and all remote edge indices on all remote ranks.
"""
struct Edge <: Entity
    unique_representation::EdgeRepresentation # Identify via node in grid
    local_edges::Vector{EdgeIndex}
end

