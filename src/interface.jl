"""
Supertype for the non-overlapping grid interface.
"""
abstract type AbstractNODGrid{sdim} <: Ferrite.AbstractGrid{sdim} end

"""
"""
compute_owner(::AbstractNODGrid, ::SharedEntity)


# TODO we might want to replace these with actual primitives, because these are MPI specific.
"""
"""
global_comm(::AbstractNODGrid)

"""
"""
interface_comm(::AbstractNODGrid)

"""
"""
global_rank(::AbstractNODGrid)

"""
"""
global_nranks(::AbstractNODGrid)

"""
"""
@inline get_shared_vertices(dgrid::AbstractNODGrid) = values(dgrid.shared_vertices)
@inline get_shared_edges(dgrid::AbstractNODGrid) = values(dgrid.shared_edges)
@inline get_shared_faces(dgrid::AbstractNODGrid) = values(dgrid.shared_faces)

@inline get_shared_vertex(dgrid::AbstractNODGrid, vi::VertexIndex) = dgrid.shared_vertices[vi]
@inline get_shared_edge(dgrid::AbstractNODGrid, ei::EdgeIndex) = dgrid.shared_edges[ei]
@inline get_shared_face(dgrid::AbstractNODGrid, fi::FaceIndex) = dgrid.shared_faces[fi]

"""
"""
@inline is_shared_vertex(dgrid::AbstractNODGrid, vi::VertexIndex) = haskey(dgrid.shared_vertices, vi)
@inline is_shared_edge(dgrid::AbstractNODGrid, ei::EdgeIndex) = haskey(dgrid.shared_edges, ei)
@inline is_shared_face(dgrid::AbstractNODGrid, fi::FaceIndex) = haskey(dgrid.shared_faces, fi)

@inline getlocalgrid(dgrid::AbstractNODGrid) = dgrid.local_grid

# Ferrite.AbstractGrid interface
@inline Ferrite.getnodes(dgrid::AbstractNODGrid) = getnodes(getlocalgrid(dgrid))
@inline Ferrite.getnodes(grid::AbstractNODGrid, v::Union{Int, Vector{Int}}) = getnodes(getlocalgrid(dgrid), v)
@inline Ferrite.getnodes(grid::AbstractNODGrid, setname::String) = getnodes(getlocalgrid(dgrid), setname)
@inline Ferrite.getnnodes(dgrid::AbstractNODGrid) = getnnodes(getlocalgrid(dgrid))

@inline Ferrite.getcells(dgrid::AbstractNODGrid) = getcells(getlocalgrid(dgrid))
@inline Ferrite.getcells(dgrid::AbstractNODGrid, v::Union{Int, Vector{Int}}) = getcells(getlocalgrid(dgrid),v)
@inline Ferrite.getcells(dgrid::AbstractNODGrid, setname::String) = getcells(getlocalgrid(dgrid),setname)
"Returns the number of cells in the `<:AbstractNODGrid`."
@inline Ferrite.getncells(dgrid::AbstractNODGrid) = getncells(getlocalgrid(dgrid))
"Returns the celltype of the `<:AbstractNODGrid`."
@inline Ferrite.getcelltype(dgrid::AbstractNODGrid) = eltype(getcells(getlocalgrid(dgrid)))
@inline Ferrite.getcelltype(dgrid::AbstractNODGrid, i::Int) = typeof(getcells(getlocalgrid(dgrid),i))

@inline Ferrite.getcellset(grid::AbstractNODGrid, setname::String) = getcellset(getlocalgrid(grid), setname)
@inline Ferrite.getcellsets(grid::AbstractNODGrid) = getcellsets(getlocalgrid(grid))

@inline Ferrite.getnodeset(grid::AbstractNODGrid, setname::String) = getnodeset(getlocalgrid(grid), setname)
@inline Ferrite.getnodesets(grid::AbstractNODGrid) = getnodeset(getlocalgrid(grid), setname)

@inline Ferrite.getfaceset(grid::AbstractNODGrid, setname::String) = getfaceset(getlocalgrid(grid), setname)
@inline Ferrite.getfacesets(grid::AbstractNODGrid) = getfaceset(getlocalgrid(grid), setname)

@inline Ferrite.getedgeset(grid::AbstractNODGrid, setname::String) = getedgeset(getlocalgrid(grid), setname)
@inline Ferrite.getedgesets(grid::AbstractNODGrid) = getedgeset(getlocalgrid(grid), setname)

@inline Ferrite.getvertexset(grid::AbstractNODGrid, setname::String) = getvertexset(getlocalgrid(grid), setname)
@inline Ferrite.getvertexsets(grid::AbstractNODGrid) = getvertexset(getlocalgrid(grid), setname)
