"""
    AbstractNODGrid{sdim} <: Ferrite.AbstractGrid{sdim}

Supertype for the non-overlapping grid interface.
"""
abstract type AbstractNODGrid{sdim} <: Ferrite.AbstractGrid{sdim} end

"""
    compute_owner(::AbstractNODGrid, ::SharedEntity)

Compute which rank (1-based index) owns the corresponding entity. 
"""
compute_owner(::AbstractNODGrid, ::SharedEntity)


# TODO we might want to replace these with actual primitives, because these are MPI specific.
"""
    global_comm(::AbstractNODGrid)

Get a communicator handle for the full grid.
"""
global_comm(::AbstractNODGrid)

"""
    interface_comm(::AbstractNODGrid)

Get a graph comminicator handle for the process boundary.
"""
interface_comm(::AbstractNODGrid)

"""
    global_rank(::AbstractNODGrid)

Get the 1-based rank of the global communicator of the grid.
"""
global_rank(::AbstractNODGrid)

"""
    global_nranks(::AbstractNODGrid)

Get the number of ranks of the global communicator of the grid.
"""
global_nranks(::AbstractNODGrid)

"""
    get_shared_vertices(::AbstractNODGrid)

Get an interable over the shared vertices.
"""
@inline get_shared_vertices(dgrid::AbstractNODGrid) = values(dgrid.shared_vertices)

"""
get_shared_edges(::AbstractNODGrid)

Get an interable over the shared edges.
"""
@inline get_shared_edges(dgrid::AbstractNODGrid) = values(dgrid.shared_edges)

"""
    get_shared_faces(::AbstractNODGrid)

Get an interable over the shared faces.
"""
@inline get_shared_faces(dgrid::AbstractNODGrid) = values(dgrid.shared_faces)

"""
    get_shared_vertex(::AbstractNODGrid, ::VertexIndex)

Get the shared vertex associated to the VertexIndex, if it exists.
"""
@inline get_shared_vertex(dgrid::AbstractNODGrid, vi::VertexIndex) = dgrid.shared_vertices[vi]

"""
    get_shared_edge(::AbstractNODGrid, ::EdgeIndex)

Get the shared edge associated to the EdgeIndex, if it exists.
"""
@inline get_shared_edge(dgrid::AbstractNODGrid, ei::EdgeIndex) = dgrid.shared_edges[ei]

"""
    get_shared_face(::AbstractNODGrid, ::FaceIndex)

Get the shared edge associated to the FaceIndex, if it exists.
"""
@inline get_shared_face(dgrid::AbstractNODGrid, fi::FaceIndex) = dgrid.shared_faces[fi]

"""
    is_shared_vertex(::AbstractNODGrid, ::VertexIndex)

Check if a VertexIndex is a shared vertex.
"""
@inline is_shared_vertex(dgrid::AbstractNODGrid, vi::VertexIndex) = haskey(dgrid.shared_vertices, vi)

"""
    is_shared_edge(::AbstractNODGrid, ::EdgeIndex)

Check if a EdgeIndex is a shared edge.
"""
@inline is_shared_edge(dgrid::AbstractNODGrid, ei::EdgeIndex) = haskey(dgrid.shared_edges, ei)

"""
    is_shared_face(::AbstractNODGrid, ::FaceIndex)

Check if a FaceIndex is a shared face.
"""
@inline is_shared_face(dgrid::AbstractNODGrid, fi::FaceIndex) = haskey(dgrid.shared_faces, fi)

"""
    getlocalgrid(::AbstractNODGrid)

Get the representative local grid containing only a vanilla local grid.
"""
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

"""
    extract_local_part!(u_ferrite::Vector, u_extension, dh::Ferrite.AbstractDofHandler)

Entry point for extensions to register a transfer function translating the solution representation of the extension 
to a Ferrite compatible vector.
"""
extract_local_part!(u_ferrite::Vector, u_extension, dh::Ferrite.AbstractDofHandler) = error("Not implemented.")
