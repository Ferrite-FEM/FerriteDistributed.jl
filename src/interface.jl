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
@inline Ferrite.getspatialdim(::AbstractNODGrid{sdim}) where {sdim} = sdim
@inline Ferrite.getnodes(dgrid::AbstractNODGrid) = getnodes(getlocalgrid(dgrid))
@inline Ferrite.getnodes(dgrid::AbstractNODGrid, v::Union{Int, Vector{Int}}) = getnodes(getlocalgrid(dgrid), v)
@inline Ferrite.getnodes(dgrid::AbstractNODGrid, setname::String) = getnodes(getlocalgrid(dgrid), setname)
@inline Ferrite.getnnodes(dgrid::AbstractNODGrid) = getnnodes(getlocalgrid(dgrid))

@inline Ferrite.getcells(dgrid::AbstractNODGrid) = getcells(getlocalgrid(dgrid))
@inline Ferrite.getcells(dgrid::AbstractNODGrid, v::Union{Int, Vector{Int}}) = getcells(getlocalgrid(dgrid),v)
@inline Ferrite.getcells(dgrid::AbstractNODGrid, setname::String) = getcells(getlocalgrid(dgrid),setname)
"Returns the number of cells in the `<:AbstractNODGrid`."
@inline Ferrite.getncells(dgrid::AbstractNODGrid) = getncells(getlocalgrid(dgrid))
"Returns the celltype of the `<:AbstractNODGrid`."
@inline Ferrite.getcelltype(dgrid::AbstractNODGrid) = eltype(getcells(getlocalgrid(dgrid)))
@inline Ferrite.getcelltype(dgrid::AbstractNODGrid, i::Int) = typeof(getcells(getlocalgrid(dgrid),i))

@inline Ferrite.getcellset(dgrid::AbstractNODGrid, setname::String) = getcellset(getlocalgrid(dgrid), setname)
@inline Ferrite.getcellsets(dgrid::AbstractNODGrid) = getcellsets(getlocalgrid(dgrid))

@inline Ferrite.getnodeset(dgrid::AbstractNODGrid, setname::String) = getnodeset(getlocalgrid(dgrid), setname)
@inline Ferrite.getnodesets(dgrid::AbstractNODGrid) = getnodesets(getlocalgrid(dgrid))

@inline Ferrite.getfacetset(dgrid::AbstractNODGrid, setname::String) = getfacetset(getlocalgrid(dgrid), setname)
@inline Ferrite.getfacetsets(dgrid::AbstractNODGrid) = getfacetsets(getlocalgrid(dgrid))

@inline Ferrite.getvertexset(dgrid::AbstractNODGrid, setname::String) = getvertexset(getlocalgrid(dgrid), setname)
@inline Ferrite.getvertexsets(dgrid::AbstractNODGrid) = getvertexsets(getlocalgrid(dgrid))

"""
    extract_local_part!(u_ferrite::Vector, u_extension, dh::Ferrite.AbstractDofHandler)

Entry point for extensions to register a transfer function translating the solution representation of the extension 
to a Ferrite compatible vector.
"""
extract_local_part!(u_ferrite::Vector, u_extension, dh::Ferrite.AbstractDofHandler) = error("Not implemented.")

"""
    gather_dof_values!(u_local::AbstractVector, u, dh::Ferrite.AbstractDofHandler) -> u_local

Gather local values from a distributed vector `u` into `u_local`, ordered by `dh`'s local
DOF indices. After calling this, `u_local[celldofs(cell)]` gives the correct element DOF
values. Extensions should implement this for their distributed vector types.
"""
function gather_dof_values! end

"""
    gather_dof_values(u, dh::Ferrite.AbstractDofHandler) -> Vector

Allocating version of [`gather_dof_values!`](@ref). Returns a `Vector` of local DOF values
ordered by `dh`'s local DOF indices, so `result[celldofs(cell)]` gives the correct element
DOF values.
"""
function gather_dof_values end
