"""
    CoverTopology(grid::AbstractGrid)

`CoverTopology` stores neighborhood information of a grid using cover semantics.
Unlike [`ExclusiveTopology`](@ref), neighbors are classified at all matching levels:
a face-sharing pair also appears in edge and vertex neighbor lists.

Uses [`ArrayOfVectorViews`] for compact, cache-friendly storage (matching Ferrite's
`ExclusiveTopology` data layout).
"""
struct CoverTopology <: AbstractTopology
    # maps a global vertex id to all cells containing the vertex
    vertex_to_cell::Ferrite.CollectionsOfViews.ArrayOfVectorViews{Int, 1}
    # index of the vector = cell id -> all other connected cells (plain Int, not CellIndex)
    cell_neighbor::Ferrite.CollectionsOfViews.ArrayOfVectorViews{Int, 1}
    # face_face_neighbor[cellid, local_face_id] -> connected face entities
    face_face_neighbor::Ferrite.CollectionsOfViews.ArrayOfVectorViews{FaceIndex, 2}
    # vertex_vertex_neighbor[cellid, local_vertex_id] -> connected vertex entities
    vertex_vertex_neighbor::Ferrite.CollectionsOfViews.ArrayOfVectorViews{VertexIndex, 2}
    # edge_edge_neighbor[cellid, local_edge_id] -> connected edge entities
    edge_edge_neighbor::Ferrite.CollectionsOfViews.ArrayOfVectorViews{EdgeIndex, 2}
    # list of unique faces in the grid given as FaceIndex
    face_skeleton::Union{Vector{FaceIndex}, Nothing}
end

function Base.show(io::IO, ::MIME"text/plain", topology::CoverTopology)
    println(io, "CoverTopology\n")
    print(io, "  Vertex neighbors: $(size(topology.vertex_vertex_neighbor))\n")
    print(io, "  Face neighbors: $(size(topology.face_face_neighbor))\n")
    println(io, "  Edge neighbors: $(size(topology.edge_edge_neighbor))")
end

# Cover-specific: find ALL matching edges between two cells (not just the first).
function _add_all_edge_neighbors_cover!(edge_buf::Ferrite.ConstructionBuffer, cell::AbstractCell, cell_id::Int, cell_neighbor::AbstractCell, cell_neighbor_id::Int)
    for (lei, edge) in enumerate(Ferrite.edges(cell))
        uniqueedge = Ferrite.sortedge_fast(edge)
        for (lei2, edge_neighbor) in enumerate(Ferrite.edges(cell_neighbor))
            uniqueedge2 = Ferrite.sortedge_fast(edge_neighbor)
            if uniqueedge == uniqueedge2
                Ferrite.push_at_index!(edge_buf, EdgeIndex(cell_neighbor_id, lei2), cell_id, lei)
            end
        end
    end
end

function CoverTopology(grid::Ferrite.AbstractGrid)
    cells = Ferrite.getcells(grid)
    nnodes = Ferrite.getnnodes(grid)
    ncells = length(cells)

    max_vertices, max_edges, max_faces = Ferrite._max_nentities_per_cell(cells)

    # Reuse Ferrite's compact builders
    vertex_to_cell = Ferrite.build_vertex_to_cell(cells; max_vertices, nnodes)
    cell_neighbor = Ferrite.build_cell_neighbor(grid, cells, vertex_to_cell; ncells)

    # Allocate ConstructionBuffers for neighbor tables
    face_buf = Ferrite.CollectionsOfViews.ConstructionBuffer(
        sizehint!(FaceIndex[], ncells * max_faces * Ferrite._getsizehint(grid, FaceIndex)),
        (ncells, max_faces), Ferrite._getsizehint(grid, FaceIndex))
    edge_buf = Ferrite.CollectionsOfViews.ConstructionBuffer(
        sizehint!(EdgeIndex[], ncells * max_edges * Ferrite._getsizehint(grid, EdgeIndex)),
        (ncells, max_edges), Ferrite._getsizehint(grid, EdgeIndex))
    vert_buf = Ferrite.CollectionsOfViews.ConstructionBuffer(
        sizehint!(VertexIndex[], ncells * max_vertices * Ferrite._getsizehint(grid, VertexIndex)),
        (ncells, max_vertices), Ferrite._getsizehint(grid, VertexIndex))

    # Cover semantics: >= thresholds (not == like ExclusiveTopology)
    for (cell_id, cell) in enumerate(cells)
        for neighbor_cell_id in cell_neighbor[cell_id]
            neighbor_cell = cells[neighbor_cell_id]
            Ferrite.getrefdim(neighbor_cell) == Ferrite.getrefdim(cell) || continue

            num_shared = Ferrite._num_shared_vertices(cell, neighbor_cell)

            if num_shared >= 1
                Ferrite._add_single_vertex_neighbor!(vert_buf, cell, cell_id, neighbor_cell, neighbor_cell_id)
            end
            if num_shared >= 2
                if Ferrite.getrefdim(cell) == 2
                    Ferrite._add_single_face_neighbor!(face_buf, cell, cell_id, neighbor_cell, neighbor_cell_id)
                elseif Ferrite.getrefdim(cell) == 3
                    _add_all_edge_neighbors_cover!(edge_buf, cell, cell_id, neighbor_cell, neighbor_cell_id)
                else
                    @error "Case not implemented."
                end
            end
            if num_shared >= 3
                Ferrite._add_single_face_neighbor!(face_buf, cell, cell_id, neighbor_cell, neighbor_cell_id)
            end
            if num_shared <= 0
                @error "Found connected elements without shared vertex... Mesh broken?"
            end
        end
    end

    face_face_neighbor = Ferrite.CollectionsOfViews.ArrayOfVectorViews(face_buf)
    edge_edge_neighbor = Ferrite.CollectionsOfViews.ArrayOfVectorViews(edge_buf)
    vertex_vertex_neighbor = Ferrite.CollectionsOfViews.ArrayOfVectorViews(vert_buf)

    return CoverTopology(vertex_to_cell, cell_neighbor, face_face_neighbor, vertex_vertex_neighbor, edge_edge_neighbor, nothing)
end

"""
    getneighborhood(top::CoverTopology, grid::AbstractGrid, cellidx::CellIndex, include_self=false)
    getneighborhood(top::CoverTopology, grid::AbstractGrid, faceidx::FaceIndex, include_self=false)
    getneighborhood(top::CoverTopology, grid::AbstractGrid, vertexidx::VertexIndex, include_self=false)
    getneighborhood(top::CoverTopology, grid::AbstractGrid, edgeidx::EdgeIndex, include_self=false)

Returns all directly connected entities of the same type using cover semantics.
If `include_self` is true, the given entity is included in the returned list.
"""
function Ferrite.getneighborhood(top::CoverTopology, grid::Ferrite.AbstractGrid, cellidx::CellIndex, include_self=false)
    patch = top.cell_neighbor[cellidx.idx]
    if include_self
        return view(push!(collect(patch), cellidx.idx), 1:(length(patch) + 1))
    else
        return patch
    end
end

function Ferrite.getneighborhood(top::CoverTopology, grid::Ferrite.AbstractGrid, faceidx::FaceIndex, include_self=false)
    neighbors = faceidx[2] <= size(top.face_face_neighbor, 2) ? top.face_face_neighbor[faceidx[1], faceidx[2]] : FaceIndex[]
    if include_self
        return view(push!(collect(neighbors), faceidx), 1:(length(neighbors) + 1))
    else
        return neighbors
    end
end

function Ferrite.getneighborhood(top::CoverTopology, grid::Ferrite.AbstractGrid, vertexidx::VertexIndex, include_self=false)
    neighbors = top.vertex_vertex_neighbor[vertexidx[1], vertexidx[2]]
    if include_self
        return view(push!(collect(neighbors), vertexidx), 1:(length(neighbors) + 1))
    else
        return neighbors
    end
end

function Ferrite.getneighborhood(top::CoverTopology, grid::Ferrite.AbstractGrid{3}, edgeidx::EdgeIndex, include_self=false)
    neighbors = top.edge_edge_neighbor[edgeidx[1], edgeidx[2]]
    if include_self
        return view(push!(collect(neighbors), edgeidx), 1:(length(neighbors) + 1))
    else
        return neighbors
    end
end
