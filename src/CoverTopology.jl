using SparseArrays

"""
    CoverTopology(cells::Vector{C}) where C <: AbstractCell

`CoverTopology` stores the intuitive neighborhood information of a grid. Here the 
neighborhood is a set of similar entities which fully cover each other.

!!! TODO move to Ferrite core?
"""
struct CoverTopology <: AbstractTopology
    # maps a global vertex id to all cells containing the vertex
    vertex_to_cell::Vector{Set{Int}}
    # index of the vector = cell id ->  all other connected cells
    cell_neighbor::Vector{EntityNeighborhood{CellIndex}}
    # face_face_neighbor[cellid,local_face_id] -> exclusive connected entities (not restricted to one entity)
    face_face_neighbor::Matrix{EntityNeighborhood{FaceIndex}}
    # vertex_vertex_neighbor[cellid,local_vertex_id] -> exclusive connected entities to the given vertex
    vertex_vertex_neighbor::Matrix{EntityNeighborhood{VertexIndex}}
    # edge_edge_neighbor[cellid,local_edge_id] -> exclusive connected entities of the given edge
    edge_edge_neighbor::Matrix{EntityNeighborhood{EdgeIndex}}
    # list of unique faces in the grid given as FaceIndex
    face_skeleton::Union{Vector{FaceIndex}, Nothing}
end

function Base.show(io::IO, ::MIME"text/plain", topology::CoverTopology)
    println(io, "CoverTopology\n")
    print(io, "  Vertex neighbors: $(size(topology.vertex_vertex_neighbor))\n")
    print(io, "  Face neighbors: $(size(topology.face_face_neighbor))\n")
    println(io, "  Edge neighbors: $(size(topology.edge_edge_neighbor))")
end

function _add_all_edge_neighbors!(edge_table, cell::C1, cell_id, cell_neighbor::C2, cell_neighbor_id) where {C1, C2}
    for (lei, edge) ∈ enumerate(edges(cell))
        uniqueedge = Ferrite.sortedge_fast(edge)
        for (lei2, edge_neighbor) ∈ enumerate(edges(cell_neighbor))
            uniqueedge2 = Ferrite.sortedge_fast(edge_neighbor)
            if uniqueedge == uniqueedge2
                push!(edge_table[cell_id, lei].neighbor_info, EdgeIndex(cell_neighbor_id, lei2))
            end
        end
    end
end

function _cover_topology_ctor(cells::Vector{C}, vertex_cell_table::Array{Set{Int}}, vertex_table, face_table, edge_table, cell_neighbor_table) where C <: AbstractCell
    for (cell_id, cell) in enumerate(cells)
        # Gather all cells which are connected via vertices
        cell_neighbor_ids = Set{Int}()
        for vertex ∈ vertices(cell)
            for vertex_cell_id ∈ vertex_cell_table[vertex]
                if vertex_cell_id != cell_id
                    push!(cell_neighbor_ids, vertex_cell_id)
                end
            end
        end
        cell_neighbor_table[cell_id] = EntityNeighborhood(CellIndex.(collect(cell_neighbor_ids)))

        # Any of the neighbors is now sorted in the respective categories
        for cell_neighbor_id ∈ cell_neighbor_ids
            # Buffer neighbor
            cell_neighbor = cells[cell_neighbor_id]
            # TODO handle mixed-dimensional case
            getdim(cell_neighbor) == getdim(cell) || continue

            num_shared_vertices = Ferrite._num_shared_vertices(cell, cell_neighbor)

            # Simplest case: Only one vertex is shared => Vertex neighbor
            if num_shared_vertices >= 1
                for (lvi, vertex) ∈ enumerate(vertices(cell))
                    for (lvi2, vertex_neighbor) ∈ enumerate(vertices(cell_neighbor))
                        if vertex_neighbor == vertex
                            push!(vertex_table[cell_id, lvi].neighbor_info, VertexIndex(cell_neighbor_id, lvi2))
                            break
                        end
                    end
                end
            end
            # Shared path
            if num_shared_vertices >= 2
                if getdim(cell) == 2
                    Ferrite._add_single_face_neighbor!(face_table, cell, cell_id, cell_neighbor, cell_neighbor_id)
                elseif getdim(cell) == 3
                    _add_all_edge_neighbors!(edge_table, cell, cell_id, cell_neighbor, cell_neighbor_id)
                else
                    @error "Case not implemented."
                end
            end
            # Shared surface
            if num_shared_vertices >= 3
                Ferrite._add_single_face_neighbor!(face_table, cell, cell_id, cell_neighbor, cell_neighbor_id)
            end
            # Broken mesh?
            if num_shared_vertices <= 0
                @error "Found connected elements without shared vertex... Mesh broken?"
            end
        end
    end
end

"""
"""
function CoverTopology(cells::Vector{C}) where C <: Ferrite.AbstractCell
    # TODO REFACTOR ME BEGIN - Redudancy with ExclusiveTopology
    # Setup the cell to vertex table
    cell_vertices_table = vertices.(cells) #needs generic interface for <: AbstractCell
    vertex_cell_table = Set{Int}[Set{Int}() for _ ∈ 1:maximum(maximum.(cell_vertices_table))]

    # Setup vertex to cell connectivity by flipping the cell to vertex table
    for (cellid, cell_vertices) in enumerate(cell_vertices_table)
        for vertex in cell_vertices
            push!(vertex_cell_table[vertex], cellid)
        end
    end

    # Compute correct matrix size
    celltype = eltype(cells)
    max_vertices = 0
    max_faces = 0
    max_edges = 0
    if isconcretetype(celltype)
        dim = getdim(cells[1])

        max_vertices = nvertices(cells[1])
        dim > 1 && (max_faces = nfaces(cells[1]))
        dim > 2 && (max_edges = nedges(cells[1]))
    else
        celltypes = Set(typeof.(cells))
        for celltype in celltypes
            celltypeidx = findfirst(x->typeof(x)==celltype,cells)
            dim = getdim(cells[celltypeidx])

            max_vertices = max(max_vertices,nvertices(cells[celltypeidx]))
            dim > 1 && (max_faces = max(max_faces, nfaces(cells[celltypeidx])))
            dim > 2 && (max_edges = max(max_edges, nedges(cells[celltypeidx])))
        end
    end

    # Setup matrices
    vertex_table = Matrix{EntityNeighborhood{VertexIndex}}(undef, length(cells), max_vertices)
    for j = 1:size(vertex_table,2)
        for i = 1:size(vertex_table,1)
            vertex_table[i,j] = EntityNeighborhood{VertexIndex}(VertexIndex[])
        end
    end
    face_table   = Matrix{EntityNeighborhood{FaceIndex}}(undef, length(cells), max_faces)
    for j = 1:size(face_table,2)
        for i = 1:size(face_table,1)
            face_table[i,j] = EntityNeighborhood{FaceIndex}(FaceIndex[])
        end
    end
    edge_table   = Matrix{EntityNeighborhood{EdgeIndex}}(undef, length(cells), max_edges)
    for j = 1:size(edge_table,2)
        for i = 1:size(edge_table,1)
            edge_table[i,j] = EntityNeighborhood{EdgeIndex}(EdgeIndex[])
        end
    end
    cell_neighbor_table = Vector{EntityNeighborhood{CellIndex}}(undef, length(cells))
    # TODO REFACTOR ME END - Redudancy with ExclusiveTopology
    
    _cover_topology_ctor(cells, vertex_cell_table, vertex_table, face_table, edge_table, cell_neighbor_table)
    
    return CoverTopology(vertex_cell_table,cell_neighbor_table,face_table,vertex_table,edge_table,nothing)
end

"""
    getneighborhood(top::CoverTopology, grid::AbstractGrid, cellidx::CellIndex, include_self=false)
    getneighborhood(top::CoverTopology, grid::AbstractGrid, faceidx::FaceIndex, include_self=false)
    getneighborhood(top::CoverTopology, grid::AbstractGrid, vertexidx::VertexIndex, include_self=false)
    getneighborhood(top::CoverTopology, grid::AbstractGrid, edgeidx::EdgeIndex, include_self=false)

Returns all directly connected entities of the same type, i.e. calling the function with a `VertexIndex` will return
a list of directly connected vertices (connected via face/edge). If `include_self` is true, the given `*Index` is included
in the returned list.

!!! warning
    This feature is highly experimental and very likely subjected to interface changes in the future.
"""
function Ferrite.getneighborhood(top::CoverTopology, grid::Ferrite.AbstractGrid, cellidx::CellIndex, include_self=false)
    patch = getcells(top.cell_neighbor[cellidx.idx])
    if include_self
        return [patch; cellidx.idx]
    else
        return patch
    end
end

function Ferrite.getneighborhood(top::CoverTopology, grid::Ferrite.AbstractGrid, faceidx::FaceIndex, include_self=false)
    data = faceidx[2] <= size(top.face_face_neighbor, 2) ? top.face_face_neighbor[faceidx[1],faceidx[2]].neighbor_info : []
    if include_self
        return [data; faceidx]
    else
        return data
    end
end

function Ferrite.getneighborhood(top::CoverTopology, grid::Ferrite.AbstractGrid, vertexidx::VertexIndex, include_self=false)
    if include_self
        return [top.vertex_vertex_neighbor[vertexidx[1], vertexidx[2]].neighbor_info; vertexidx]
    else
        return top.vertex_vertex_neighbor[vertexidx[1], vertexidx[2]].neighbor_info
    end
end

function Ferrite.getneighborhood(top::CoverTopology, grid::Ferrite.AbstractGrid{3}, edgeidx::EdgeIndex, include_self=false)
    if include_self
        return [top.edge_edge_neighbor[edgeidx[1], edgeidx[2]].neighbor_info; edgeidx]
    else
        return top.edge_edge_neighbor[edgeidx[1], edgeidx[2]].neighbor_info
    end
end


CoverTopology(grid::Ferrite.AbstractGrid) = CoverTopology(getcells(grid))
