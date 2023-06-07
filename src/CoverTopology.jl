using SparseArrays

"""
    CoverTopology(cells::Vector{C}) where C <: AbstractCell

`CoverTopology` stores the intuitive neighborhood information of a grid. Here the 
neighborhood is a set of similar entities which fully cover each other.

!!! TODO move to Ferrite core?
"""
struct CoverTopology <: Ferrite.AbstractTopology
    # maps a global vertex id to all cells containing the vertex
    vertex_to_cell::Dict{Int,Vector{Int}}
    # index of the vector = cell id ->  all other connected cells
    cell_neighbor::Vector{Ferrite.EntityNeighborhood{CellIndex}}
    # face_neighbor[cellid,local_face_id] -> exclusive connected entities (not restricted to one entity)
    face_neighbor::SparseMatrixCSC{Ferrite.EntityNeighborhood,Int}
    # vertex_neighbor[cellid,local_vertex_id] -> exclusive connected entities to the given vertex
    vertex_neighbor::SparseMatrixCSC{Ferrite.EntityNeighborhood,Int}
    # edge_neighbor[cellid,local_edge_id] -> exclusive connected entities of the given edge
    edge_neighbor::SparseMatrixCSC{Ferrite.EntityNeighborhood,Int}
    # list of unique faces in the grid given as FaceIndex
    face_skeleton::Vector{FaceIndex}
end

"""
"""
function CoverTopology(cells::Vector{C}) where C <: Ferrite.AbstractCell
    cell_vertices_table = Ferrite.vertices.(cells) #needs generic interface for <: AbstractCell
    vertex_cell_table = Dict{Int,Vector{Int}}()

    for (cellid, cell_nodes) in enumerate(cell_vertices_table)
       for node in cell_nodes
            if haskey(vertex_cell_table, node)
                push!(vertex_cell_table[node], cellid)
            else
                vertex_cell_table[node] = [cellid]
            end
        end
    end

    I_face = Int[]; J_face = Int[]; V_face = Ferrite.EntityNeighborhood[]
    I_edge = Int[]; J_edge = Int[]; V_edge = Ferrite.EntityNeighborhood[]
    I_vertex = Int[]; J_vertex = Int[]; V_vertex = Ferrite.EntityNeighborhood[]
    cell_neighbor_table = Vector{Ferrite.EntityNeighborhood{CellIndex}}(undef, length(cells))

    for (cellid, cell) in enumerate(cells)
        #cell neighborhood
        cell_neighbors = getindex.((vertex_cell_table,), cell_vertices_table[cellid]) # cell -> vertex -> cell
        cell_neighbors = unique(reduce(vcat,cell_neighbors)) # non unique list initially
        filter!(x->x!=cellid, cell_neighbors) # get rid of self neighborhood
        cell_neighbor_table[cellid] = Ferrite.EntityNeighborhood(CellIndex.(cell_neighbors))

        for neighbor in cell_neighbors
            neighbor_local_ids = findall(x->x in cell.nodes, cells[neighbor].nodes)
            cell_local_ids = findall(x->x in cells[neighbor].nodes, cell.nodes)
            # cells are connected via exactly one vertex
            if length(cell_local_ids) == 1
                Ferrite._vertex_neighbor!(V_vertex, I_vertex, J_vertex, cellid, cell, neighbor_local_ids, neighbor, cells[neighbor])
            # cells are only connected via exactly one face
            elseif length(cell_local_ids) == Ferrite.nvertices_on_face(cell, 1)
                Ferrite._face_neighbor!(V_face, I_face, J_face, cellid, cell, neighbor_local_ids, neighbor, cells[neighbor])
                # Add edges on face
                if Ferrite.getdim(cell) > 2
                    for cell_edge_nodes ∈ Ferrite.edges(cell)
                        neighbor_edge_local_ids = findall(x->x ∈ cell_edge_nodes, cells[neighbor].nodes)
                        if length(neighbor_edge_local_ids) == Ferrite.nvertices_on_edge(cell, 1)
                            Ferrite._edge_neighbor!(V_edge, I_edge, J_edge, cellid, cell, neighbor_edge_local_ids, neighbor, cells[neighbor])
                        end
                    end
                end
                # Add vertices on face
                for cell_vertex_node ∈ Ferrite.vertices(cell)
                    neighbor_vertex_local_id = findall(x->x == cell_vertex_node, cells[neighbor].nodes)
                    if length(neighbor_vertex_local_id) == 1
                        Ferrite._vertex_neighbor!(V_vertex, I_vertex, J_vertex, cellid, cell, neighbor_vertex_local_id, neighbor, cells[neighbor])
                    end
                end
            # cells are only connected via exactly one edge
            elseif Ferrite.getdim(cell) > 2 && length(cell_local_ids) == Ferrite.nvertices_on_edge(cell, 1)
                _edge_neighbor!(V_edge, I_edge, J_edge, cellid, cell, neighbor_local_ids, neighbor, cells[neighbor])
                # Add vertices on edge
                for cell_vertex_nodes ∈ Ferrite.vertices(cell)
                    neighbor_vertex_local_id = findall(x->x == cell_vertex_nodes, cells[neighbor].nodes)
                    if length(neighbor_vertex_local_id) == 1
                        Ferrite._vertex_neighbor!(V_vertex, I_vertex, J_vertex, cellid, cell, neighbor_vertex_local_id, neighbor, cells[neighbor])
                    end
                end
            end
        end
    end

    face_neighbor = sparse(I_face,J_face,V_face)
    vertex_neighbor = sparse(I_vertex,J_vertex,V_vertex)
    edge_neighbor = sparse(I_edge,J_edge,V_edge)

    # Face Skeleton
    face_skeleton_global = Set{NTuple}()
    face_skeleton_local = Vector{FaceIndex}()
    fs_length = length(face_skeleton_global)
    for (cellid,cell) in enumerate(cells)
        for (local_face_id,face) in enumerate(Ferrite.faces(cell))
            push!(face_skeleton_global, first(Ferrite.sortface(face)))
            fs_length_new = length(face_skeleton_global)
            if fs_length != fs_length_new
                push!(face_skeleton_local, FaceIndex(cellid,local_face_id))
                fs_length = fs_length_new
            end
        end
    end
    return CoverTopology(vertex_cell_table,cell_neighbor_table,face_neighbor,vertex_neighbor,edge_neighbor,face_skeleton_local)
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
    data = faceidx[2] <= size(top.face_neighbor, 2) ? top.face_neighbor[faceidx[1],faceidx[2]].neighbor_info : []
    if include_self
        return [data; faceidx]
    else
        return data
    end
end

function Ferrite.getneighborhood(top::CoverTopology, grid::Ferrite.AbstractGrid, vertexidx::VertexIndex, include_self=false)
    if include_self
        return [top.vertex_neighbor[vertexidx[1], vertexidx[2]].neighbor_info; vertexidx]
    else
        return top.vertex_neighbor[vertexidx[1], vertexidx[2]].neighbor_info
    end
end

function Ferrite.getneighborhood(top::CoverTopology, grid::Ferrite.AbstractGrid{3}, edgeidx::EdgeIndex, include_self=false)
    if include_self
        return [top.edge_neighbor[edgeidx[1], edgeidx[2]].neighbor_info; edgeidx]
    else
        return top.edge_neighbor[edgeidx[1], edgeidx[2]].neighbor_info
    end
end


CoverTopology(grid::Ferrite.AbstractGrid) = CoverTopology(getcells(grid))
