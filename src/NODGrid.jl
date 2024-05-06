"""
    NODGrid{dim,C<:AbstractCell,T<:Real} <: AbstractNODGrid{dim}

Scalable non-overlapping distributed grid. This data structure is composed of a
`local_grid` and full topological information on the process boundary, i.e.
how vertices, edges and faces are connectedted between processes.

!!! todo
    PartitionedArrays.jl ready constructor via extension
"""
mutable struct NODGrid{dim,C<:AbstractCell,T<:Real} <: AbstractNODGrid{dim}
    # Dense comminicator on the grid
    grid_comm::MPI.Comm
    # Sparse communicator along the shared vertex neighbors
    # We only need this one because the vertices induce the edge and face neighbors.
    interface_comm::MPI.Comm
    # Here we store the full local grid
    local_grid::Grid{dim,C,T}
    # Local copies of the shared entities of the form (local index, (process id in grid_comm, remote index))
    # The entities consistently contain their *Index, because faces and edges are not materialized. 
    shared_vertices::Dict{VertexRepresentation, SharedVertex}
    shared_edges::Dict{EdgeRepresentation, SharedEdge}
    shared_faces::Dict{FaceRepresentation, SharedFace}
end

"""
    global_comm(::NODGrid)

Global dense communicator of the distributed grid.
"""
@inline global_comm(dgrid::NODGrid) = dgrid.grid_comm

"""
    interface_comm(::NODGrid)

Graph communicator for shared vertices. Guaranteed to be derived from the communicator 
returned by @global_comm .
"""
@inline interface_comm(dgrid::NODGrid) = dgrid.interface_comm

"""
    global_rank(::NODGrid)

Get the rank on the global communicator of the distributed grid.
"""
@inline global_rank(dgrid::NODGrid) =  MPI.Comm_rank(global_comm(dgrid))+1

"""
    global_nranks(::NODGrid)

Get the number of ranks on the global communicator of the distributed grid.
"""
@inline global_nranks(dgrid::NODGrid) =  MPI.Comm_size(global_comm(dgrid))

"""
    NODGrid(grid_comm::MPI.Comm, grid_to_distribute::Grid{dim,C,T})

Construct a non-overlapping distributed grid from a grid with the SFC induced by the element ordering on a specified MPI communicator.
It is assumed that this function is called with exactly the same grid on each MPI process in the communicator.
"""
function NODGrid(grid_comm::MPI.Comm, grid_to_distribute::Grid{dim,C,T}, alg = PartitioningAlgorithm.SFC()) where {dim,C,T}
    grid_topology = CoverTopology(grid_to_distribute)
    nparts = MPI.Comm_size(grid_comm)
    partitioning = create_partitioning(grid_to_distribute, grid_topology, nparts, alg)
    NODGrid(grid_comm, grid_to_distribute, grid_topology, partitioning)
end

"""
    NODGrid(grid_comm::MPI.Comm, grid_to_distribute::Grid{dim,C,T}, grid_topology::CoverTopology, partitioning::Vector{<:Integer})

Construct a non-overlapping distributed grid from a grid with given topology and partitioning on a specified MPI communicator.
"""
function NODGrid(grid_comm::MPI.Comm, global_grid::Grid{dim,C,T}, grid_topology::CoverTopology, partitioning::Vector{<:Integer}) where {dim,C,T}
    n_cells_global = getncells(global_grid)
    @assert n_cells_global > 0 "Please provide a non-empty input mesh."

    partmin,partmax = extrema(partitioning)
    @assert partmin > 0
    @assert partmax <= MPI.Comm_size(grid_comm)

    my_rank = MPI.Comm_rank(grid_comm)+1

    # Start extraction of local grid
    # 1. Extract local cells
    local_cells = getcells(global_grid)[[i for i ∈ 1:n_cells_global if partitioning[i] == my_rank]]
    @assert length(local_cells) > 0 # Cannot handle empty partitions yet

    # 2. Find unique nodes
    local_node_index_set = Set{Int}()
    for cell ∈ local_cells
        for global_node_idx ∈ cell.nodes # @TODO abstraction
            push!(local_node_index_set, global_node_idx)
        end
    end

    # 3. Build a map for global to local node indices
    next_local_node_idx = 1
    global_to_local_node_map = Dict{Int,Int}()
    for global_node_idx ∈ local_node_index_set
        global_to_local_node_map[global_node_idx] = next_local_node_idx
        next_local_node_idx += 1
    end

    # 4. Extract local nodes
    local_nodes = Vector{Node{dim,T}}(undef,length(local_node_index_set))
    begin
        global_nodes = getnodes(global_grid)
        for global_node_idx ∈ local_node_index_set
            local_node_idx = global_to_local_node_map[global_node_idx]
            local_nodes[local_node_idx] = global_nodes[global_node_idx]
        end
    end

    # 5. Transform cell indices
    for local_cell_idx ∈ 1:length(local_cells)
        local_cells[local_cell_idx] = C(map(global_node_idx -> global_to_local_node_map[global_node_idx], local_cells[local_cell_idx].nodes))
    end

    # 6. Extract sets
    # @TODO deduplicate the code. We should be able to merge each of these into a macro or function.
    # We build this map now, so we avoid the communication later.
    global_to_local_cell_map = Dict{Int,Dict{Int,Int}}()
    for rank ∈ 1:MPI.Comm_size(grid_comm)
        global_to_local_cell_map[rank] = Dict{Int,Int}()
        next_local_cell_idx = 1
        for global_cell_idx ∈ 1:n_cells_global
            if partitioning[global_cell_idx] == rank
                global_to_local_cell_map[rank][global_cell_idx] = next_local_cell_idx
                next_local_cell_idx += 1
            end
        end
    end

    cellsets = Dict{String,Set{Int}}()
    for key ∈ keys(global_grid.cellsets)
        cellsets[key] = Set{Int}() # create empty set, so it does not crash during assembly
        for global_cell_idx ∈ global_grid.cellsets[key]
            if haskey(global_to_local_cell_map[my_rank], global_cell_idx)
                push!(cellsets[key], global_to_local_cell_map[my_rank][global_cell_idx])
            end
        end
    end

    nodesets = Dict{String,Set{Int}}()
    for key ∈ keys(global_grid.nodesets)
        nodesets[key] = Set{Int}() # create empty set, so it does not crash during assembly
        for global_node_idx ∈ global_grid.nodesets[key]
            if haskey(global_to_local_node_map, global_node_idx)
                push!(nodesets[key], global_to_local_node_map[global_node_idx])
            end
        end
    end

    facesets = Dict{String,Set{FaceIndex}}()
    for key ∈ keys(global_grid.facesets)
        facesets[key] = Set{FaceIndex}() # create empty set, so it does not crash during assembly
        for (global_cell_idx, i) ∈ global_grid.facesets[key]
            if haskey(global_to_local_cell_map[my_rank], global_cell_idx)
                push!(facesets[key], FaceIndex(global_to_local_cell_map[my_rank][global_cell_idx], i))
            end
        end
    end

    edgesets = Dict{String,Set{EdgeIndex}}()
    for key ∈ keys(global_grid.edgesets)
        edgesets[key] = Set{EdgeIndex}() # create empty set, so it does not crash during assembly
        for (global_cell_idx, i) ∈ global_grid.edgesets[key]
            if haskey(global_to_local_cell_map[my_rank], global_cell_idx)
                push!(edgesets[key], EdgeIndex(global_to_local_cell_map[my_rank][global_cell_idx], i))
            end
        end
    end

    vertexsets = Dict{String,Set{VertexIndex}}()
    for key ∈ keys(global_grid.vertexsets)
        vertexsets[key] = Set{VertexIndex}() # create empty set, so it does not crash during assembly
        for (global_cell_idx, i) ∈ global_grid.vertexsets[key]
            if haskey(global_to_local_cell_map[my_rank], global_cell_idx)
                push!(vertexsets[key], VertexIndex(global_to_local_cell_map[my_rank][global_cell_idx], i))
            end
        end
    end

    # Invert lookup table for fast queries
    local_to_global_cell_map = zeros(Int, length(local_cells))
    for (global_cellidx, local_cellidx) in global_to_local_cell_map[my_rank]
        local_to_global_cell_map[local_cellidx] = global_cellidx
    end
    # Debug check
    for (local_cellidx, global_cellidx) in enumerate(local_to_global_cell_map)
        @assert global_cellidx != 0 "$local_cellidx not mapped to any global element. Aborting. (R$my_rank)"
    end 

    local_grid = Grid(
        local_cells,
        local_nodes,
        cellsets=cellsets,
        nodesets=nodesets,
        facesets=facesets,
        edgesets=edgesets,
        vertexsets=vertexsets
    )

    # We use these to efficiently determine the unique vertices, faces and edges
    shared_vertices = Dict{VertexRepresentation,SharedVertex}()
    shared_edges    = Dict{EdgeRepresentation,SharedEdge}()
    shared_faces    = Dict{FaceRepresentation,SharedFace}()
    # TODO rewrite more efficiently by looping over the local boundary and check for the codim 1 entity if the global grid has an associated neighboring element
    for (local_cell_idx,local_cell) ∈ enumerate(getcells(local_grid))
        global_cell_idx = local_to_global_cell_map[local_cell_idx]
        global_cell = getcells(global_grid, global_cell_idx)

        # Vertex
        for (vi, local_vertex_node) ∈ enumerate(Ferrite.vertices(local_cell))
            # If we have already visited the vertex we can just skip
            local_vid = VertexRepresentation(local_vertex_node)
            haskey(shared_vertices, local_vid) && continue

            # Note that by construction the cells in the global and local grid share the same orientation
            cell_vertex = VertexIndex(global_cell_idx, vi)
            # Stores local shared vertex index [1, num_local_shared_vertices] -> VertexIndex in local grid
            local_vertices = Vector{VertexIndex}()
            # Stores for each neighboring rank local shared vertex index [1, num_local_shared_vertices_towards_rank] -> VertexIndex in local grid on remote
            remote_vertices = Dict{Int,Vector{VertexIndex}}()
            for global_neighbor_vertex ∈ getneighborhood(grid_topology, global_grid, cell_vertex, true)
                # Unpack VertexIndex
                (global_cell_neighbor_idx, neighbor_vi) = global_neighbor_vertex
                # Get rank of associated element
                neighbor_rank = partitioning[global_cell_neighbor_idx]
                # Store whether the neighbor is remote or not
                if neighbor_rank != my_rank
                    if !haskey(remote_vertices,neighbor_rank)
                        remote_vertices[neighbor_rank] = Vector(undef,0)
                    end
                    Ferrite.@debug println("Detected shared vertex $local_vid remote neighbor $global_neighbor_vertex on $neighbor_rank (R$my_rank)")
                    push!(remote_vertices[neighbor_rank], VertexIndex(global_to_local_cell_map[neighbor_rank][global_cell_neighbor_idx], neighbor_vi))
                else
                    Ferrite.@debug println("Detected shared vertex $local_vid local neighbor $global_neighbor_vertex (R$my_rank)")
                    push!(local_vertices, VertexIndex(global_to_local_cell_map[my_rank][global_cell_neighbor_idx], neighbor_vi))
                end
            end

            # Just store store the information if there is some actual remote neighbor
            if length(remote_vertices) > 0
                shared_vertices[local_vid] = SharedVertex(local_vid, local_vertices, remote_vertices)
            else
                # local vertex: do nothing
            end
        end

        # Face
        if dim > 1
            # If we have already visited the face we can just skip
            for (fi, local_face_nodes) ∈ enumerate(Ferrite.faces(local_cell))
                # The face should also just have one real face neighbor in the topology
                global_neighbor_faces = getneighborhood(grid_topology, global_grid, FaceIndex(global_cell_idx, fi), false)
                length(global_neighbor_faces) == 0 && continue # True boundary
                @assert length(global_neighbor_faces) == 1 "Face topology broken! (R$my_rank)"

                # If we hit the same shared face twice in a local grid, then the grid must be broken, because the shared faces must be on the boundary and hence just associated to one local cell!
                local_fid = FaceRepresentation(local_face_nodes)
                @assert !haskey(shared_faces, local_fid) "Grid topology broken. Boundary face with multiple elements attached detected."

                # Unpack face
                (global_cell_neighbor_idx, neighbor_fi) = global_neighbor_faces[1]
                neighbor_rank = partitioning[global_cell_neighbor_idx]
                if neighbor_rank != my_rank
                    # Construct local information for current and remote rank
                    Ferrite.@debug println("Detected shared face $local_fid neighbor $(global_neighbor_faces[1]) on $neighbor_rank (R$my_rank)")
                    lfi = FaceIndex(local_cell_idx, fi)
                    rfi = Dict(Pair(neighbor_rank, FaceIndex(global_to_local_cell_map[neighbor_rank][global_cell_neighbor_idx], neighbor_fi)))
                    shared_faces[local_fid] = SharedFace(local_fid, lfi, rfi)
                else
                    # local face: do nothing
                end
            end
        end

        # Edge
        if dim > 2
            for (ei, local_edge_nodes) ∈ enumerate(Ferrite.edges(local_cell))
                # If we have already visited the edge we can just skip
                local_eid = EdgeRepresentation(local_edge_nodes)
                haskey(shared_edges, local_eid) && continue

                # Note that by construction the cells in the global and local grid share the same orientation
                cell_edge = EdgeIndex(global_cell_idx, ei)
                # Stores local shared edge index [1, num_local_shared_edges] -> EdgeIndex in local grid
                local_edges = Vector{EdgeIndex}()
                # Stores for each neighboring rank local shared edge index [1, num_local_shared_edges_towards_rank] -> EdgeIndex in local grid on remote
                remote_edges = Dict{Int,Vector{EdgeIndex}}()
                for global_neighbor_edge ∈ getneighborhood(grid_topology, global_grid, cell_edge, true)
                    # Unpack edge 
                    (global_cell_neighbor_idx, neighbor_ei) = global_neighbor_edge
                    neighbor_rank = partitioning[global_cell_neighbor_idx]
                    # Store whether the neighbor is remote or not
                    if neighbor_rank != my_rank
                        if !haskey(remote_edges, global_neighbor_edge)
                            remote_edges[neighbor_rank] = Vector(undef,0)
                        end
                        Ferrite.@debug println("Detected shared edge $local_eid remote neighbor $global_neighbor_edge on $neighbor_rank (R$my_rank)")
                        push!(remote_edges[neighbor_rank], EdgeIndex(global_to_local_cell_map[neighbor_rank][global_cell_neighbor_idx], neighbor_ei))
                    else
                        Ferrite.@debug println("Detected shared edge $local_eid local neighbor $global_neighbor_edge (R$my_rank)")
                        push!(local_edges, EdgeIndex(global_to_local_cell_map[my_rank][global_cell_neighbor_idx], neighbor_ei))
                    end
                end

                if length(remote_edges) > 0
                    shared_edges[local_eid] = SharedEdge(local_eid, local_edges, remote_edges)
                end
            end
        end
    end

    # Neighborhood graph
    neighbors_set = Set{Cint}()
    for (vi, sv) ∈ shared_vertices
        for (rank, vvi) ∈ sv.remote_vertices
            push!(neighbors_set, rank)
        end
    end
    # Adjust ranks back to to C index convention
    dest = collect(neighbors_set).-1
    degree = length(dest)
    interface_comm = MPI.Dist_graph_create(grid_comm, Cint[my_rank-1], Cint[degree], Cint.(dest))

    return NODGrid(grid_comm,interface_comm,local_grid,shared_vertices,shared_edges,shared_faces)
end


# Here we define the entity ownership by the process sharing an entity with lowest rank in the grid communicator.
function compute_owner(dgrid::AbstractNODGrid, shared_entity::SharedEntity)
    my_rank = MPI.Comm_rank(global_comm(dgrid))+1 # Shift rank up by 1 to match Julia's indexing convention
    return minimum([my_rank; [remote_rank for (remote_rank, _) ∈ remote_entities(shared_entity)]])
end

"""
    generate_nod_grid(comm::MPI.Comm, args...)

Helper to directly generate non-overlapping distributed grids, designed to replace the call to [`generate_grid`](@ref) for use in distributed environments.
"""
function generate_nod_grid(comm::MPI.Comm, args...; partitioning_alg=PartitioningAlgorithm.SFC())
    full_grid = generate_grid(args...)
    return NODGrid(comm, full_grid, partitioning_alg)
end
