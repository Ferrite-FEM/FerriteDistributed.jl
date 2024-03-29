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
    shared_vertices::Dict{VertexIndex,SharedVertex}
    shared_edges::Dict{EdgeIndex,SharedEdge}
    shared_faces::Dict{FaceIndex,SharedFace}
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
function NODGrid(grid_comm::MPI.Comm, grid_to_distribute::Grid{dim,C,T}, grid_topology::CoverTopology, partitioning::Vector{<:Integer}) where {dim,C,T}
    n_cells_global = getncells(grid_to_distribute)
    @assert n_cells_global > 0 "Please provide a non-empty input mesh."

    partmin,partmax = extrema(partitioning)
    @assert partmin > 0
    @assert partmax <= MPI.Comm_size(grid_comm)

    my_rank = MPI.Comm_rank(grid_comm)+1

    # Start extraction of local grid
    # 1. Extract local cells
    local_cells = getcells(grid_to_distribute)[[i for i ∈ 1:n_cells_global if partitioning[i] == my_rank]]
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
        global_nodes = getnodes(grid_to_distribute)
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
    for key ∈ keys(grid_to_distribute.cellsets)
        cellsets[key] = Set{Int}() # create empty set, so it does not crash during assembly
        for global_cell_idx ∈ grid_to_distribute.cellsets[key]
            if haskey(global_to_local_cell_map[my_rank], global_cell_idx)
                push!(cellsets[key], global_to_local_cell_map[my_rank][global_cell_idx])
            end
        end
    end

    nodesets = Dict{String,Set{Int}}()
    for key ∈ keys(grid_to_distribute.nodesets)
        nodesets[key] = Set{Int}() # create empty set, so it does not crash during assembly
        for global_node_idx ∈ grid_to_distribute.nodesets[key]
            if haskey(global_to_local_node_map, global_node_idx)
                push!(nodesets[key], global_to_local_node_map[global_node_idx])
            end
        end
    end

    facesets = Dict{String,Set{FaceIndex}}()
    for key ∈ keys(grid_to_distribute.facesets)
        facesets[key] = Set{FaceIndex}() # create empty set, so it does not crash during assembly
        for (global_cell_idx, i) ∈ grid_to_distribute.facesets[key]
            if haskey(global_to_local_cell_map[my_rank], global_cell_idx)
                push!(facesets[key], FaceIndex(global_to_local_cell_map[my_rank][global_cell_idx], i))
            end
        end
    end

    edgesets = Dict{String,Set{EdgeIndex}}()
    for key ∈ keys(grid_to_distribute.edgesets)
        edgesets[key] = Set{EdgeIndex}() # create empty set, so it does not crash during assembly
        for (global_cell_idx, i) ∈ grid_to_distribute.edgesets[key]
            if haskey(global_to_local_cell_map[my_rank], global_cell_idx)
                push!(edgesets[key], EdgeIndex(global_to_local_cell_map[my_rank][global_cell_idx], i))
            end
        end
    end

    vertexsets = Dict{String,Set{VertexIndex}}()
    for key ∈ keys(grid_to_distribute.vertexsets)
        vertexsets[key] = Set{VertexIndex}() # create empty set, so it does not crash during assembly
        for (global_cell_idx, i) ∈ grid_to_distribute.vertexsets[key]
            if haskey(global_to_local_cell_map[my_rank], global_cell_idx)
                push!(vertexsets[key], VertexIndex(global_to_local_cell_map[my_rank][global_cell_idx], i))
            end
        end
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

    shared_vertices = Dict{VertexIndex,SharedVertex}()
    shared_edges = Dict{EdgeIndex,SharedEdge}()
    shared_faces = Dict{FaceIndex,SharedFace}()
    for (global_cell_idx,global_cell) ∈ enumerate(getcells(grid_to_distribute))
        if partitioning[global_cell_idx] == my_rank
            # Vertex
            for (i, _) ∈ enumerate(Ferrite.vertices(global_cell))
                cell_vertex = VertexIndex(global_cell_idx, i)
                remote_vertices = Dict{Int,Vector{VertexIndex}}()
                for other_vertex ∈ getneighborhood(grid_topology, grid_to_distribute, cell_vertex, true)
                    (global_cell_neighbor_idx, j) = other_vertex
                    other_rank = partitioning[global_cell_neighbor_idx]
                    if other_rank != my_rank
                        if Ferrite.toglobal(grid_to_distribute,cell_vertex) == Ferrite.toglobal(grid_to_distribute,other_vertex)
                            if !haskey(remote_vertices,other_rank)
                                remote_vertices[other_rank] = Vector(undef,0)
                            end
                            Ferrite.@debug println("Detected shared vertex $cell_vertex neighbor $other_vertex (R$my_rank)")
                            push!(remote_vertices[other_rank], VertexIndex(global_to_local_cell_map[other_rank][global_cell_neighbor_idx], j))
                        end
                    end
                end

                if length(remote_vertices) > 0
                    idx = VertexIndex(global_to_local_cell_map[my_rank][global_cell_idx], i)
                    shared_vertices[idx] = SharedVertex(idx, remote_vertices)
                end
            end

            # Face
            if dim > 1
                for (i, _) ∈ enumerate(Ferrite.faces(global_cell))
                    cell_face = FaceIndex(global_cell_idx, i)
                    remote_faces = Dict{Int,Vector{FaceIndex}}()
                    for other_face ∈ getneighborhood(grid_topology, grid_to_distribute, cell_face, true)
                        (global_cell_neighbor_idx, j) = other_face
                        other_rank = partitioning[global_cell_neighbor_idx]
                        if other_rank != my_rank
                            if Ferrite.toglobal(grid_to_distribute,cell_face) == Ferrite.toglobal(grid_to_distribute,other_face)
                                if !haskey(remote_faces,other_rank)
                                    remote_faces[other_rank] = Vector(undef,0)
                                end
                                Ferrite.@debug println("Detected shared face $cell_face neighbor $other_face (R$my_rank)")
                                push!(remote_faces[other_rank], FaceIndex(global_to_local_cell_map[other_rank][global_cell_neighbor_idx], j))
                            end
                        end
                    end

                    if length(remote_faces) > 0
                        idx = FaceIndex(global_to_local_cell_map[my_rank][global_cell_idx], i)
                        shared_faces[idx] = SharedFace(idx, remote_faces)
                    end
                end
            end

            # Edge
            if dim > 2
                for (i, _) ∈ enumerate(Ferrite.edges(global_cell))
                    cell_edge = EdgeIndex(global_cell_idx, i)
                    remote_edges = Dict{Int,Vector{EdgeIndex}}()
                    for other_edge ∈ getneighborhood(grid_topology, grid_to_distribute, cell_edge, true)
                        (global_cell_neighbor_idx, j) = other_edge
                        other_rank = partitioning[global_cell_neighbor_idx]
                        if other_rank != my_rank
                            if Ferrite.toglobal(grid_to_distribute,cell_edge) == Ferrite.toglobal(grid_to_distribute,other_edge)
                                if !haskey(remote_edges,other_edge)
                                    remote_edges[other_rank] = Vector(undef,0)
                                end
                                Ferrite.@debug println("Detected shared edge $cell_edge neighbor $other_edge (R$my_rank)")
                                push!(remote_edges[other_rank], EdgeIndex(global_to_local_cell_map[other_rank][global_cell_neighbor_idx], j))
                            end
                        end
                    end

                    if length(remote_edges) > 0
                        idx = EdgeIndex(global_to_local_cell_map[my_rank][global_cell_idx], i)
                        shared_edges[idx] = SharedEdge(idx, remote_edges)
                    end
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
