function Ferrite.create_sparsity_pattern(::Type{<:PSparseMatrix}, dh::Ferrite.AbstractDofHandler, ch::Union{ConstraintHandler,Nothing}=nothing; kwargs...)
    error("Not implemented.")
end

"""
    _nod_to_oag_perm(dh)

Compute a permutation mapping NODDofHandler local dof indices to OwnAndGhostIndices
local indices. OwnAndGhostIndices always places owned dofs first (1:n_own) then ghost
dofs (n_own+1:n), while NODDofHandler interleaves them. This permutation is needed
to correctly index into PVector/PSparseMatrix local values using dof handler indices.
"""
function _nod_to_oag_perm(dh)
    my_rank = global_rank(getglobalgrid(dh))
    nldofs = num_local_dofs(dh)
    n_own = count(==(my_rank), dh.ldof_to_rank)
    perm = Vector{Int}(undef, nldofs)
    own_i = 0
    ghost_i = 0
    for i in 1:nldofs
        if dh.ldof_to_rank[i] == my_rank
            own_i += 1
            perm[i] = own_i
        else
            ghost_i += 1
            perm[i] = n_own + ghost_i
        end
    end
    return perm
end

"""
    extract_local_part!(u_local::AbstractVector, u::PVector, dh::NODDofHandler)

Gather local values from a distributed `PVector` into `u_local`, ordered by `dh`'s local
DOF indices. After calling this, `u_local[celldofs(cell)]` gives the correct element DOF
values.
"""
function FerriteDistributed.extract_local_part!(u_local::Vector, u::PVector, dh::FerriteDistributed.NODDofHandler)
    perm = _nod_to_oag_perm(dh)
    map(local_values(u)) do u_oag
        for i in eachindex(perm)
            u_local[i] = u_oag[perm[i]]
        end
    end
    return u_local
end

"""
Simplest partitioned assembler in COO format to obtain a PSparseMatrix and a PVector.
"""
struct COOAssembler{T}
    I::Vector{Int}
    J::Vector{Int}
    V::Vector{T}

    cols
    rows
    f::PVector

    👻remotes
    dh
    perm::Vector{Int}  # NODDofHandler local index → OAG local index

    # TODO PartitionedArrays backend as additional input arg
    # TODO fix type
    function COOAssembler{T}(dh) where {T}
        ldof_to_gdof = dh.ldof_to_gdof
        ldof_to_rank = dh.ldof_to_rank
        nldofs = num_local_dofs(dh)
        ngdofs = num_global_dofs(dh)
        dgrid = getglobalgrid(dh)
        dim = Ferrite.getspatialdim(dgrid)

        I = Int[]
        J = Int[]
        V = T[]
        sizehint!(I, nldofs)
        sizehint!(J, nldofs)
        sizehint!(V, nldofs)

        # @TODO the code below can be massively simplified by introducing a ghost layer to the
        #       distributed grid, which can efficiently precompute some of the values below.
        comm = global_comm(dgrid)
        np = MPI.Comm_size(comm)
        my_rank = global_rank(dgrid)

        Ferrite.@debug println("starting assembly... (R$my_rank)")

        # Neighborhood graph
        # @TODO cleanup old code below and use graph primitives instead.
        (source_len, destination_len, _) = MPI.Dist_graph_neighbors_count(interface_comm(dgrid))
        sources = Vector{Cint}(undef, source_len)
        destinations = Vector{Cint}(undef, destination_len)
        MPI.Dist_graph_neighbors!(interface_comm(dgrid), sources, destinations)

        # Adjust to Julia index convention
        sources .+= 1
        destinations .+= 1

        Ferrite.@debug println("Neighborhood | $sources | $destinations (R$my_rank)")

        # Invert the relations to clarify the code
        source_index = Dict{Cint, Int}()
        for (i,remote_rank) ∈ enumerate(sources)
            source_index[remote_rank] = i
        end
        destination_index = Dict{Int, Cint}()
        for (i,remote_rank) ∈ enumerate(destinations)
            destination_index[remote_rank] = i
        end

        # Note: We assume a symmetric neighborhood for now... this may not be true in general.
        # neighbors = MPIData(Int32.(sources), comm, (np,))

        # Extract locally owned dofs
        ltdof_indices = ldof_to_rank.==my_rank
        ltdof_to_gdof = ldof_to_gdof[ltdof_indices]

        Ferrite.@debug println("ltdof_to_gdof $ltdof_to_gdof (R$my_rank)")
        Ferrite.@debug println("ldof_to_gdof $ldof_to_gdof (R$my_rank)")
        Ferrite.@debug println("ldof_to_rank $ldof_to_rank (R$my_rank)")

        # Process owns rows of owned dofs. The process also may write to some remote dofs,
        # which correspond to non-owned share entities. Here we construct the rows for the
        # distributed matrix.
        # We decide for row (i.e. test function) ownership, because it the image of
        # SpMV is process local.
        row_own_to_global = ldof_to_gdof[ltdof_indices]
        row_ghost_mask = .!ltdof_indices
        row_ghost_to_global = ldof_to_gdof[row_ghost_mask]
        row_ghost_to_owner = Int32.(ldof_to_rank[row_ghost_mask])
        row_g2o = Dict{Int,Int32}()
        for (g, o) in zip(ldof_to_gdof, Int32.(ldof_to_rank))
            row_g2o[g] = o
        end
        row_own = OwnIndices(ngdofs, Int32(my_rank), row_own_to_global)
        row_ghost = GhostIndices(ngdofs, row_ghost_to_global, row_ghost_to_owner)
        row_partition_indices = OwnAndGhostIndices(row_own, row_ghost, gid -> get(row_g2o, gid, Int32(0)))
        rows = MPIArray(row_partition_indices, comm, (np,))

        Ferrite.@debug println("rows done (R$my_rank)")

        # For the locally visible columns we also have to take into account that remote
        # processes will write their data in some of these, because their remotely
        # owned trial functions overlap with the locally owned test functions.
        ghost_dof_to_global = Int[]
        ghost_dof_rank = Int32[]

        # ------------ Ghost dof synchronization ----------
        # Prepare sending ghost dofs to neighbors 👻
        #@TODO move relevant parts into dof handler
        #@TODO communication can be optimized by deduplicating entries in, and compressing the following arrays
        #@TODO reorder communication by field to eliminate need for `ghost_dof_field_index_to_send`
        ghost_dof_to_send = [Int[] for i ∈ 1:destination_len] # global dof id
        ghost_rank_to_send = [Int[] for i ∈ 1:destination_len] # rank of dof
        ghost_dof_field_index_to_send = [Int[] for i ∈ 1:destination_len]
        ghost_dof_owner = [Int[] for i ∈ 1:destination_len] # corresponding owner
        ghost_dof_pivot_to_send = [Int[] for i ∈ 1:destination_len] # corresponding dof to interact with
        for (pivot_vertex, pivot_shared_vertex) ∈ dgrid.shared_vertices
            # Start by searching shared entities which are not owned
            pivot_vertex_owner_rank = compute_owner(dgrid, pivot_shared_vertex)
            pivot_cell_idx = pivot_vertex[1]

            if my_rank != pivot_vertex_owner_rank
                sender_slot = destination_index[pivot_vertex_owner_rank]

                Ferrite.@debug println("$pivot_vertex may require synchronization (R$my_rank)")
                # Note: We have to send ALL dofs on the element to the remote.
                cell_dofs_upper_bound = (pivot_cell_idx == getncells(dh.grid)) ? length(dh.cell_dofs) : dh.cell_dofs_offset[pivot_cell_idx+1]
                cell_dofs = dh.cell_dofs[dh.cell_dofs_offset[pivot_cell_idx]:cell_dofs_upper_bound]

                for (field_idx, field_name) in zip(1:num_fields(dh), getfieldnames(dh))
                    !has_vertex_dofs(dh, field_idx, pivot_vertex) && continue
                    pivot_vertex_dofs = vertex_dofs(dh, field_idx, pivot_vertex)

                    for d ∈ 1:dh.field_dims[field_idx]
                        Ferrite.@debug println("  adding dof $(pivot_vertex_dofs[d]) to ghost sync synchronization on slot $sender_slot (R$my_rank)")

                        # Extract dofs belonging to the current field
                        #cell_field_dofs = cell_dofs[dof_range(dh, field_name)]
                        #for cell_field_dof ∈ cell_field_dofs
                        for cell_dof ∈ cell_dofs
                            append!(ghost_dof_pivot_to_send[sender_slot], ldof_to_gdof[pivot_vertex_dofs[d]])
                            append!(ghost_dof_to_send[sender_slot], ldof_to_gdof[cell_dof])
                            append!(ghost_rank_to_send[sender_slot], ldof_to_rank[cell_dof])
                            append!(ghost_dof_field_index_to_send[sender_slot], field_idx)
                        end
                    end
                end
            end
        end

        if dim > 1
            for (pivot_face, pivot_shared_face) ∈ dgrid.shared_faces
                # Start by searching shared entities which are not owned
                pivot_face_owner_rank = compute_owner(dgrid, pivot_shared_face)
                pivot_cell_idx = pivot_face[1]

                if my_rank != pivot_face_owner_rank
                    sender_slot = destination_index[pivot_face_owner_rank]

                    Ferrite.@debug println("$pivot_face may require synchronization (R$my_rank)")
                    # Note: We have to send ALL dofs on the element to the remote.
                    cell_dofs_upper_bound = (pivot_cell_idx == getncells(dh.grid)) ? length(dh.cell_dofs) : dh.cell_dofs_offset[pivot_cell_idx+1]
                    cell_dofs = dh.cell_dofs[dh.cell_dofs_offset[pivot_cell_idx]:cell_dofs_upper_bound]

                    for (field_idx, field_name) in zip(1:num_fields(dh), getfieldnames(dh))
                        !has_face_dofs(dh, field_idx, pivot_face) && continue
                        pivot_face_dofs = face_dofs(dh, field_idx, pivot_face)

                        for d ∈ 1:dh.field_dims[field_idx]
                            Ferrite.@debug println("  adding dof $(pivot_face_dofs[d]) to ghost sync synchronization on slot $sender_slot (R$my_rank)")

                            # Extract dofs belonging to the current field
                            #cell_field_dofs = cell_dofs[dof_range(dh, field_name)]
                            #for cell_field_dof ∈ cell_field_dofs
                            for cell_dof ∈ cell_dofs
                                append!(ghost_dof_pivot_to_send[sender_slot], ldof_to_gdof[pivot_face_dofs[d]])
                                append!(ghost_dof_to_send[sender_slot], ldof_to_gdof[cell_dof])
                                append!(ghost_rank_to_send[sender_slot], ldof_to_rank[cell_dof])
                                append!(ghost_dof_field_index_to_send[sender_slot], field_idx)
                            end
                        end
                    end
                end
            end
        end

        if dim > 2
            for (pivot_edge, pivot_shared_edge) ∈ dgrid.shared_edges
                # Start by searching shared entities which are not owned
                pivot_edge_owner_rank = compute_owner(dgrid, pivot_shared_edge)
                pivot_cell_idx = pivot_edge[1]

                if my_rank != pivot_edge_owner_rank
                    sender_slot = destination_index[pivot_edge_owner_rank]

                    Ferrite.@debug println("$pivot_edge may require synchronization (R$my_rank)")
                    # Note: We have to send ALL dofs on the element to the remote.
                    cell_dofs_upper_bound = (pivot_cell_idx == getncells(dh.grid)) ? length(dh.cell_dofs) : dh.cell_dofs_offset[pivot_cell_idx+1]
                    cell_dofs = dh.cell_dofs[dh.cell_dofs_offset[pivot_cell_idx]:cell_dofs_upper_bound]

                    for (field_idx, field_name) in zip(1:num_fields(dh), getfieldnames(dh))
                        !has_edge_dofs(dh, field_idx, pivot_edge) && continue
                        pivot_edge_dofs = edge_dofs(dh, field_idx, pivot_edge)

                        for d ∈ 1:dh.field_dims[field_idx]
                            Ferrite.@debug println("  adding dof $(pivot_edge_dofs[d]) to ghost sync synchronization on slot $sender_slot (R$my_rank)")
                            # Extract dofs belonging to the current field
                            #cell_field_dofs = cell_dofs[dof_range(dh, field_name)]
                            #for cell_field_dof ∈ cell_field_dofs
                            for cell_dof ∈ cell_dofs
                                append!(ghost_dof_pivot_to_send[sender_slot], ldof_to_gdof[pivot_edge_dofs[d]])
                                append!(ghost_dof_to_send[sender_slot], ldof_to_gdof[cell_dof])
                                append!(ghost_rank_to_send[sender_slot], ldof_to_rank[cell_dof])
                                append!(ghost_dof_field_index_to_send[sender_slot], field_idx)
                            end
                        end
                    end
                end
            end
        end

        ghost_send_buffer_lengths = Int[length(i) for i ∈ ghost_dof_to_send]
        ghost_recv_buffer_lengths = zeros(Int, destination_len)
        MPI.Neighbor_alltoall!(UBuffer(ghost_send_buffer_lengths,1), UBuffer(ghost_recv_buffer_lengths,1), interface_comm(dgrid));
        Ferrite.@debug for (i,ghost_recv_buffer_length) ∈ enumerate(ghost_recv_buffer_lengths)
            println("receiving $ghost_recv_buffer_length ghosts from $(sources[i])  (R$my_rank)")
        end

        # Communicate ghost information 👻
        # @TODO coalesce communication
        ghost_send_buffer_dofs = reduce(vcat, ghost_dof_to_send; init=Int[])
        ghost_recv_buffer_dofs = zeros(Int, sum(ghost_recv_buffer_lengths))
        MPI.Neighbor_alltoallv!(VBuffer(ghost_send_buffer_dofs,ghost_send_buffer_lengths), VBuffer(ghost_recv_buffer_dofs,ghost_recv_buffer_lengths), interface_comm(dgrid))

        ghost_send_buffer_fields = reduce(vcat, ghost_dof_field_index_to_send; init=Int[])
        ghost_recv_buffer_fields = zeros(Int, sum(ghost_recv_buffer_lengths))
        MPI.Neighbor_alltoallv!(VBuffer(ghost_send_buffer_fields,ghost_send_buffer_lengths), VBuffer(ghost_recv_buffer_fields,ghost_recv_buffer_lengths), interface_comm(dgrid))

        ghost_send_buffer_ranks = reduce(vcat, ghost_rank_to_send; init=Int[])
        ghost_recv_buffer_ranks = zeros(Int, sum(ghost_recv_buffer_lengths))
        MPI.Neighbor_alltoallv!(VBuffer(ghost_send_buffer_ranks,ghost_send_buffer_lengths), VBuffer(ghost_recv_buffer_ranks,ghost_recv_buffer_lengths), interface_comm(dgrid))

        ghost_send_buffer_dofs_piv = reduce(vcat, ghost_dof_pivot_to_send; init=Int[])
        ghost_recv_buffer_dofs_piv = zeros(Int, sum(ghost_recv_buffer_lengths))
        MPI.Neighbor_alltoallv!(VBuffer(ghost_send_buffer_dofs_piv,ghost_send_buffer_lengths), VBuffer(ghost_recv_buffer_dofs_piv,ghost_recv_buffer_lengths), interface_comm(dgrid))

        # Reconstruct source ranks
        ghost_recv_buffer_source_ranks = Int[]
        for (source_idx, recv_len) ∈ enumerate(ghost_recv_buffer_lengths)
            append!(ghost_recv_buffer_source_ranks, ones(recv_len)*sources[source_idx])
        end

        Ferrite.@debug println("received $ghost_recv_buffer_dofs with owners $ghost_recv_buffer_ranks (R$my_rank)")

        unique_ghosts_dr = sort(unique(first,zip(ghost_recv_buffer_dofs,ghost_recv_buffer_ranks)))
        # unzip manually and make sure we do not add duplicate entries to our columns
        for (dof,rank) ∈ unique_ghosts_dr
            if rank != my_rank && dof ∉ ldof_to_gdof
                push!(ghost_dof_to_global, dof)
                push!(ghost_dof_rank, rank)
            end
        end

        # ------------- Construct rows and cols of distributed matrix --------
        all_local_cols = Int[ldof_to_gdof; ghost_dof_to_global]
        all_local_col_ranks = Int32[ldof_to_rank; ghost_dof_rank]
        Ferrite.@debug println("all_local_cols $all_local_cols (R$my_rank)")
        Ferrite.@debug println("all_local_col_ranks $all_local_col_ranks (R$my_rank)")

        col_own_mask = all_local_col_ranks .== my_rank
        col_own_to_global = all_local_cols[col_own_mask]
        col_ghost_to_global = all_local_cols[.!col_own_mask]
        col_ghost_to_owner = all_local_col_ranks[.!col_own_mask]
        col_g2o = Dict{Int,Int32}()
        for (g, o) in zip(all_local_cols, all_local_col_ranks)
            col_g2o[g] = o
        end
        col_own = OwnIndices(ngdofs, Int32(my_rank), col_own_to_global)
        col_ghost = GhostIndices(ngdofs, col_ghost_to_global, col_ghost_to_owner)
        col_partition_indices = OwnAndGhostIndices(col_own, col_ghost, gid -> get(col_g2o, gid, Int32(0)))
        cols = MPIArray(col_partition_indices, comm, (np,))

        Ferrite.@debug println("cols and rows constructed (R$my_rank)")
        f = pzeros(rows)
        Ferrite.@debug println("f constructed (R$my_rank)")

        👻remotes = zip(ghost_recv_buffer_dofs_piv, ghost_recv_buffer_dofs, ghost_recv_buffer_ranks,ghost_recv_buffer_fields)
        Ferrite.@debug println("👻remotes $👻remotes (R$my_rank)")

        perm = _nod_to_oag_perm(dh)

        return new(I, J, V, cols, rows, f, 👻remotes, dh, perm)
    end
end

# TODO fix type
# Ferrite.start_assemble(dh::FerriteMPI.DistributedDofHandler, _::MPIArray) = COOAssembler{Float64}(dh)
Ferrite.start_assemble(dh, _::MPIArray) = COOAssembler{Float64}(dh)

@propagate_inbounds function Ferrite.assemble!(a::COOAssembler{T}, edof::AbstractVector{Int}, Ke::AbstractMatrix{T}) where {T}
    n_dofs = length(edof)
    append!(a.V, Ke)
    @inbounds for j in 1:n_dofs
        append!(a.I, edof)
        for i in 1:n_dofs
            push!(a.J, edof[j])
        end
    end
end

@propagate_inbounds function Ferrite.assemble!(a::COOAssembler{T}, dofs::AbstractVector{Int}, Ke::AbstractMatrix{T}, fe::AbstractVector{T}) where {T}
    Ferrite.assemble!(a, dofs, Ke)
    mapped_dofs = a.perm[dofs]
    map(local_values(a.f)) do f_local
        Ferrite.assemble!(f_local, mapped_dofs, fe)
    end
end

function Ferrite.end_assemble(assembler::COOAssembler{T}) where {T}
    dgrid = getglobalgrid(assembler.dh)
    comm = global_comm(dgrid)
    np = MPI.Comm_size(comm)
    my_rank = global_rank(dgrid)

    # --------------------- Add ghost entries in IJ 👻 --------------------
    I = map(i->assembler.dh.ldof_to_gdof[i], assembler.I)
    J = map(j->assembler.dh.ldof_to_gdof[j], assembler.J)
    V = map(v->v, assembler.V)

    # Fix ghost layer 👻! Note that the locations for remote processes to write their
    # data into are missing up to this point.
    # TODO here still the interaction between fields is missing...
    for (i, (pivot_dof, global_ghost_dof, ghost_owner_rank, ghost_field_idx)) ∈ enumerate(assembler.👻remotes)
        for dᵢ ∈ 1:1#assembler.dh.field_dims[ghost_field_idx]
            for dⱼ ∈ 1:1#assembler.dh.field_dims[ghost_field_idx]
                push!(I, pivot_dof+dᵢ-1)
                push!(J, global_ghost_dof+dⱼ-1)
                push!(V, 0.0)
            end
        end
    end

    Ferrite.@debug println("I=$(I) (R$my_rank)")
    Ferrite.@debug println("J=$(J) (R$my_rank)")
    K = PartitionedArrays.psparse(
        MPIArray(I, comm, (np,)),
        MPIArray(J, comm, (np,)),
        MPIArray(V, comm, (np,)),
        assembler.rows, assembler.cols;
        split_format=Val(false)
    ) |> fetch

    # psparse already assembles K by default (assemble=Val(true))
    PartitionedArrays.assemble!(assembler.f) |> wait

    return K, assembler.f
end
