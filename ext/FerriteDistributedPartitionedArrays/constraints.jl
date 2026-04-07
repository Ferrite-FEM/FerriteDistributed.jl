function Ferrite.meandiag(K::PartitionedArrays.PSparseMatrix)
    # Get local portion of z
    z_pa = map(local_values(K)) do K_local
        z = zero(eltype(K_local))
        for i in 1:size(K_local, 1)
            z += abs(K_local[i, i])
        end
        return z;
    end
    return PartitionedArrays.sum(z_pa) / size(K, 1)
end

"""
Poor man's Dirichlet BC application for PartitionedArrays. :)
    TODO integrate with constraints.
"""
function Ferrite.apply_zero!(K::PartitionedArrays.PSparseMatrix, f::PartitionedArrays.PVector, ch::ConstraintHandler)
    perm = _nod_to_oag_perm(ch.dh)
    pdofs = perm[ch.prescribed_dofs]

    map(local_values(f)) do f_local
        f_local[pdofs] .= 0.0
    end

    map(local_values(K), local_values(f)) do K_local, f_local
        for cdof in pdofs
            K_local[cdof, :] .= 0.0
            K_local[:, cdof] .= 0.0
            K_local[cdof, cdof] = 1.0
        end
    end
end

"""
Poor man's Dirichlet BC application for PartitionedArrays. :)
    TODO integrate with constraints.
    TODO optimize.
"""
function Ferrite.apply!(K::PartitionedArrays.PSparseMatrix, f::PartitionedArrays.PVector, ch::ConstraintHandler)
    perm = _nod_to_oag_perm(ch.dh)
    pdofs = perm[ch.prescribed_dofs]

    # Start by substracting the inhomogeneous solution from the right hand side
    u_constrained = PartitionedArrays.pzeros(K.col_partition)
    map(local_values(u_constrained)) do u_local
        u_local[pdofs] .= ch.inhomogeneities
    end
    f .-= K*u_constrained

    m = Ferrite.meandiag(K)

    comm = f.index_partition.comm

    # Then fix the RHS
    map(local_values(f), f.index_partition) do f_local, partition
        # Note: RHS only non-zero for owned RHS entries
        lto = local_to_owner(partition)
        f_local[pdofs] .= ch.inhomogeneities .* map(p -> p == MPI.Comm_rank(comm)+1, lto[pdofs]) * m
    end

    # Zero out locally visible rows and columns
    map(local_values(K)) do K_local
        for cdof ∈ pdofs
            K_local[cdof, :] .= 0.0
            K_local[:, cdof] .= 0.0
            K_local[cdof, cdof] = m
        end
    end

    # Zero out columns associated to the ghost dofs constrained on a remote process
    # TODO optimize. If we assume that the sparsity pattern is symmetric, then we can constrain
    #      via the column information of the matrix.

    # Step 1: Send out all local ghosts to all other processes...
        partition = K.col_partition.item
        remote_ghost_ldofs = ghost_to_local(partition)
        lto = local_to_owner(partition)
        remote_ghost_parts = lto[remote_ghost_ldofs]
        ltg = local_to_global(partition)
        remote_ghost_gdofs = ltg[remote_ghost_ldofs]

    my_rank = MPI.Comm_rank(comm)+1
    buffer_sizes_send = zeros(Cint, MPI.Comm_size(comm))
    buffer_sizes_recv = Vector{Cint}(undef, MPI.Comm_size(comm))
        for part ∈ remote_ghost_parts
            buffer_sizes_send[part] += 1
        end
    MPI.Alltoall!(UBuffer(buffer_sizes_send, 1), UBuffer(buffer_sizes_recv, 1), comm)
    Ferrite.@debug println("Got $buffer_sizes_recv (R$my_rank)")

    remote_ghosts_recv = Vector{Int}(undef, sum(buffer_sizes_recv))
    MPI.Alltoallv!(VBuffer(remote_ghost_gdofs, buffer_sizes_send), VBuffer(remote_ghosts_recv, buffer_sizes_recv), comm)
    Ferrite.@debug println("Got $remote_ghosts_recv (R$my_rank)")

    # Step 2: Union with all locally constrained dofs
    Ferrite.@debug println("$my_rank : Step 2....")
    remote_ghosts_constrained_send = copy(remote_ghosts_recv)
    prescribed_gdofs = ltg[pdofs]
    for (i, remote_ghost_dof) ∈ enumerate(remote_ghosts_recv)
        remote_ghosts_constrained_send[i] = remote_ghost_dof ∈ prescribed_gdofs
    end

    # Step 3: Send trash back
    Ferrite.@debug println("$my_rank : Step 3....")
    remote_ghosts_constrained_recv = Vector{Int}(undef, sum(buffer_sizes_send))
    MPI.Alltoallv!(VBuffer(remote_ghosts_constrained_send, buffer_sizes_recv), VBuffer(remote_ghosts_constrained_recv, buffer_sizes_send), comm)

    Ferrite.@debug println("$my_rank : remote constraints on $(remote_ghost_ldofs[remote_ghosts_constrained_recv .== 1])")

    # Step 4: Constrain remaining columns
    map(local_values(K)) do K_local
        for cdof ∈ remote_ghost_ldofs[remote_ghosts_constrained_recv .== 1]
            K_local[:, cdof] .= 0.0
        end
    end
end
