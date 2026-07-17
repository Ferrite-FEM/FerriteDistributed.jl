# Hypre to Ferrite vector
function FerriteDistributed.extract_local_part!(u::Vector{T}, uh::HYPREVector, dh::FerriteDistributed.NODDofHandler) where {T}
    # Copy the owned dofs from HYPRE to Julia. They are ordered by ascending global dof
    # number and hence also by ascending local dof number.
    uj = Vector{Float64}(undef, num_local_true_dofs(dh))
    copy!(uj, uh)

    dgrid = getglobalgrid(dh)
    my_rank = global_rank(dgrid)

    next_dof = 1
    for (ldof, rank) in enumerate(dh.ldof_to_rank)
        if rank == my_rank
            u[ldof] = uj[next_dof]
            next_dof += 1
        end
    end

    # Request the values of the non-owned dofs from their owners, which are always
    # interface neighbors.
    ic = InterfaceCommunicator(dgrid)
    request_send = empty_send_buffers(Int, ic)
    ghost_ldofs = empty_send_buffers(Int, ic)
    for (ldof, rank) in enumerate(dh.ldof_to_rank)
        rank == my_rank && continue
        slot = ic.destination_index[rank]
        push!(request_send[slot], dh.ldof_to_gdof[ldof])
        push!(ghost_ldofs[slot], ldof)
    end
    request_recv = exchange(ic, request_send)

    # Answer the requests. The neighborhood is symmetric, so the replies can be routed
    # back through the same channel.
    gdof_offset = first(local_dof_range(dh)) - 1
    reply_send = empty_send_buffers(Float64, ic)
    for (si, gdofs) in enumerate(request_recv)
        reply_send[ic.destination_index[ic.sources[si]]] = Float64[uj[gdof - gdof_offset] for gdof in gdofs]
    end
    reply_recv = exchange(ic, reply_send)

    for (si, values) in enumerate(reply_recv)
        ldofs = ghost_ldofs[ic.destination_index[ic.sources[si]]]
        @assert length(values) == length(ldofs)
        for (ldof, value) in zip(ldofs, values)
            u[ldof] = value
        end
    end

    return u
end
