"""
    NODDofHandler(grid::AbstractNODGrid)

Construct a `NODDofHandler` based on `grid`.

Distributed version of [`DofHandler`](@docs). 

Supports:
- `Grid`s with a single concrete cell type.
- One or several fields on the whole domaine.

!!! todo
    Update to new dof management interface
"""
struct NODDofHandler{dim,T,G<:AbstractNODGrid{dim}} <: Ferrite.AbstractDofHandler
    field_names::Vector{Symbol}
    field_dims::Vector{Int}
    # TODO: field_interpolations can probably be better typed: We should at least require
    #       all the interpolations to have the same dimension and reference shape
    field_interpolations::Vector{Interpolation}
    bc_values::Vector{Ferrite.BCValues{T}} # TODO: BcValues is created/handeld by the constrainthandler, so this can be removed
    cell_dofs::Vector{Int}
    cell_dofs_offset::Vector{Int}
    closed::Ferrite.ScalarWrapper{Bool}
    grid::G
    ndofs::Ferrite.ScalarWrapper{Int}

    ldof_to_gdof::Vector{Int}
    ldof_to_rank::Vector{Int32}
end
# NOTE - REDUNDANT
Ferrite.getfieldnames(dh::NODDofHandler) = dh.field_names
# NOTE - REDUNDANT
Ferrite.nnodes_per_cell(dh::NODDofHandler, cell::Int=1) = Ferrite.nnodes_per_cell(get_grid(dh), cell)
# NOTE - REDUNDANT
function Ferrite.add!(dh::NODDofHandler, name::Symbol, dim::Int, ip::Interpolation=default_interpolation(getcelltype(get_grid(dh))))
    @assert !Ferrite.isclosed(dh)
    @assert !in(name, dh.field_names)
    push!(dh.field_names, name)
    push!(dh.field_dims, dim)
    push!(dh.field_interpolations, ip)
    return dh
end
# NOTE - REDUNDANT
function Ferrite.add!(dh::NODDofHandler, name::Symbol, ip::Interpolation)
    return Ferrite.add!(dh, name, 1, ip)
end
# NOTE - REDUNDANT
function Ferrite.add!(dh::NODDofHandler, name::Symbol, ip::VectorizedInterpolation{vdim}) where {vdim}
    return Ferrite.add!(dh, name, vdim, ip)
end
# NOTE - REDUNDANT
Ferrite.ndofs_per_cell(dh::NODDofHandler, cell::Int=1) = dh.cell_dofs_offset[cell+1] - dh.cell_dofs_offset[cell]
# NOTE - REDUNDANT
function find_field(dh::NODDofHandler, field_name::Symbol)
    j = findfirst(i->i == field_name, dh.field_names)
    j === nothing && error("could not find field :$field_name in DofHandler (existing fields: $(getfieldnames(dh)))")
    return j
end
# NOTE - REDUNDANT
function celldofs!(global_dofs::Vector{Int}, dh::NODDofHandler, i::Int)
    @assert Ferrite.isclosed(dh)
    @assert length(global_dofs) == Ferrite.ndofs_per_cell(dh, i)
    unsafe_copyto!(global_dofs, 1, dh.cell_dofs, dh.cell_dofs_offset[i], length(global_dofs))
    return global_dofs
end
# NOTE - REDUNDANT
function field_offset(dh::NODDofHandler, field_idx::Int)
    offset = 0
    for i in 1:field_idx-1
        offset += getnbasefunctions(Ferrite.getfieldinterpolation(dh,i))::Int
    end
    return offset
end
# NOTE - REDUNDANT
function field_offset(dh::NODDofHandler, field_name::Symbol)
    field_idx = findfirst(i->i == field_name, getfieldnames(dh))
    field_idx === nothing && error("did not find field $field_name")
    return field_offset(dh,field_idx)
end
# NOTE - REDUNDANT
function Ferrite.celldofs!(global_dofs::Vector{Int}, dh::NODDofHandler, i::Int)
    @assert Ferrite.isclosed(dh)
    @assert length(global_dofs) == Ferrite.ndofs_per_cell(dh, i)
    unsafe_copyto!(global_dofs, 1, dh.cell_dofs, dh.cell_dofs_offset[i], length(global_dofs))
    return global_dofs
end
# NOTE - REDUNDANT
function getfielddim(dh::NODDofHandler, field_name::Symbol) 
    field_idx = findfirst(i->i == field_name, getfieldnames(dh))
    field_idx === nothing && error("did not find field $field_name")
    return getfielddim(dh, field_idx)
end
# NOTE - REDUNDANT
function Ferrite.getfieldinterpolation(dh::NODDofHandler, field_idx::Int)
    ip = dh.field_interpolations[field_idx]
    return ip
end
# NOTE - REDUNDANT
Ferrite.getfielddim(dh::NODDofHandler, field_idx::Int) = dh.field_dims[field_idx]
# NOTE - REDUNDANT
function Ferrite.dof_range(dh::NODDofHandler, field_idx::Int)
    offset = field_offset(dh, field_idx)
    n_field_dofs = getnbasefunctions(Ferrite.getfieldinterpolation(dh, field_idx))::Int
    return (offset+1):(offset+n_field_dofs)
end
# NOTE - REDUNDANT
num_fields(dh::NODDofHandler) = length(dh.field_names)

"""
Compute the global dof range of the dofs owned by the calling process. It is guaranteed to be continuous.
"""
function local_dof_range(dh::NODDofHandler)
    my_rank = global_rank(getglobalgrid(dh))
    ltdofs = dh.ldof_to_gdof[dh.ldof_to_rank .== my_rank]
    return minimum(ltdofs):maximum(ltdofs)
end

"""
Construct the correct distributed dof handler from a given distributed grid.
"""
function Ferrite.DofHandler(grid::AbstractNODGrid{dim}) where {dim}
    isconcretetype(getcelltype(grid)) || error("Grid includes different celltypes. DistributedMixedDofHandler not implemented yet.")
    NODDofHandler(Symbol[], Int[], Interpolation[], Ferrite.BCValues{Float64}[], Int[], Int[], Ferrite.ScalarWrapper(false), grid, Ferrite.ScalarWrapper(-1), Int[], Int32[])
end

# NOTE - REDUNDANT
function Base.show(io::IO, ::MIME"text/plain", dh::NODDofHandler)
    println(io, "NODDofHandler")
    println(io, "  Fields:")
    for i in 1:num_fields(dh)
        println(io, "    ", repr(dh.field_names[i]), ", interpolation: ", dh.field_interpolations[i],", dim: ", dh.field_dims[i])
    end
    if !Ferrite.isclosed(dh)
        print(io, "  Not closed!")
    else
        println(io, "  Dofs per cell: ", Ferrite.ndofs_per_cell(dh))
        print(io, "  Total local dofs: ", ndofs(dh))
    end
end

# NOTE - REDUNDANT
Ferrite.getdim(dh::NODDofHandler{dim}) where {dim} = dim 

getlocalgrid(dh::NODDofHandler) = getlocalgrid(dh.grid)
getglobalgrid(dh::NODDofHandler) = dh.grid

# Compat layer against serial code
Ferrite.get_grid(dh::NODDofHandler) = getlocalgrid(dh)

# TODO problem here is that the reorder has to be synchronized. We also cannot arbitrary reorder dofs, 
# because some distributed matrix data structures have strict requirements on the orderings.
Ferrite.renumber!(dh::NODDofHandler, perm::AbstractVector{<:Integer}) = error("Not implemented.")

"""
TODO fix for shells
"""
function compute_dof_ownership(dh::NODDofHandler)
    dgrid = getglobalgrid(dh)
    my_rank = global_rank(dgrid)

    dof_owner = Vector{Int}(undef,ndofs(dh))
    fill!(dof_owner, my_rank)

    for sv ∈ get_shared_vertices(dgrid)
        lvi = sv.local_idx
        for field_idx in 1:num_fields(dh)
            if has_vertex_dofs(dh, field_idx, lvi)
                local_dofs = vertex_dofs(dh, field_idx, lvi)
                dof_owner[local_dofs] .= compute_owner(dgrid, sv)
            end
        end
    end

    for sf ∈ get_shared_faces(dgrid)
        lfi = sf.local_idx
        for field_idx in 1:num_fields(dh)
            if has_face_dofs(dh, field_idx, lfi)
                local_dofs = face_dofs(dh, field_idx, lfi)
                dof_owner[local_dofs] .= compute_owner(dgrid, sf)
            end
        end
    end

    for se ∈ get_shared_edges(dgrid)
        lei = se.local_idx
        for field_idx in 1:num_fields(dh)
            if has_edge_dofs(dh, field_idx, lei)
                local_dofs = edge_dofs(dh, field_idx, lei)
                dof_owner[local_dofs] .= compute_owner(dgrid, se)
            end
        end
    end

    return dof_owner
end

"""
Compute the number of dofs owned by the current process.
"""
num_local_true_dofs(dh::NODDofHandler) = sum(dh.ldof_to_rank .== global_rank(getglobalgrid(dh)))

"""
Compute the number of dofs visible to the current process.
"""
num_local_dofs(dh::NODDofHandler) = length(dh.ldof_to_gdof)

"""
Compute the number of dofs in the global system.
"""
num_global_dofs(dh::NODDofHandler) = MPI.Allreduce(num_local_true_dofs(dh), MPI.SUM, global_comm(getglobalgrid(dh)))

"""
Renumber the dofs in local ordering to their corresponding global numbering.

TODO: Refactor for MixedDofHandler integration
"""
function local_to_global_numbering(dh::NODDofHandler{dim}) where {dim}
    dgrid = getglobalgrid(dh)

    # MPI rank starting with 1 to match Julia's index convention
    my_rank = global_rank(dgrid)

    local_to_global = Vector{Int}(undef,ndofs(dh))
    fill!(local_to_global,0) # 0 is the invalid index!
    # Start by numbering local dofs only from 1:#local_dofs

    # Lookup for synchronization in the form (Remote Rank,Shared Entity)
    # @TODO replace dict with vector and tie to MPI neighborhood graph of the mesh
    vertices_send = Dict{Int,Vector{VertexIndex}}()
    n_vertices_recv = Dict{Int,Int}()
    faces_send = Dict{Int,Vector{FaceIndex}}()
    n_faces_recv = Dict{Int,Int}()
    edges_send = Dict{Int,Vector{EdgeIndex}}()
    edges_recv = Dict{Int,Vector{EdgeIndex}}()

    # We start by assigning a local dof to all owned entities.
    # An entity is owned if:
    # 1. *All* topological neighbors are on the local process
    # 2. If the rank of the local process it lower than the rank of *all* topological neighbors
    # A topological neighbor in this context is hereby defined per entity:
    # * vertex: All elements whose vertex is the vertex in question
    # * cell: Just the cell itself
    # * All other entities: All cells for which one of the corresponding entities interior intersects 
    #                       with the interior of the entity in question.
    next_local_idx = 1
    for (ci, cell) in enumerate(getcells(get_grid(dh)))
        Ferrite.@debug println("cell #$ci (R$my_rank)")
        for field_idx in 1:num_fields(dh)
            Ferrite.@debug println("  field: $(dh.field_names[field_idx]) (R$my_rank)")
            interpolation_info = Ferrite.InterpolationInfo(Ferrite.getfieldinterpolation(dh, field_idx))
            if interpolation_info.nvertexdofs[1] > 0
                for (vi,vertex) in enumerate(Ferrite.vertices(cell))
                    Ferrite.@debug println("    vertex#$vertex (R$my_rank)")
                    lvi = VertexIndex(ci,vi)
                    # Dof is owned if it is local or if my rank is the smallest in the neighborhood
                    if !is_shared_vertex(dgrid, lvi) || (compute_owner(dgrid, get_shared_vertex(dgrid, lvi)) == my_rank)
                        # Update dof assignment
                        dof_local_indices = vertex_dofs(dh, field_idx, lvi)
                        if local_to_global[dof_local_indices[1]] == 0
                            for d in 1:getfielddim(dh, field_idx)
                                Ferrite.@debug println("      mapping vertex dof#$dof_local_indices[d] to $next_local_idx (R$my_rank)")
                                local_to_global[dof_local_indices[d]] = next_local_idx
                                next_local_idx += 1
                            end
                        else
                            for d in 1:getfielddim(dh, field_idx)
                                Ferrite.@debug println("      vertex dof#$(dof_local_indices[d]) already mapped to $(local_to_global[dof_local_indices[d]]) (R$my_rank)")
                            end
                        end
                    end

                    # Update shared vertex lookup table
                    if is_shared_vertex(dgrid, lvi)
                        master_rank = my_rank
                        remote_vertex_dict = remote_entities(get_shared_vertex(dgrid, lvi))
                        for master_rank_new ∈ keys(remote_vertex_dict)
                            master_rank = min(master_rank, master_rank_new)
                        end
                        for (remote_rank, svs) ∈ remote_vertex_dict
                            if master_rank == my_rank # I own the dof - we have to send information
                                if !haskey(vertices_send,remote_rank)
                                    vertices_send[remote_rank] = Vector{VertexIndex}()
                                end
                                Ferrite.@debug println("      prepare sending vertex #$(lvi) to $remote_rank (R$my_rank)")
                                for i ∈ svs
                                    push!(vertices_send[remote_rank],lvi)
                                end
                            elseif master_rank == remote_rank  # dof is owned by remote - we have to receive information
                                if !haskey(n_vertices_recv,remote_rank)
                                    n_vertices_recv[remote_rank] = length(svs)
                                else
                                    n_vertices_recv[remote_rank] += length(svs)
                                end
                                Ferrite.@debug println("      prepare receiving vertex #$(lvi) from $remote_rank (R$my_rank)")
                            end
                        end
                    end
                end
            end

            if dim > 2 # edges only in 3D
                if interpolation_info.nedgedofs[1] > 0
                    for (ei,edge) in enumerate(Ferrite.edges(cell))
                        Ferrite.@debug println("    edge#$edge (R$my_rank)")
                        lei = EdgeIndex(ci,ei)
                        # Dof is owned if it is local or if my rank is the smallest in the neighborhood
                        if !is_shared_edge(dgrid, lei) || (compute_owner(dgrid, get_shared_edge(dgrid, lei)) == my_rank)
                            # Update dof assignment
                            dof_local_indices = edge_dofs(dh, field_idx, lei)
                            if local_to_global[dof_local_indices[1]] == 0
                                for d in 1:getfielddim(dh, field_idx)
                                    Ferrite.@debug println("      mapping edge dof#$(dof_local_indices[d]) to $next_local_idx (R$my_rank)")
                                    local_to_global[dof_local_indices[d]] = next_local_idx
                                    next_local_idx += 1
                                end
                            else
                                for d in 1:getfielddim(dh, field_idx)
                                    Ferrite.@debug println("      edge dof#$(dof_local_indices[d]) already mapped to $(local_to_global[dof_local_indices[d]]) (R$my_rank)")
                                end
                            end
                        end

                        # Update shared edge lookup table
                        if is_shared_edge(dgrid, lei)
                            master_rank = my_rank
                            remote_edge_dict = remote_entities(get_shared_edge(dgrid, lei))
                            for master_rank_new ∈ keys(remote_edge_dict)
                                master_rank = min(master_rank, master_rank_new)
                            end
                            for (remote_rank, svs) ∈ remote_edge_dict
                                if master_rank == my_rank # I own the dof - we have to send information
                                    if !haskey(edges_send,remote_rank)
                                        edges_send[remote_rank] = EdgeIndex[]
                                    end
                                    Ferrite.@debug println("      prepare sending edge #$(lei) to $remote_rank (R$my_rank)")
                                    for i ∈ svs
                                        push!(edges_send[remote_rank], lei)
                                    end
                                elseif master_rank == remote_rank  # dof is owned by remote - we have to receive information
                                    if !haskey(edges_recv,remote_rank)
                                        edges_recv[remote_rank] = EdgeIndex[]
                                    end
                                    push!(edges_recv[remote_rank], lei)
                                    Ferrite.@debug println("      prepare receiving edge #$(lei) from $remote_rank (R$my_rank)")
                                end
                            end
                        end
                    end
                end
            end

            if interpolation_info.nfacedofs[1] > 0 && (interpolation_info.reference_dim == dim)
                for (fi,face) in enumerate(Ferrite.faces(cell))
                    Ferrite.@debug println("    face#$face (R$my_rank)")
                    lfi = FaceIndex(ci,fi)
                    # Dof is owned if it is local or if my rank is the smallest in the neighborhood
                    if !is_shared_face(dgrid, lfi) || (compute_owner(dgrid, get_shared_face(dgrid, lfi)) == my_rank)
                        # Update dof assignment
                        dof_local_indices = face_dofs(dh, field_idx, lfi)
                        if local_to_global[dof_local_indices[1]] == 0
                            for d in 1:getfielddim(dh, field_idx)
                                Ferrite.@debug println("      mapping face dof#$(dof_local_indices[d]) to $next_local_idx (R$my_rank)")
                                local_to_global[dof_local_indices[d]] = next_local_idx
                                next_local_idx += 1
                            end
                        else
                            for d in 1:getfielddim(dh, field_idx)
                                Ferrite.@debug println("      face dof#$(dof_local_indices[d]) already mapped to $(local_to_global[dof_local_indices[d]]) (R$my_rank)")
                            end
                        end
                    end

                    # Update shared face lookup table
                    if is_shared_face(dgrid, lfi)
                        master_rank = my_rank
                        remote_face_dict = remote_entities(get_shared_face(dgrid, lfi))
                        for master_rank_new ∈ keys(remote_face_dict)
                            master_rank = min(master_rank, master_rank_new)
                        end
                        for (remote_rank, svs) ∈ remote_face_dict
                            if master_rank == my_rank # I own the dof - we have to send information
                                if !haskey(faces_send,remote_rank)
                                    faces_send[remote_rank] = FaceIndex[]
                                end
                                Ferrite.@debug println("      prepare sending face #$(lfi) to $remote_rank (R$my_rank)")
                                for i ∈ svs
                                    push!(faces_send[remote_rank],lfi)
                                end
                            elseif master_rank == remote_rank  # dof is owned by remote - we have to receive information
                                if !haskey(n_faces_recv,remote_rank)
                                    n_faces_recv[remote_rank] = length(svs)
                                else
                                    n_faces_recv[remote_rank] += length(svs)
                                end
                                Ferrite.@debug println("      prepare receiving face #$(lfi) from $remote_rank (R$my_rank)")
                            end
                        end
                    end
                end # face loop
            end

            if interpolation_info.ncelldofs > 0 # always distribute new dofs for cell
                Ferrite.@debug println("    cell#$ci")
                if interpolation_info.ncelldofs > 0
                    # Update dof assignment
                    dof_local_indices = cell_dofs(dh, field_idx, ci)
                    if local_to_global[dof_local_indices[1]] == 0
                        for d in 1:getfielddim(dh, field_idx)
                            Ferrite.@debug println("      mapping cell dof#$(dof_local_indices[d]) to $next_local_idx (R$my_rank)")
                            local_to_global[dof_local_indices[d]] = next_local_idx
                            next_local_idx += 1
                        end
                    else
                        for d in 1:getfielddim(dh, field_idx)
                            # Should never happen...
                            Ferrite.@debug println("      WARNING! cell dof#$(dof_local_indices[d]) already mapped to $(local_to_global[dof_local_indices[d]]) (R$my_rank)")
                        end
                    end
                end # cell loop
            end
        end # field loop
    end

    #
    num_true_local_dofs = next_local_idx-1
    Ferrite.@debug println("#true local dofs $num_true_local_dofs (R$my_rank)")

    # @TODO optimize the following synchronization with MPI line graph topology 
    # and allgather
    # Set true local indices
    local_offset = 0
    if my_rank > 1
        local_offset = MPI.Recv(Int, global_comm(dgrid); source=my_rank-1-1)
    end
    if my_rank < MPI.Comm_size(global_comm(dgrid))
        MPI.Send(local_offset+num_true_local_dofs, global_comm(dgrid); dest=my_rank+1-1)
    end
    Ferrite.@debug println("#shifted local dof range $(local_offset+1):$(local_offset+num_true_local_dofs) (R$my_rank)")

    # Shift assigned local dofs (dofs with value >0) into the global range
    # At this point in the algorithm the dofs with value 0 are the dofs owned of neighboring processes
    for i ∈ 1:length(local_to_global)
        if local_to_global[i] != 0
            local_to_global[i] += local_offset
        end
    end

    # Sync non-owned dofs with neighboring processes.
    # TODO: Use MPI graph primitives to simplify this code
    # TODO: Simplify with dimension-agnostic code...
    for sending_rank ∈ 1:MPI.Comm_size(global_comm(dgrid))
        if my_rank == sending_rank
            for remote_rank ∈ 1:MPI.Comm_size(global_comm(dgrid))
                if haskey(vertices_send, remote_rank)
                    n_vertices = length(vertices_send[remote_rank])
                    Ferrite.@debug println("Sending $n_vertices vertices to rank $remote_rank (R$my_rank)")
                    remote_cells = Array{Int64}(undef,n_vertices)
                    remote_cell_vis = Array{Int64}(undef,n_vertices)
                    next_buffer_idx = 1
                    for lvi ∈ vertices_send[remote_rank]
                        sv = dgrid.shared_vertices[lvi]
                        @assert haskey(sv.remote_vertices, remote_rank)
                        for (cvi, llvi) ∈ sv.remote_vertices[remote_rank][1:1] # Just don't ask :)
                            remote_cells[next_buffer_idx] = cvi
                            remote_cell_vis[next_buffer_idx] = llvi
                            next_buffer_idx += 1
                        end
                    end
                    MPI.Send(remote_cells, global_comm(dgrid); dest=remote_rank-1)
                    MPI.Send(remote_cell_vis, global_comm(dgrid); dest=remote_rank-1)
                    for field_idx ∈ 1:num_fields(dh)
                        next_buffer_idx = 1
                        ip = Ferrite.getfieldinterpolation(dh, field_idx)
                        if nvertexdofs(ip) == 0
                            Ferrite.@debug println("Skipping send vertex on field $(dh.field_names[field_idx]) (R$my_rank)")
                            continue
                        end
                        # fill correspondence array
                        corresponding_global_dofs = Array{Int64}(undef,n_vertices)
                        for vertex ∈ vertices_send[remote_rank]
                            if has_vertex_dofs(dh, field_idx, vertex)
                                # We just put the first dof into the array to reduce communication
                                vdofs = vertex_dofs(dh, field_idx, vertex)
                                corresponding_global_dofs[next_buffer_idx] = local_to_global[vdofs[1]]
                            end
                            next_buffer_idx += 1
                        end
                        MPI.Send(corresponding_global_dofs, global_comm(dgrid); dest=remote_rank-1)
                    end
                end

                if haskey(faces_send, remote_rank)
                    n_faces = length(faces_send[remote_rank])
                    Ferrite.@debug println("Sending $n_faces faces to rank $remote_rank (R$my_rank)")
                    remote_cells = Array{Int64}(undef,n_faces)
                    remote_cell_vis = Array{Int64}(undef,n_faces)
                    next_buffer_idx = 1
                    for lvi ∈ faces_send[remote_rank]
                        sv = dgrid.shared_faces[lvi]
                        @assert haskey(sv.remote_faces, remote_rank)
                        for (cvi, llvi) ∈ sv.remote_faces[remote_rank][1:1] # Just don't ask :)
                            remote_cells[next_buffer_idx] = cvi
                            remote_cell_vis[next_buffer_idx] = llvi 
                            next_buffer_idx += 1
                        end
                    end
                    MPI.Send(remote_cells, global_comm(dgrid); dest=remote_rank-1)
                    MPI.Send(remote_cell_vis, global_comm(dgrid); dest=remote_rank-1)
                    for field_idx ∈ 1:num_fields(dh)
                        next_buffer_idx = 1
                        ip = Ferrite.getfieldinterpolation(dh, field_idx)
                        if nfacedofs(ip) == 0
                            Ferrite.@debug println("Skipping send faces on field $(dh.field_names[field_idx]) (R$my_rank)")
                            continue
                        end
                        # fill correspondence array
                        corresponding_global_dofs = Array{Int64}(undef,n_faces)
                        for face ∈ faces_send[remote_rank]
                            if has_face_dofs(dh, field_idx, face)
                                # We just put the first dof into the array to reduce communication
                                fdofs = face_dofs(dh, field_idx, face)
                                corresponding_global_dofs[next_buffer_idx] = local_to_global[fdofs[1]]
                            end
                            next_buffer_idx += 1
                        end
                        MPI.Send(corresponding_global_dofs, global_comm(dgrid); dest=remote_rank-1)
                    end
                end

                if haskey(edges_send, remote_rank)
                    # Well .... that some hotfix straight outta hell.
                    edges_send_unique_set = Set{Tuple{Int,Int}}()
                    edges_send_unique = Set{EdgeIndex}()
                    for lei ∈ edges_send[remote_rank]
                        edge = Ferrite.toglobal(dgrid, lei)
                        if edge ∉ edges_send_unique_set
                            push!(edges_send_unique_set, edge)
                            push!(edges_send_unique, lei)
                        end
                    end
                    n_edges = length(edges_send_unique)
                    Ferrite.@debug println("Sending $n_edges edges to rank $remote_rank (R$my_rank)")
                    remote_cells = Array{Int64}(undef,n_edges)
                    remote_cell_vis = Array{Int64}(undef,n_edges)
                    next_buffer_idx = 1
                    for lvi ∈ edges_send_unique
                        sv = dgrid.shared_edges[lvi]
                        @assert haskey(sv.remote_edges, remote_rank)
                        for (cvi, llvi) ∈ sv.remote_edges[remote_rank][1:1] # Just don't ask :)
                            remote_cells[next_buffer_idx] = cvi
                            remote_cell_vis[next_buffer_idx] = llvi 
                            next_buffer_idx += 1
                        end
                    end
                    MPI.Send(remote_cells, global_comm(dgrid); dest=remote_rank-1)
                    MPI.Send(remote_cell_vis, global_comm(dgrid); dest=remote_rank-1)
                    for field_idx ∈ 1:num_fields(dh)
                        next_buffer_idx = 1
                        ip = Ferrite.getfieldinterpolation(dh, field_idx)
                        if nedgedofs(ip) == 0
                            Ferrite.@debug println("Skipping send edges on field $(dh.field_names[field_idx]) (R$my_rank)")
                            continue
                        end
                        # fill correspondence array
                        corresponding_global_dofs = Array{Int64}(undef,n_edges)
                        for edge ∈ edges_send_unique
                            if has_edge_dofs(dh, field_idx, edge)
                                # We just put the first dof into the array to reduce communication
                                edofs = edge_dofs(dh, field_idx, edge)
                                corresponding_global_dofs[next_buffer_idx] = local_to_global[edofs[1]]
                            end
                            next_buffer_idx += 1
                        end
                        MPI.Send(corresponding_global_dofs, global_comm(dgrid); dest=remote_rank-1)
                    end
                end
            end
        else
            if haskey(n_vertices_recv, sending_rank)
                n_vertices = n_vertices_recv[sending_rank]
                Ferrite.@debug println("Receiving $n_vertices vertices from rank $sending_rank (R$my_rank)")
                local_cells = Array{Int64}(undef,n_vertices)
                local_cell_vis = Array{Int64}(undef,n_vertices)
                MPI.Recv!(local_cells, global_comm(dgrid); source=sending_rank-1)
                MPI.Recv!(local_cell_vis, global_comm(dgrid); source=sending_rank-1)
                for field_idx in 1:num_fields(dh)
                    ip = Ferrite.getfieldinterpolation(dh, field_idx)
                    if nvertexdofs(ip) == 0
                        Ferrite.@debug println("  Skipping recv of vertices on field $(dh.field_names[field_idx]) (R$my_rank)")
                        continue
                    end
                    corresponding_global_dofs = Array{Int64}(undef,n_vertices)
                    MPI.Recv!(corresponding_global_dofs, global_comm(dgrid); source=sending_rank-1)
                    for (cdi,vertex) ∈ enumerate(VertexIndex.(zip(local_cells,local_cell_vis)))
                        if has_vertex_dofs(dh, field_idx, vertex)
                            vdofs = vertex_dofs(dh, field_idx, vertex)
                            for d in 1:getfielddim(dh, field_idx)
                                local_to_global[vdofs[d]] = corresponding_global_dofs[cdi]+d-1
                                Ferrite.@debug println("  Updating field $(dh.field_names[field_idx]) vertex $vertex to $(corresponding_global_dofs[cdi]+d-1) (R$my_rank)")
                            end
                        else
                            Ferrite.@debug println("  Skipping recv on field $(dh.field_names[field_idx]) vertex $vertex (R$my_rank)")
                        end
                    end
                end
            end

            if haskey(n_faces_recv, sending_rank)
                n_faces = n_faces_recv[sending_rank]
                Ferrite.@debug println("Receiving $n_faces faces from rank $sending_rank (R$my_rank)")
                local_cells = Array{Int64}(undef,n_faces)
                local_cell_vis = Array{Int64}(undef,n_faces)
                MPI.Recv!(local_cells, global_comm(dgrid); source=sending_rank-1)
                MPI.Recv!(local_cell_vis, global_comm(dgrid); source=sending_rank-1)
                for field_idx in 1:num_fields(dh)
                    ip = Ferrite.getfieldinterpolation(dh, field_idx)
                    if nfacedofs(ip) == 0
                        Ferrite.@debug println("  Skipping recv of faces on field $(dh.field_names[field_idx]) (R$my_rank)")
                        continue
                    end
                    corresponding_global_dofs = Array{Int64}(undef,n_faces)
                    MPI.Recv!(corresponding_global_dofs, global_comm(dgrid); source=sending_rank-1)
                    for (cdi,face) ∈ enumerate(FaceIndex.(zip(local_cells,local_cell_vis)))
                        if has_face_dofs(dh, field_idx, face)
                            fdofs = face_dofs(dh, field_idx, face)
                            for d in 1:getfielddim(dh, field_idx)
                                local_to_global[fdofs[d]] = corresponding_global_dofs[cdi]+d-1
                                Ferrite.@debug println("  Updating field $(dh.field_names[field_idx]) face $face to $(corresponding_global_dofs[cdi]) (R$my_rank)")
                            end
                        else
                            Ferrite.@debug println("  Skipping recv on field $(dh.field_names[field_idx]) face $face (R$my_rank)")
                        end
                    end
                end
            end

            if haskey(edges_recv, sending_rank)
                edges_recv_unique_set = Set{Tuple{Int,Int}}()
                for lei ∈ edges_recv[sending_rank]
                    edge = Ferrite.toglobal(dgrid, lei)
                    push!(edges_recv_unique_set, edge)
                end
                n_edges = length(edges_recv_unique_set)
                Ferrite.@debug println("Receiving $n_edges edges from rank $sending_rank (R$my_rank)")
                local_cells = Array{Int64}(undef,n_edges)
                local_cell_vis = Array{Int64}(undef,n_edges)
                MPI.Recv!(local_cells, global_comm(dgrid); source=sending_rank-1)
                MPI.Recv!(local_cell_vis, global_comm(dgrid); source=sending_rank-1)
                for field_idx in 1:num_fields(dh)
                    ip = Ferrite.getfieldinterpolation(dh, field_idx)
                    if nedgedofs(ip) == 0
                        Ferrite.@debug println("  Skipping recv on field $(dh.field_names[field_idx]) (R$my_rank)")
                        continue
                    end
                    corresponding_global_dofs = Array{Int64}(undef,n_edges)
                    MPI.Recv!(corresponding_global_dofs, global_comm(dgrid); source=sending_rank-1)
                    Ferrite.@debug println("   Received $corresponding_global_dofs edge dofs from $sending_rank (R$my_rank)")
                    for (cdi,edge) ∈ enumerate(EdgeIndex.(zip(local_cells,local_cell_vis)))
                        if has_edge_dofs(dh, field_idx, edge)
                            edofs = edge_dofs(dh, field_idx, edge)
                            for d in 1:getfielddim(dh, field_idx)
                                local_to_global[edofs[d]] = corresponding_global_dofs[cdi]+d-1
                                Ferrite.@debug println("  Updating field $(dh.field_names[field_idx]) edge $edge to $(corresponding_global_dofs[cdi]) (R$my_rank)")
                            end
                        else
                            Ferrite.@debug println("  Skipping recv on field $(dh.field_names[field_idx]) edge $edge (R$my_rank)")
                        end
                    end
                end
            end
        end
    end

    # Postcondition: All local dofs need a corresponding global dof!
    Ferrite.@debug println("Local to global mapping: $local_to_global (R$my_rank)")
    @assert findfirst(local_to_global .== 0) === nothing

    return local_to_global
end

function Ferrite.close!(dh::NODDofHandler)
    # We could merge these functions into an optimized one if we want.
    Ferrite.__close!(dh)
    append!(dh.ldof_to_rank, compute_dof_ownership(dh))
    append!(dh.ldof_to_gdof, local_to_global_numbering(dh))
    return dh
end

# NOTE - REDUNDANT
function Ferrite.__close!(dh::NODDofHandler{dim}) where {dim}
    @assert !Ferrite.isclosed(dh)

    # `vertexdict` keeps track of the visited vertices. We store the global vertex
    # number and the first dof we added to that vertex.
    vertexdicts = [Dict{Int,Int}() for _ in 1:num_fields(dh)]

    # `edgedict` keeps track of the visited edges, this will only be used for a 3D problem
    # An edge is determined from two vertices, but we also need to store the direction
    # of the first edge we encounter and add dofs too. When we encounter the same edge
    # the next time we check if the direction is the same, otherwise we reuse the dofs
    # in the reverse order
    edgedicts = [Dict{Tuple{Int,Int},Tuple{Int,Bool}}() for _ in 1:num_fields(dh)]

    # `facedict` keeps track of the visited faces. We only need to store the first dof we
    # added to the face; if we encounter the same face again we *always* reverse the order
    # In 2D a face (i.e. a line) is uniquely determined by 2 vertices, and in 3D a
    # face (i.e. a surface) is uniquely determined by 3 vertices.
    facedicts = [Dict{NTuple{dim,Int},Int}() for _ in 1:num_fields(dh)]

    # celldofs are never shared between different cells so there is no need
    # for a `celldict` to keep track of which cells we have added dofs too.

    # We create the `InterpolationInfo` structs with precomputed information for each
    # interpolation since that allows having the cell loop as the outermost loop,
    # and the interpolation loop inside without using a function barrier
    interpolation_infos = Ferrite.InterpolationInfo[]
    for interpolation in dh.field_interpolations
        # push!(dh.interpolation_info, InterpolationInfo(interpolation))
        push!(interpolation_infos, Ferrite.InterpolationInfo(interpolation))
    end

    # not implemented yet: more than one facedof per face in 3D
    # dim == 3 && @assert(!any(x->x.nfacedofs > 1, interpolation_infos))

    nextdof = 1 # next free dof to distribute
    push!(dh.cell_dofs_offset, 1) # dofs for the first cell start at 1

    # loop over all the cells, and distribute dofs for all the fields
    for (ci, cell) in enumerate(getcells(get_grid(dh)))
        @debug println("cell #$ci")
        for fi in 1:num_fields(dh)
            interpolation_info = interpolation_infos[fi]
            @debug println("  field: $(dh.field_names[fi])")
            if interpolation_info.nvertexdofs[1] > 0
                for vertex in Ferrite.vertices(cell)
                    @debug println("    vertex#$vertex")
                    token = Base.ht_keyindex2!(vertexdicts[fi], vertex)
                    if token > 0 # haskey(vertexdicts[fi], vertex) # reuse dofs
                        reuse_dof = vertexdicts[fi].vals[token] # vertexdicts[fi][vertex]
                        for d in 1:dh.field_dims[fi]
                            @debug println("      reusing dof #$(reuse_dof + (d-1))")
                            push!(dh.cell_dofs, reuse_dof + (d-1))
                        end
                    else # token <= 0, distribute new dofs
                        for vertexdof in 1:interpolation_info.nvertexdofs[1]
                            Base._setindex!(vertexdicts[fi], nextdof, vertex, -token) # vertexdicts[fi][vertex] = nextdof
                            for d in 1:dh.field_dims[fi]
                                @debug println("      adding dof#$nextdof")
                                push!(dh.cell_dofs, nextdof)
                                nextdof += 1
                            end
                        end
                    end
                end # vertex loop
            end
            if dim == 3 # edges only in 3D
                if interpolation_info.nedgedofs[1] > 0
                    for edge in Ferrite.edges(cell)
                        sedge, dir = Ferrite.sortedge(edge)
                        @debug println("    edge#$sedge dir: $(dir)")
                        token = Base.ht_keyindex2!(edgedicts[fi], sedge)
                        if token > 0 # haskey(edgedicts[fi], sedge), reuse dofs
                            startdof, olddir = edgedicts[fi].vals[token] # edgedicts[fi][sedge] # first dof for this edge (if dir == true)
                            for edgedof in (dir.regular == olddir ? (1:interpolation_info.nedgedofs[1]) : (interpolation_info.nedgedofs[1]:-1:1))
                                for d in 1:dh.field_dims[fi]
                                    reuse_dof = startdof + (d-1) + (edgedof-1)*dh.field_dims[fi]
                                    @debug println("      reusing dof#$(reuse_dof)")
                                    push!(dh.cell_dofs, reuse_dof)
                                end
                            end
                        else # token <= 0, distribute new dofs
                            Base._setindex!(edgedicts[fi], (nextdof, dir.regular), sedge, -token) # edgedicts[fi][sedge] = (nextdof, dir),  store only the first dof for the edge
                            for edgedof in 1:interpolation_info.nedgedofs[1]
                                for d in 1:dh.field_dims[fi]
                                    @debug println("      adding dof#$nextdof")
                                    push!(dh.cell_dofs, nextdof)
                                    nextdof += 1
                                end
                            end
                        end
                    end # edge loop
                end
            end
            if interpolation_info.nfacedofs[1] > 0 && (interpolation_info.reference_dim == dim)
                for face in Ferrite.faces(cell)
                    sface, dir = Ferrite.sortface(face) # TODO: faces(cell) may as well just return the sorted list
                    @debug println("    face#$sface")
                    token = Base.ht_keyindex2!(facedicts[fi], sface)
                    if token > 0 # haskey(facedicts[fi], sface), reuse dofs
                        startdof = facedicts[fi].vals[token] # facedicts[fi][sface]
                        for facedof in interpolation_info.nfacedofs[1]:-1:1 # always reverse (YOLO)
                            for d in 1:dh.field_dims[fi]
                                reuse_dof = startdof + (d-1) + (facedof-1)*dh.field_dims[fi]
                                @debug println("      reusing dof#$(reuse_dof)")
                                push!(dh.cell_dofs, reuse_dof)
                            end
                        end
                    else # distribute new dofs
                        Base._setindex!(facedicts[fi], nextdof, sface, -token)# facedicts[fi][sface] = nextdof,  store the first dof for this face
                        for facedof in 1:interpolation_info.nfacedofs[1]
                            for d in 1:dh.field_dims[fi]
                                @debug println("      adding dof#$nextdof")
                                push!(dh.cell_dofs, nextdof)
                                nextdof += 1
                            end
                        end
                    end
                end # face loop
            end
            if interpolation_info.ncelldofs > 0 # always distribute new dofs for cell
                @debug println("    cell#$ci")
                for celldof in 1:interpolation_info.ncelldofs
                    for d in 1:dh.field_dims[fi]
                        @debug println("      adding dof#$nextdof")
                        push!(dh.cell_dofs, nextdof)
                        nextdof += 1
                    end
                end # cell loop
            end
        end # field loop
        # push! the first index of the next cell to the offset vector
        push!(dh.cell_dofs_offset, length(dh.cell_dofs)+1)
    end # cell loop
    dh.ndofs[] = maximum(dh.cell_dofs)
    dh.closed[] = true

    return dh, vertexdicts, edgedicts, facedicts
end

# NOTE - REDUNDANT
function Ferrite.reinit!(cc::CellCache{<:Any,<:Ferrite.AbstractGrid,<:NODDofHandler}, i::Int)
    cc.cellid[] = i
    if cc.flags.nodes
        cellnodes!(cc.nodes, cc.grid, i)
    end
    if cc.flags.coords
        cellcoords!(cc.coords, cc.grid, i)
    end
    if cc.dh !== nothing && cc.flags.dofs
        celldofs!(cc.dofs, cc.dh, i)
    end
    return cc
end

# NOTE - REDUNDANT
function Ferrite.CellCache(dh::NODDofHandler{dim}, flags::UpdateFlags=UpdateFlags()) where {dim}
    N = Ferrite.nnodes_per_cell(get_grid(dh))
    nodes = zeros(Int, N)
    coords = zeros(Vec{dim, get_coordinate_eltype(get_grid(dh))}, N)
    n = ndofs_per_cell(dh)
    celldofs = zeros(Int, n)
    return Ferrite.CellCache(flags, get_grid(dh), ScalarWrapper(-1), nodes, coords, dh, celldofs)
end

# NOTE - REDUNDANT
function Ferrite.evaluate_at_grid_nodes(dh::NODDofHandler, u::Vector, fieldname::Symbol)
    return Ferrite._evaluate_at_grid_nodes(dh, u, fieldname)
end

# NOTE - REDUNDANT?
function Ferrite._evaluate_at_grid_nodes(dh::NODDofHandler, u::Vector{T}, fieldname::Symbol, ::Val{vtk}=Val(false)) where {T, vtk}
    # Make sure the field exists
    fieldname ∈ getfieldnames(dh) || error("Field $fieldname not found.")
    # Figure out the return type (scalar or vector)
    field_idx = find_field(dh, fieldname)
    ip = getfieldinterpolation(dh, field_idx)
    RT = ip isa ScalarInterpolation ? T : Vec{n_components(ip),T}
    if vtk
        # VTK output of solution field (or L2 projected scalar data)
        n_c = n_components(ip)
        vtk_dim = n_c == 2 ? 3 : n_c # VTK wants vectors padded to 3D
        data = fill(NaN * zero(T), vtk_dim, getnnodes(dh.grid))
    else
        # Just evalutation at grid nodes
        data = fill(NaN * zero(RT), getnnodes(dh.grid))
    end
    # Loop over the fieldhandlers
    # for fh in dh.fieldhandlers
    fh = FieldHandler([Field(dh.field_names[i], dh.field_interpolations[i]) for i ∈ 1:length(dh.field_names)], Set(1:getncells(get_grid(dh)))) # TODO REMOVE THIS HOTFIX
        # Check if this fh contains this field, otherwise continue to the next
        field_idx = Ferrite.find_field(fh, fieldname)
        # field_idx === nothing && continue
        # Set up CellValues with the local node coords as quadrature points
        CT = getcelltype(dh.grid, first(fh.cellset))
        ip_geo = default_interpolation(CT)
        local_node_coords = reference_coordinates(ip_geo)
        qr = QuadratureRule{getrefshape(ip)}(zeros(length(local_node_coords)), local_node_coords)
        ip = getfieldinterpolation(fh, field_idx)
        if ip isa VectorizedInterpolation
            # TODO: Remove this hack when embedding works...
            cv = CellValues(qr, ip.ip, ip_geo)
        else
            cv = CellValues(qr, ip, ip_geo)
        end
        drange = dof_range(fh, fieldname)
        # Function barrier
        Ferrite._evaluate_at_grid_nodes!(data, dh, fh, u, cv, drange, RT)
    # end
    return data
end

# NOTE - REDUNDANT
function Ferrite._evaluate_at_grid_nodes!(data::Union{Vector,Matrix}, dh::NODDofHandler, fh::FieldHandler,
        u::Vector{T}, cv::CellValues, drange::UnitRange, ::Type{RT}) where {T, RT}
    ue = zeros(T, length(drange))
    # TODO: Remove this hack when embedding works...
    if RT <: Vec && cv isa CellValues{<:ScalarInterpolation}
        uer = reinterpret(RT, ue)
    else
        uer = ue
    end
    for cell in CellIterator(dh, fh.cellset)
        # Note: We are only using the shape functions: no reinit!(cv, cell) necessary
        @assert getnquadpoints(cv) == length(cell.nodes)
        for (i, I) in pairs(drange)
            ue[i] = u[cell.dofs[I]]
        end
        for (qp, nodeid) in pairs(cell.nodes)
            val = function_value(cv, qp, uer)
            if data isa Matrix # VTK
                data[1:length(val), nodeid] .= val
                data[(length(val)+1):end, nodeid] .= 0 # purge the NaN
            else
                data[nodeid] = val
            end
        end
    end
    return data
end

# TODO REMOVEME - This is legacy code ported to make the example work as is on the distributed assembly PR
function Ferrite.add!(ch::ConstraintHandler{DH}, dbc::Dirichlet) where {DH <: NODDofHandler}
    if length(dbc.faces) == 0
        @warn("adding Dirichlet Boundary Condition to set containing 0 entities")
    end
    celltype = getcelltype(get_grid(ch.dh))
    @assert isconcretetype(celltype)

    # Extract stuff for the field
    field_idx = find_field(ch.dh, dbc.field_name) # throws if name not found
    interpolation = getfieldinterpolation(ch.dh, field_idx)
    if interpolation isa VectorizedInterpolation
        interpolation = interpolation.ip
    end
    field_dim = getfielddim(ch.dh, field_idx)

    if !all(c -> 0 < c <= field_dim, dbc.components)
        error("components $(dbc.components) not within range of field :$(dbc.field_name) ($(field_dim) dimension(s))")
    end

    # Empty components means constrain them all
    isempty(dbc.components) && append!(dbc.components, 1:field_dim)

    if eltype(dbc.faces)==Int #Special case when dbc.faces is a nodeset
        bcvalue = BCValues(interpolation, default_interpolation(celltype), FaceIndex) #Not used by node bcs, but still have to pass it as an argument
    else
        bcvalue = BCValues(interpolation, default_interpolation(celltype), eltype(dbc.faces))
    end
    _add!(ch, dbc, dbc.faces, interpolation, field_dim, field_offset(ch.dh, dbc.field_name), bcvalue)

    return ch
end


function _add!(ch::ConstraintHandler, dbc::Dirichlet, bcfaces::Set{Index}, interpolation::Interpolation, field_dim::Int, offset::Int, bcvalue::BCValues, cellset::Set{Int}=Set{Int}(1:getncells(get_grid(ch.dh)))) where {Index<:BoundaryIndex}
    local_face_dofs, local_face_dofs_offset =
        _local_face_dofs_for_bc(interpolation, field_dim, dbc.components, offset, boundaryfunction(eltype(bcfaces)))
    copy!(dbc.local_face_dofs, local_face_dofs)
    copy!(dbc.local_face_dofs_offset, local_face_dofs_offset)

    # loop over all the faces in the set and add the global dofs to `constrained_dofs`
    constrained_dofs = Int[]
    cc = CellCache(ch.dh, UpdateFlags(; nodes=false, coords=false, dofs=true))
    for (cellidx, faceidx) in bcfaces
        if cellidx ∉ cellset
            delete!(dbc.faces, Index(cellidx, faceidx))
            continue # skip faces that are not part of the cellset
        end
        reinit!(cc, cellidx)
        r = local_face_dofs_offset[faceidx]:(local_face_dofs_offset[faceidx+1]-1)
        append!(constrained_dofs, cc.dofs[local_face_dofs[r]]) # TODO: for-loop over r and simply push! to ch.prescribed_dofs
        @debug println("adding dofs $(cc.dofs[local_face_dofs[r]]) to dbc")
    end

    # save it to the ConstraintHandler
    push!(ch.dbcs, dbc)
    push!(ch.bcvalues, bcvalue)
    for d in constrained_dofs
        Ferrite.add_prescribed_dof!(ch, d, NaN, nothing)
    end
    return ch
end

function _local_face_dofs_for_bc(interpolation, field_dim, components, offset, boundaryfunc::F=faces) where F
    @assert issorted(components)
    local_face_dofs = Int[]
    local_face_dofs_offset = Int[1]
    for (_, face) in enumerate(boundaryfunc(interpolation))
        for fdof in face, d in 1:field_dim
            if d in components
                push!(local_face_dofs, (fdof-1)*field_dim + d + offset)
            end
        end
        push!(local_face_dofs_offset, length(local_face_dofs) + 1)
    end
    return local_face_dofs, local_face_dofs_offset
end
