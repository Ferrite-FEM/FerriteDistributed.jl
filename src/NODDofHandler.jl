# Entity kinds used to tag messages in the dof synchronization protocol.
const VERTEX_KIND = 1
const EDGE_KIND = 2
const FACE_KIND = 3

"""
    EntityDofInfo

Entity dof bookkeeping produced by `Ferrite.__close!` on the local dof handler. The dicts
are indexed by the (global) field index and keyed by local node ids of the entity, mapping
to the first local dof of the (entity, field) dof block. `nentitydofs[kind, field]` is the
size of such a block (including all components).
"""
struct EntityDofInfo
    vertexdicts::Vector{Vector{Int}}
    edgedicts::Vector{Dict{NTuple{2,Int},Int}}
    facedicts::Vector{Dict{NTuple{3,Int},Int}}
    nentitydofs::Matrix{Int} # 3 × nfields
end

"""
    NODDofHandler(grid::AbstractNODGrid)

Distributed dof handler for non-overlapping distributed grids. Wraps a standard
`Ferrite.DofHandler` acting on the local grid and extends it with a global numbering
of the dofs across all processes.

Fields are added with `add!` or via `SubDofHandler`s, exactly as for the serial
`DofHandler`, and the construction is finalized with `close!`. Fields restricted to
subdomains are fully supported; ownership and global numbering of dofs on the process
boundary are negotiated per (entity, field) pair.
"""
mutable struct NODDofHandler{dim,G<:AbstractNODGrid{dim},LDH<:Ferrite.DofHandler{dim}} <: Ferrite.AbstractDofHandler
    const ldh::LDH
    const grid::G
    const ldof_to_gdof::Vector{Int}
    const ldof_to_rank::Vector{Int32}
    entity_dofs::Union{EntityDofInfo,Nothing}
    gdof_offset::Int # number of dofs owned by lower ranks
end

function NODDofHandler(grid::AbstractNODGrid)
    ldh = Ferrite.DofHandler(getlocalgrid(grid))
    return NODDofHandler(ldh, grid, Int[], Int32[], nothing, -1)
end

"""
Construct the correct distributed dof handler from a given distributed grid.
"""
Ferrite.DofHandler(grid::AbstractNODGrid) = NODDofHandler(grid)

# Forward the structural `DofHandler` fields to the local dof handler so that generic
# Ferrite code (e.g. the ConstraintHandler machinery) works on the wrapper.
const _LDH_FORWARDED_FIELDS = (:subdofhandlers, :field_names, :cell_dofs, :cell_dofs_offset, :cell_to_subdofhandler, :closed, :ndofs)

function Base.getproperty(dh::NODDofHandler, s::Symbol)
    s in _LDH_FORWARDED_FIELDS && return getproperty(getfield(dh, :ldh), s)
    return getfield(dh, s)
end
Base.propertynames(dh::NODDofHandler) = (fieldnames(typeof(dh))..., _LDH_FORWARDED_FIELDS...)

"Get the dof handler of the local grid."
getlocaldofhandler(dh::NODDofHandler) = getfield(dh, :ldh)
getlocalgrid(dh::NODDofHandler) = getlocalgrid(getfield(dh, :grid))
getglobalgrid(dh::NODDofHandler) = getfield(dh, :grid)

# Compat layer against serial code
Ferrite.get_grid(dh::NODDofHandler) = getlocalgrid(dh)
Ferrite.getspatialdim(::NODDofHandler{dim}) where {dim} = dim

entity_dof_info(dh::NODDofHandler) = getfield(dh, :entity_dofs)::EntityDofInfo

num_fields(dh::NODDofHandler) = length(Ferrite.getfieldnames(dh))

# Forwards to the local dof handler
function Ferrite.add!(dh::NODDofHandler, name::Symbol, ip::Interpolation)
    Ferrite.add!(getlocaldofhandler(dh), name, ip)
    return dh
end
Ferrite.SubDofHandler(dh::NODDofHandler, cellset::Ferrite.AbstractVecOrSet{Int}) = Ferrite.SubDofHandler(getlocaldofhandler(dh), cellset)
Ferrite.ndofs_per_cell(dh::NODDofHandler) = Ferrite.ndofs_per_cell(getlocaldofhandler(dh))
Ferrite.ndofs_per_cell(dh::NODDofHandler, cell::Int) = Ferrite.ndofs_per_cell(getlocaldofhandler(dh), cell)
# Both methods needed to avoid dispatch ambiguities with the Ferrite counterparts
Ferrite.celldofs!(global_dofs::Vector{Int}, dh::NODDofHandler, i::Int) = Ferrite.celldofs!(global_dofs, getlocaldofhandler(dh), i)
Ferrite.celldofs!(global_dofs::AbstractVector{Int}, dh::NODDofHandler, i::Int) = Ferrite.celldofs!(global_dofs, getlocaldofhandler(dh), i)
Ferrite.celldofs(dh::NODDofHandler, i::Int) = Ferrite.celldofs(getlocaldofhandler(dh), i)
Ferrite.getfieldnames(dh::NODDofHandler) = Ferrite.getfieldnames(getlocaldofhandler(dh))
Ferrite.find_field(dh::NODDofHandler, field_name::Symbol) = Ferrite.find_field(getlocaldofhandler(dh), field_name)
Ferrite.getfieldinterpolation(dh::NODDofHandler, field_idx) = Ferrite.getfieldinterpolation(getlocaldofhandler(dh), field_idx)
Ferrite.dof_range(dh::NODDofHandler, name::Symbol) = Ferrite.dof_range(getlocaldofhandler(dh), name)
Ferrite.n_components(dh::NODDofHandler, field) = Ferrite.n_components(getlocaldofhandler(dh), field)
Ferrite.CellCache(dh::NODDofHandler, flags::UpdateFlags=UpdateFlags()) = Ferrite.CellCache(getlocaldofhandler(dh), flags)
function Ferrite.CellIterator(dh::NODDofHandler, set::Union{Ferrite.IntegerCollection,Nothing}=nothing, flags::UpdateFlags=UpdateFlags())
    return Ferrite.CellIterator(getlocaldofhandler(dh), set, flags)
end
Ferrite.CellIterator(dh::NODDofHandler, flags::UpdateFlags) = Ferrite.CellIterator(getlocaldofhandler(dh), nothing, flags)
Ferrite.evaluate_at_grid_nodes(dh::NODDofHandler, u::AbstractVector, fieldname::Symbol) = Ferrite.evaluate_at_grid_nodes(getlocaldofhandler(dh), u, fieldname)
Ferrite._evaluate_at_grid_nodes(dh::NODDofHandler, u::AbstractVector, fieldname::Symbol, vtk=Val(false)) = Ferrite._evaluate_at_grid_nodes(getlocaldofhandler(dh), u, fieldname, vtk)

# TODO problem here is that the reorder has to be synchronized. We also cannot arbitrarily
# reorder dofs, because some distributed matrix data structures have strict requirements
# on the orderings.
Ferrite.renumber!(::NODDofHandler, ::AbstractVector{<:Integer}) = error("Not implemented.")

# A boundary set may not intersect the local part of the domain; unlike in the serial
# case this is not an error.
function Ferrite.add!(ch::ConstraintHandler{<:NODDofHandler}, dbc::Dirichlet)
    isempty(dbc.facets) && return ch
    return invoke(Ferrite.add!, Tuple{ConstraintHandler, Dirichlet}, ch, dbc)
end

function Base.show(io::IO, mime::MIME"text/plain", dh::NODDofHandler)
    println(io, "NODDofHandler (rank $(global_rank(getglobalgrid(dh))) of $(global_nranks(getglobalgrid(dh)))) wrapping:")
    show(io, mime, getlocaldofhandler(dh))
end

"""
Compute the global dof range of the dofs owned by the calling process. It is guaranteed to
be continuous. Empty (but correctly positioned) if the process does not own any dofs.
"""
function local_dof_range(dh::NODDofHandler)
    offset = getfield(dh, :gdof_offset)
    return (offset + 1):(offset + num_local_true_dofs(dh))
end

"""
Compute the number of dofs owned by the current process.
"""
num_local_true_dofs(dh::NODDofHandler) = count(==(global_rank(getglobalgrid(dh))), dh.ldof_to_rank)

"""
Compute the number of dofs visible to the current process.
"""
num_local_dofs(dh::NODDofHandler) = length(dh.ldof_to_gdof)

"""
Compute the number of dofs in the global system.
"""
num_global_dofs(dh::NODDofHandler) = MPI.Allreduce(num_local_true_dofs(dh), MPI.SUM, global_comm(getglobalgrid(dh)))

# ------------------------------------------------------------------------------
#                          Entity dof queries
# ------------------------------------------------------------------------------

function _entity_key(grid::Ferrite.AbstractGrid, kind::Int, cell::Int, idx::Int)
    cellgeo = getcells(grid, cell)
    kind == VERTEX_KIND && return Ferrite.vertices(cellgeo)[idx]
    kind == EDGE_KIND && return Ferrite.sortedge_fast(Ferrite.edges(cellgeo)[idx])
    return Ferrite.sortface_fast(Ferrite.faces(cellgeo)[idx])
end

function _entity_first_dof(ed::EntityDofInfo, kind::Int, key, field_idx::Int)
    kind == VERTEX_KIND && return ed.vertexdicts[field_idx][key::Int]
    kind == EDGE_KIND && return get(ed.edgedicts[field_idx], key::NTuple{2,Int}, 0)
    return get(ed.facedicts[field_idx], key::NTuple{3,Int}, 0)
end

function _entity_dofs(dh::NODDofHandler, kind::Int, field_idx::Int, cell::Int, idx::Int)
    ed = entity_dof_info(dh)
    key = _entity_key(get_grid(dh), kind, cell, idx)
    first_dof = _entity_first_dof(ed, kind, key, field_idx)
    first_dof == 0 && return 1:0
    return first_dof:(first_dof + ed.nentitydofs[kind, field_idx] - 1)
end

"""
Compute the dofs belonging to a given vertex of a given field.
"""
vertex_dofs(dh::NODDofHandler, field_idx::Int, vertex::VertexIndex) = _entity_dofs(dh, VERTEX_KIND, field_idx, vertex[1], vertex[2])

"""
Compute the dofs belonging to the interior of a given edge of a given field.
"""
edge_dofs(dh::NODDofHandler, field_idx::Int, edge::EdgeIndex) = _entity_dofs(dh, EDGE_KIND, field_idx, edge[1], edge[2])

"""
Compute the dofs belonging to the interior of a given face of a given field.
"""
face_dofs(dh::NODDofHandler, field_idx::Int, face::FaceIndex) = _entity_dofs(dh, FACE_KIND, field_idx, face[1], face[2])

"""
Compute the dofs belonging to the interior of a given cell of a given field.
"""
function cell_dofs(dh::NODDofHandler, field_idx::Int, cell::Int)
    sdh_idx = dh.cell_to_subdofhandler[cell]
    sdh_idx == 0 && return 1:0
    sdh = dh.subdofhandlers[sdh_idx]
    lidx = Ferrite._find_field(sdh, dh.field_names[field_idx])
    lidx === nothing && return 1:0
    info = Ferrite.InterpolationInfo(sdh.field_interpolations[lidx])
    n = info.nvolumedofs * info.n_copies
    n == 0 && return 1:0
    r = Ferrite.dof_range(sdh, lidx)
    return Ferrite.celldofs(dh, cell)[r[(end - n + 1):end]]
end

has_vertex_dofs(dh::NODDofHandler, field_idx::Int, vertex::VertexIndex) = !isempty(vertex_dofs(dh, field_idx, vertex))
has_edge_dofs(dh::NODDofHandler, field_idx::Int, edge::EdgeIndex) = !isempty(edge_dofs(dh, field_idx, edge))
has_face_dofs(dh::NODDofHandler, field_idx::Int, face::FaceIndex) = !isempty(face_dofs(dh, field_idx, face))
has_cell_dofs(dh::NODDofHandler, field_idx::Int, cell::Int) = !isempty(cell_dofs(dh, field_idx, cell))

# ------------------------------------------------------------------------------
#                     Distributed dof ownership and numbering
# ------------------------------------------------------------------------------

function Ferrite.close!(dh::NODDofHandler)
    ldh = getlocaldofhandler(dh)
    _, vertexdicts, edgedicts, facedicts = Ferrite.__close!(ldh)
    nranks = global_nranks(getglobalgrid(dh))
    counts = _entity_dof_counts(ldh, nranks)
    setfield!(dh, :entity_dofs, EntityDofInfo(vertexdicts, edgedicts, facedicts, counts))
    _distribute_global_dofs!(dh)
    return dh
end

# Number of dofs per (entity kind, field) block, including all components. The counts must
# agree between all SubDofHandlers sharing a field, otherwise dof blocks on subdomain
# interfaces are ambiguous.
function _entity_dof_counts(ldh::Ferrite.DofHandler, nranks::Int)
    nfields = length(ldh.field_names)
    counts = fill(-1, 3, nfields)
    for sdh in ldh.subdofhandlers
        for (lidx, name) in pairs(sdh.field_names)
            gidx = findfirst(==(name), ldh.field_names)::Int
            info = Ferrite.InterpolationInfo(sdh.field_interpolations[lidx])
            nv = _uniform_entity_count(info.nvertexdofs, name)
            ne = _uniform_entity_count(info.nedgedofs, name)
            nf = _uniform_entity_count(info.nfacedofs, name)
            if nranks > 1 && (nv > 1 || ne > 1)
                error("Field :$name has more than one dof per vertex or edge; the orientation handling required for this is not implemented for distributed grids.")
            end
            c = (nv * info.n_copies, ne * info.n_copies, nf * info.n_copies)
            for kind in 1:3
                if counts[kind, gidx] == -1
                    counts[kind, gidx] = c[kind]
                elseif counts[kind, gidx] != c[kind]
                    error("Field :$name has different numbers of entity dofs in different SubDofHandlers; this is not supported for distributed grids.")
                end
            end
        end
    end
    replace!(counts, -1 => 0)
    return counts
end

function _uniform_entity_count(v::Vector{Int}, name::Symbol)
    m = maximum(v; init=0)
    all(x -> x == 0 || x == m, v) || error("Field :$name has a varying number of dofs per entity; this is not supported for distributed grids.")
    return m
end

# All ranks must agree on field indices during communication, but ranks may only know a
# subset of the fields (a subdomain might not intersect every partition). Build a canonical
# field name list, identical on all ranks.
function _synchronize_field_names(dh::NODDofHandler)
    comm = global_comm(getglobalgrid(dh))
    buf = Vector{UInt8}(codeunits(join(string.(Ferrite.getfieldnames(dh)), '\n')))
    lengths = MPI.Allgather(Int32(length(buf)), comm)
    vbuf = MPI.VBuffer(Vector{UInt8}(undef, sum(lengths)), lengths)
    MPI.Allgatherv!(buf, vbuf, comm)
    names = Symbol[]
    offset = 0
    for l in lengths
        chunk = String(vbuf.data[(offset + 1):(offset + l)])
        offset += l
        for s in split(chunk, '\n'; keepempty=false)
            sym = Symbol(s)
            sym in names || push!(names, sym)
        end
    end
    return names
end

# Merge the (cell, idx)-keyed shared entities of the grid by their physical key (the local
# node ids) and record one remote representative (cell, idx) per remote rank.
function _merge_shared_entities(kind::Int, shared_entities, grid::Ferrite.AbstractGrid, ::Type{K}) where {K}
    merged = Dict{K,Dict{Int,NTuple{2,Int}}}()
    for se in shared_entities
        (cell, idx) = se.local_idx
        key = _entity_key(grid, kind, cell, idx)::K
        remotes = get!(Dict{Int,NTuple{2,Int}}, merged, key)
        for (rank, remote_idxs) in remote_entities(se)
            haskey(remotes, rank) && continue
            ri = first(remote_idxs)
            remotes[rank] = (ri[1], ri[2])
        end
    end
    return merged
end

function _distribute_global_dofs!(dh::NODDofHandler)
    dgrid = getglobalgrid(dh)
    grid = getlocalgrid(dh)
    my_rank = global_rank(dgrid)
    nfields = num_fields(dh)
    nldofs = ndofs(dh)
    ed = entity_dof_info(dh)

    ldof_to_rank = getfield(dh, :ldof_to_rank)
    ldof_to_gdof = getfield(dh, :ldof_to_gdof)
    resize!(ldof_to_rank, nldofs)
    fill!(ldof_to_rank, my_rank)
    resize!(ldof_to_gdof, nldofs)
    fill!(ldof_to_gdof, 0)

    ic = InterfaceCommunicator(dgrid)

    canonical_names = _synchronize_field_names(dh)
    local_names = Ferrite.getfieldnames(dh)
    local_to_canonical = Int[findfirst(==(name), canonical_names) for name in local_names]
    canonical_to_local = Int[something(findfirst(==(name), local_names), 0) for name in canonical_names]

    shared = (
        _merge_shared_entities(VERTEX_KIND, get_shared_vertices(dgrid), grid, Int),
        _merge_shared_entities(EDGE_KIND, get_shared_edges(dgrid), grid, NTuple{2,Int}),
        _merge_shared_entities(FACE_KIND, get_shared_faces(dgrid), grid, NTuple{3,Int}),
    )

    # Round 1: announce to all sharing neighbors which fields have dofs on each shared
    # entity. Entities are addressed in the receiver's local (cell, idx) indexing and
    # fields by their canonical index.
    presence_send = empty_send_buffers(Int, ic)
    for kind in 1:3, (key, remotes) in shared[kind]
        for field_idx in 1:nfields
            _entity_first_dof(ed, kind, key, field_idx) == 0 && continue
            cfield = local_to_canonical[field_idx]
            for (rank, (rcell, ridx)) in remotes
                push!(presence_send[ic.destination_index[rank]], kind, rcell, ridx, cfield)
            end
        end
    end
    presence_recv = exchange(ic, presence_send)

    # presence[kind][(key, cfield)] -> remote ranks with dofs for cfield on the entity
    presence = (
        Dict{Tuple{Int,Int},Vector{Int}}(),
        Dict{Tuple{NTuple{2,Int},Int},Vector{Int}}(),
        Dict{Tuple{NTuple{3,Int},Int},Vector{Int}}(),
    )
    for (si, buf) in enumerate(presence_recv)
        rank = ic.sources[si]
        for i in 1:4:length(buf)
            kind, cell, idx, cfield = buf[i], buf[i + 1], buf[i + 2], buf[i + 3]
            if kind == VERTEX_KIND
                key = _entity_key(grid, kind, cell, idx)::Int
                push!(get!(Vector{Int}, presence[VERTEX_KIND], (key, cfield)), rank)
            elseif kind == EDGE_KIND
                key = _entity_key(grid, kind, cell, idx)::NTuple{2,Int}
                push!(get!(Vector{Int}, presence[EDGE_KIND], (key, cfield)), rank)
            else
                key = _entity_key(grid, kind, cell, idx)::NTuple{3,Int}
                push!(get!(Vector{Int}, presence[FACE_KIND], (key, cfield)), rank)
            end
        end
    end

    # The owner of an (entity, field) dof block is the lowest rank with dofs for the field
    # on the entity. Note that this can differ between fields on the same entity when
    # fields are restricted to subdomains.
    _block_owner(kind, key, cfield) = minimum(get(presence[kind], (key, cfield), Int[]); init=my_rank)

    for kind in 1:3, (key, _) in shared[kind]
        for field_idx in 1:nfields
            first_dof = _entity_first_dof(ed, kind, key, field_idx)
            first_dof == 0 && continue
            owner = _block_owner(kind, key, local_to_canonical[field_idx])
            for d in first_dof:(first_dof + ed.nentitydofs[kind, field_idx] - 1)
                ldof_to_rank[d] = owner
            end
        end
    end

    # Assign global numbers to the owned dofs, in ascending local dof order, shifted by the
    # total number of dofs owned by lower ranks.
    num_owned = count(==(my_rank), ldof_to_rank)
    offset = MPI.Exscan(num_owned, +, global_comm(dgrid))
    my_rank == 1 && (offset = 0)
    setfield!(dh, :gdof_offset, offset)
    next_gdof = offset
    for d in 1:nldofs
        if ldof_to_rank[d] == my_rank
            next_gdof += 1
            ldof_to_gdof[d] = next_gdof
        end
    end

    # Round 2: the owner of each shared (entity, field) dof block sends the global number
    # of the block to all neighbors with dofs for the field on the entity.
    sync_send = empty_send_buffers(Int, ic)
    for kind in 1:3, (key, remotes) in shared[kind]
        for field_idx in 1:nfields
            first_dof = _entity_first_dof(ed, kind, key, field_idx)
            first_dof == 0 && continue
            cfield = local_to_canonical[field_idx]
            _block_owner(kind, key, cfield) == my_rank || continue
            n = ed.nentitydofs[kind, field_idx]
            for rank in get(presence[kind], (key, cfield), Int[])
                (rcell, ridx) = remotes[rank]
                push!(sync_send[ic.destination_index[rank]], kind, rcell, ridx, cfield, ldof_to_gdof[first_dof], n)
            end
        end
    end
    sync_recv = exchange(ic, sync_send)

    for buf in sync_recv
        for i in 1:6:length(buf)
            kind, cell, idx, cfield, gdof, n = buf[i], buf[i + 1], buf[i + 2], buf[i + 3], buf[i + 4], buf[i + 5]
            field_idx = canonical_to_local[cfield]
            key = _entity_key(grid, kind, cell, idx)
            first_dof = _entity_first_dof(ed, kind, key, field_idx)
            @assert first_dof != 0
            n == ed.nentitydofs[kind, field_idx] || error("Dof block size mismatch on the process boundary for field :$(canonical_names[cfield]); the interpolations on the neighboring processes are incompatible.")
            for k in 0:(n - 1)
                ldof_to_gdof[first_dof + k] = gdof + k
            end
        end
    end

    @assert all(!=(0), ldof_to_gdof)
    return dh
end
