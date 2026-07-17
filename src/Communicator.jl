"""
    InterfaceCommunicator(dgrid::AbstractNODGrid)

Communication channel along the process boundary of a non-overlapping distributed grid.
Wraps the neighborhood graph communicator of the grid and provides sparse all-to-all
exchange primitives between interface neighbors.

The neighborhood graph is required to be symmetric, so replies to received messages can
be routed back via `destination_index[sources[i]]`.
"""
struct InterfaceCommunicator
    comm::MPI.Comm
    sources::Vector{Int}       # 1-based ranks we receive from
    destinations::Vector{Int}  # 1-based ranks we send to
    source_index::Dict{Int,Int}
    destination_index::Dict{Int,Int}
end

function InterfaceCommunicator(dgrid::AbstractNODGrid)
    comm = interface_comm(dgrid)
    (source_len, destination_len, _) = MPI.Dist_graph_neighbors_count(comm)
    sources = Vector{Cint}(undef, source_len)
    destinations = Vector{Cint}(undef, destination_len)
    MPI.Dist_graph_neighbors!(comm, sources, destinations)
    sources = Int.(sources) .+ 1
    destinations = Int.(destinations) .+ 1
    issetequal(sources, destinations) || error("The interface communicator requires a symmetric neighborhood graph.")
    return InterfaceCommunicator(
        comm, sources, destinations,
        Dict{Int,Int}(r => i for (i, r) in enumerate(sources)),
        Dict{Int,Int}(r => i for (i, r) in enumerate(destinations)),
    )
end

nsources(ic::InterfaceCommunicator) = length(ic.sources)
ndestinations(ic::InterfaceCommunicator) = length(ic.destinations)

"""
    empty_send_buffers(T, ic::InterfaceCommunicator)

Create one empty `Vector{T}` per destination rank, for use with [`exchange`](@ref).
"""
empty_send_buffers(::Type{T}, ic::InterfaceCommunicator) where {T} = [T[] for _ in 1:ndestinations(ic)]

"""
    exchange(ic::InterfaceCommunicator, send_data::Vector{Vector{T}}) -> Vector{Vector{T}}

Sparse all-to-all exchange between interface neighbors. `send_data[i]` is sent to rank
`ic.destinations[i]` and the result contains one buffer per rank in `ic.sources`.
The order of the elements within each buffer is preserved.
"""
function exchange(ic::InterfaceCommunicator, send_data::Vector{Vector{T}}) where {T}
    @assert length(send_data) == ndestinations(ic)
    send_lengths = Cint[length(b) for b in send_data]
    recv_lengths = Vector{Cint}(undef, nsources(ic))
    MPI.Neighbor_alltoall!(UBuffer(send_lengths, 1), UBuffer(recv_lengths, 1), ic.comm)
    send_buffer = reduce(vcat, send_data; init = T[])
    recv_buffer = Vector{T}(undef, sum(recv_lengths; init = 0))
    MPI.Neighbor_alltoallv!(VBuffer(send_buffer, send_lengths), VBuffer(recv_buffer, recv_lengths), ic.comm)
    recv_data = Vector{Vector{T}}(undef, nsources(ic))
    offset = 0
    for (i, l) in enumerate(recv_lengths)
        recv_data[i] = recv_buffer[(offset + 1):(offset + l)]
        offset += l
    end
    return recv_data
end
