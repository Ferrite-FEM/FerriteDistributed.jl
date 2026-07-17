"""
Module containing the code for distributed assembly via PartitionedArrays.jl
"""
module FerriteDistributedPartitionedArrays

using FerriteDistributed
import FerriteDistributed: getglobalgrid, num_global_dofs, num_local_dofs, global_comm, global_rank, num_fields,
    InterfaceCommunicator, exchange, empty_send_buffers
using MPI
using PartitionedArrays
using Base: @propagate_inbounds

include("FerriteDistributedPartitionedArrays/assembler.jl")
include("FerriteDistributedPartitionedArrays/constraints.jl")

function __init__()
    @info "FerriteDistributedPartitionedArrays extension loaded."
end

end # module FerriteDistributedPartitionedArrays
