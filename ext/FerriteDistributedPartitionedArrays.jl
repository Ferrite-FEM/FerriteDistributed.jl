"""
Module containing the code for distributed assembly via PartitionedArrays.jl
"""
module FerriteDistributedPartitionedArrays

using FerriteDistributed
import FerriteDistributed: getglobalgrid, num_global_dofs, num_local_true_dofs, num_local_dofs, global_comm, interface_comm, global_rank, compute_owner, remote_entities,
    num_fields, getfieldnames
using MPI
using PartitionedArrays
using Base: @propagate_inbounds

include("FerriteDistributedPartitionedArrays/assembler.jl")
include("FerriteDistributedPartitionedArrays/constraints.jl")
include("FerriteDistributedPartitionedArrays/export-vtk.jl")

function __init__()
    @info "FerriteDistributedPartitionedArrays extension loaded."
end

end # module FerriteDistributedPartitionedArrays
