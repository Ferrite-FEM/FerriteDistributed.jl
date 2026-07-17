"""
Module containing the code for distributed assembly via HYPRE.jl
"""
module FerriteDistributedHYPREAssembly

using FerriteDistributed
import FerriteDistributed: getglobalgrid, num_local_true_dofs, global_comm, global_rank,
    InterfaceCommunicator, exchange, empty_send_buffers, local_dof_range
using MPI
using HYPRE
using Base: @propagate_inbounds

include("FerriteDistributedHYPREAssembly/assembler.jl")
include("FerriteDistributedHYPREAssembly/conversion.jl")

function __init__()
    @info "FerriteDistributed HYPRE extension loaded."
end

end # module FerriteHYPRE
