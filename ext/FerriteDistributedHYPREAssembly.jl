"""
Module containing the code for distributed assembly via HYPRE.jl
"""
module FerriteDistributedHYPREAssembly

using FerriteDistributed
# TODO remove me. These are merely hotfixes to split the extensions trasiently via an internal API.
import FerriteDistributed: getglobalgrid, num_local_true_dofs, num_local_dofs, global_comm, interface_comm, global_rank, compute_owner, remote_entities
using MPI
using HYPRE
using Base: @propagate_inbounds

include("FerriteDistributedHYPREAssembly/assembler.jl")
include("FerriteDistributedHYPREAssembly/conversion.jl")

function __init__()
    @info "FerriteHYPRE extension loaded."
end

end # module FerriteHYPRE
