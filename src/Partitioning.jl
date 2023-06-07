########################
###   Partitioning   ###
########################

# Namespace for partitioning algorithm selection for the purpose of only exporting the PartitioningAlgorithm name.
# Idea blatantly stolen from Frederik Ekre (Ferrite.jl/src/Dofs/DofRenumbering.jl)
module PartitioningAlgorithm
    """
        PartitioningAlgorithm.SFC

    Simplest space-filling curve based partitioning. Might not scale pretty well.
    """
    struct SFC end

    """
        PartitioningAlgorithm.Metis

    Partitioning algorithm from the extension `FerriteDistributedMetisPartitioning`.
    """
    struct Metis
        alg::Symbol
    end
end # module PartitioningAlgorithm

"""
    create_partitioning(::Ferrite.AbstractGrid, ::Ferrite.AbstractTopology, nparts::Int, partition_alg)::Vector{Int}

Internal entry point for the construction of Ferrite.jl grid partitionings.
Creates a vector where the index corresponds to the cell index of the input grid and the value is the 1-based rank.
"""
create_partitioning(::Ferrite.AbstractGrid, ::Ferrite.AbstractTopology, nparts::Int, partition_alg)

# TODO better balancing.
function create_partitioning(grid::Ferrite.AbstractGrid, ::Ferrite.AbstractTopology, nparts::Int, ::PartitioningAlgorithm.SFC)
    chunk_size = getncells(grid)÷nparts
    parts = [((i-1) ÷ chunk_size) + 1 for i ∈ 1:getncells(grid)]
    parts[(end-rem(getncells(grid), nparts)):end] .= nparts # Assign fat tail for now.
    return parts
end
