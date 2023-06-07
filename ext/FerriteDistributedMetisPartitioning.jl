
module FerriteDistributedMetisPartitioning

# This extension requires modules as type parameters
# https://github.com/JuliaLang/julia/pull/47749
if VERSION >= v"1.10.0-DEV.90"

using FerriteDistributed, Metis
import FerriteDistributed: CoverTopology, PartitioningAlgorithm

struct MetisPartitioning <: PartitioningAlgorithm.Ext{Metis}
    alg::Symbol
end

PartitioningAlgorithm.Ext{Metis}(alg::Symbol) = MetisPartitioning(alg)

function FerriteDistributed.create_partitioning(grid::Grid{dim,C,T}, grid_topology::CoverTopology, n_partitions::Int, partition_alg::PartitioningAlgorithm.Ext{Metis}) where {dim,C,T}
    n_cells_global = getncells(grid)
    @assert n_cells_global > 0

    if n_partitions == 1
        return ones(Metis.idx_t, n_cells_global)
    end

    # Set up the element connectivity graph
    xadj = Vector{Metis.idx_t}(undef, n_cells_global+1)
    xadj[1] = 1
    adjncy = Vector{Metis.idx_t}(undef, 0)
    @inbounds for i in 1:n_cells_global
        n_neighbors = 0
        for neighbor ∈ getneighborhood(grid_topology, grid, CellIndex(i))
            push!(adjncy, neighbor)
            n_neighbors += 1
        end
        xadj[i+1] = xadj[i] + n_neighbors
    end

    # Generate a partitioning
    return Metis.partition(
        Metis.Graph(
            Metis.idx_t(n_cells_global),
            xadj,
            adjncy
        ),
        n_partitions;
        alg=partition_alg.alg
    )
end

end # VERSION check

end # module FerriteDistributedMetisPartitioning
