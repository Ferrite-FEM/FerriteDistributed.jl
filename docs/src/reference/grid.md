```@meta
CurrentModule = FerriteDistributed
DocTestSetup = :(using FerriteDistributed)
```

# NODGrid

```@docs
NODGrid
NODGrid(::MPI.Comm, ::Grid)
NODGrid(::MPI.Comm, ::Grid, ::CoverTopology, ::Vector{Int})
FerriteDistributed.SharedEntity
FerriteDistributed.SharedVertex
FerriteDistributed.SharedFace
FerriteDistributed.SharedEdge
FerriteDistributed.remote_entities
FerriteDistributed.global_comm
FerriteDistributed.interface_comm
FerriteDistributed.global_rank
FerriteDistributed.global_nranks
```
