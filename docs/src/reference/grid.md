```@meta
CurrentModule = FerriteNODGrid
DocTestSetup = :(using FerriteNODGrid)
```

# NODGrid

```@docs
NODGrid
NODGrid(::MPI.Comm, ::Grid)
NODGrid(::MPI.Comm, ::Grid, ::CoverTopology, ::Vector{Int})
FerriteNODGrid.SharedEntity
FerriteNODGrid.SharedVertex
FerriteNODGrid.SharedFace
FerriteNODGrid.SharedEdge
FerriteNODGrid.remote_entities
FerriteNODGrid.global_comm
FerriteNODGrid.interface_comm
FerriteNODGrid.global_rank
FerriteNODGrid.global_nranks
```
