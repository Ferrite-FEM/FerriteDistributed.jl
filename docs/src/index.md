```@meta
DocTestSetup = :(using FerriteNODGrid)
```

# FerriteNODGrid.jl

Non-overlapping distributed grids for [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl).

!!! note

    This package is still experimental and may break!

## Usage

The usage of this package is straight forward. Either generate the non-overlapping distributed grid directly via

```julia
using MPI
using FerriteNODGrid
MPI.Init()
dgrid = generate_nod_grid(MPI.COMM_WORLD, Quadrilateral, (100, 100))
```

or distribute your own grid

```julia
using MPI
using Ferrite, FerriteNODGrid
MPI.Init()
grid = ...
dgrid = NODGrid(MPI.COMM_WORLD, grid, FerriteNODGrid.PartitioningAlgorithm.SFC())
```

for more details on the calls please consult the API docs.

## Debugging Information

The debug mode for this package is tied to the [debug mode of Ferrite](https://ferrite-fem.github.io/Ferrite.jl/stable/#Debugging-Information).
Hence, it can be turned on via

```julia
using Ferrite
Ferrite.debug_mode()
```

followed by restarting the Julia process. It can be turned off again by calling

```julia
using Ferrite
Ferrite.debug_mode(enable=false)
```

also followed by restarting the Julia process.
