module FerriteDistributed

using Reexport
@reexport using Ferrite, MPI

# Common stuff from core
import Ferrite: get_coordinate_eltype, ScalarWrapper, @debug,
    nnodes_per_cell, n_components, get_grid, getdim,
    BoundaryIndex, FaceIndex, EdgeIndex, CellIndex, VertexIndex,
    AbstractTopology, EntityNeighborhood,
    AbstractCell, boundaryfunction, faces, edges, vertices, nvertices, nfaces, nedges,
    cellnodes!, cellcoords!,
    getfieldnames, getfieldinterpolation, default_interpolation,
    reference_coordinates, value, getrefshape, dof_range, getfielddim,
    BCValues

include("utils.jl")

include("CoverTopology.jl")

include("SharedEntity.jl")

include("interface.jl")

include("Partitioning.jl")

include("NODGrid.jl")
include("NODDofHandler.jl")

include("VTK.jl")

include("exports.jl")

end
