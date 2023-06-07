module FerriteDistributed

using Reexport
@reexport using Ferrite, MPI

# Common stuff from core
import Ferrite: get_coordinate_eltype, ScalarWrapper,
    nnodes_per_cell, n_components,
    FaceIndex, EdgeIndex, CellIndex, VertexIndex,
    faces, edges, vertices,
    cellnodes!, cellcoords!,
    getfieldnames, getfieldinterpolation, default_interpolation,
    reference_coordinates, value, getrefshape, dof_range, getfielddim

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
