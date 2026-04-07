module FerriteDistributed

using Reexport
@reexport using Ferrite, MPI
using OrderedCollections: OrderedSet

# Common stuff from core
import Ferrite: get_coordinate_eltype, @debug,
    nnodes_per_cell, n_components, get_grid, getrefdim, getspatialdim,
    BoundaryIndex, FaceIndex, FacetIndex, EdgeIndex, CellIndex, VertexIndex,
    AbstractTopology,
    AbstractCell, boundaryfunction, faces, edges, vertices, nvertices, nfaces, nedges,
    cellnodes!, cellcoords!,
    getfieldnames, getfieldinterpolation, geometric_interpolation,
    reference_coordinates, getrefshape, dof_range,
    BCValues

import Ferrite: WriteVTK
using WriteVTK: MeshCell, pvtk_grid, vtk_point_data, vtk_cell_data

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
