# ------------------------------------
#         GRID UTILITY BLOCK
# ------------------------------------

"""
    toglobal(grid::AbstractGrid, vertexidx::FaceIndex) -> Int
    toglobal(grid::AbstractGrid, vertexidx::Vector{FaceIndex}) -> Vector{Tuple{Int}}
This function takes the local face representation (a `FaceIndex`) and looks up the unique global id (a tuple of `Int`).
"""
Ferrite.toglobal(grid::Ferrite.AbstractGrid,faceidx::Ferrite.FaceIndex) = Ferrite.sortface(faces(getcells(grid,faceidx[1]))[faceidx[2]])
Ferrite.toglobal(grid::Ferrite.AbstractGrid,faceidx::Vector{FaceIndex}) = unique(Ferrite.toglobal.((grid,),faceidx))

"""
    toglobal(grid::AbstractGrid, vertexidx::EdgeIndex) -> Int
    toglobal(grid::AbstractGrid, vertexidx::Vector{EdgeIndex}) -> Vector{Tuple{Int}}
This function takes the local face representation (an `EdgeIndex`) and looks up the unique global id (a tuple of `Int`).
"""
Ferrite.toglobal(grid::Ferrite.AbstractGrid,edgeidx::Ferrite.EdgeIndex) = Ferrite.sortedge(edges(getcells(grid,edgeidx[1]))[edgeidx[2]])[1]
Ferrite.toglobal(grid::Ferrite.AbstractGrid,edgeidx::Vector{Ferrite.EdgeIndex}) = unique(toglobal.((grid,),edgeidx))
