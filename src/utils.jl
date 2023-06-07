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


# ------------------------------------
#         DOF UTILITY BLOCK
# ------------------------------------

#TODO move this into ferrite core
has_cell_dofs(dh::Ferrite.AbstractDofHandler, field_idx::Int, cell::Int) = ncelldofs(Ferrite.getfieldinterpolation(dh, field_idx)) > 0
has_vertex_dofs(dh::Ferrite.AbstractDofHandler, field_idx::Int, vertex::VertexIndex) = nvertexdofs(Ferrite.getfieldinterpolation(dh, field_idx)) > 0
has_edge_dofs(dh::Ferrite.AbstractDofHandler, field_idx::Int, edge::EdgeIndex) = nedgedofs(Ferrite.getfieldinterpolation(dh, field_idx)) > 0
has_face_dofs(dh::Ferrite.AbstractDofHandler, field_idx::Int, face::FaceIndex) = nfacedofs(Ferrite.getfieldinterpolation(dh, field_idx)) > 0

# entity_dofs(dh::Ferrite.AbstractDofHandler, field_idx::Int, vertex::Int) = vertexdicts[field_idx][vertex]
# entity_dofs(dh::Ferrite.AbstractDofHandler, field_idx::Int, edge::Tuple{Int,Int}) = edgedicts[field_idx][edge]
# entity_dofs(dh::Ferrite.AbstractDofHandler, field_idx::Int, face::NTuple{dim,Int}) where {dim} = facedicts[field_idx][face]

"""
Compute the dofs belonging to a given cell of a given field.
"""
function cell_dofs(dh::Ferrite.AbstractDofHandler, field_idx::Int, cell::Int)
    ip = Ferrite.getfieldinterpolation(dh, field_idx)
    fdim = getfielddim(dh, field_idx)
    nentitydofs = fdim*ncelldofs(ip)
    totaldofs = fdim*getnbasefunctions(ip)
    ldofs = dof_range(dh, field_idx)[(totaldofs-nentitydofs+1):totaldofs]
    return celldofs(dh, cell)[ldofs]
end

nvertexdofs(ip::Interpolation) = length(Ferrite.vertexdof_indices(ip)[1])
nfacedofs(ip::Interpolation) = length(Ferrite.facedof_interior_indices(ip)[1])
nedgedofs(ip::Interpolation) = length(Ferrite.edgedof_interior_indices(ip)[1])
ncelldofs(ip::Interpolation) = length(Ferrite.celldof_interior_indices(ip)[1])

"""
Compute the dofs belonging to a given vertex of a given field.
"""
function vertex_dofs(dh::Ferrite.AbstractDofHandler, field_idx::Int, vertex::VertexIndex)
    ip = Ferrite.getfieldinterpolation(dh, field_idx)
    nvdofs = nvertexdofs(ip)
    nvdofs == 0 && return Int[]
    fdim = getfielddim(dh, field_idx)
    cell,local_vertex_index = vertex
    cell_geo = getcells(getgrid(dh), cell)
    nvertices = length(vertices(cell_geo))
    nentitydofs = fdim*nvdofs*nvertices
    ldofr = dof_range(dh, field_idx)[1:nentitydofs]
    vdofs = Ferrite.celldofs(dh, cell)[ldofr]
    return reshape(vdofs, (fdim,nvertices))[:, local_vertex_index]
end

"""
Compute the dofs belonging to a given edge of a given field.
"""
function edge_dofs(dh::Ferrite.AbstractDofHandler, field_idx::Int, edge::EdgeIndex)
    ip = Ferrite.getfieldinterpolation(dh, field_idx)
    nedofs = nedgedofs(ip)
    nedofs == 0 && return Int[]
    nvdofs = nvertexdofs(ip)
    fdim = getfielddim(dh, field_idx)
    cell,local_edge_index = edge
    cell_geo = getcells(getgrid(dh), cell)
    nedges_on_cell = length(edges(cell_geo))
    nvertices_on_cell = length(vertices(cell_geo))
    nentitydofs = fdim*nedofs*nedges_on_cell
    offset = fdim*nvdofs*nvertices_on_cell
    edge_dofrange = dof_range(dh, field_idx)[(offset+1):(offset+nentitydofs)]
    lodal_edgedofs = Ferrite.celldofs(dh, cell)[edge_dofrange]
    return reshape(lodal_edgedofs, (fdim,nedges_on_cell))[:, local_edge_index]
end

"""
Compute the dofs belonging to a given face of a given field.
"""
function face_dofs(dh::Ferrite.AbstractDofHandler, field_idx::Int, face::FaceIndex)
    ip = Ferrite.getfieldinterpolation(dh, field_idx)
    dim = getdim(dh)
    nfdofs = nfacedofs(ip)
    nfdofs == 0 && return Int[]
    nvdofs = nvertexdofs(ip)
    fdim = getfielddim(dh, field_idx)
    cell,local_face_index = face
    cell_geo = getcells(getgrid(dh), cell)
    nfaces_on_cell = length(faces(cell_geo))
    nvertices_on_cell = length(vertices(cell_geo))
    nentitydofs = fdim*nfacedofs(ip)*nfaces_on_cell
    offset = fdim*nvdofs*nvertices_on_cell
    if dim > 2
        nedges_on_cell = length(edges(cell_geo))
        nedofs = nedgedofs(ip)
        offset += fdim*nedofs*nedges_on_cell
    end
    face_dofrange = dof_range(dh, field_idx)[(offset+1):(offset+nentitydofs)]
    local_facedofs = Ferrite.celldofs(dh, cell)[face_dofrange]
    return reshape(local_facedofs, (fdim,nfaces_on_cell))[:, local_face_index]
end
