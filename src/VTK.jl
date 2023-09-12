"""
vtk_grid(::AbstractString, ::AbstractNODGrid{dim}; compress::Bool=true)

Store the grid as a PVTK file.
"""
function WriteVTK.vtk_grid(filename::AbstractString, dgrid::AbstractNODGrid{dim}; compress::Bool=true) where {dim}
    part   = global_rank(dgrid)
    nparts = global_nranks(dgrid)
    cls = MeshCell[]
    for cell in getcells(dgrid)
        celltype = Ferrite.cell_to_vtkcell(typeof(cell))
        push!(cls, MeshCell(celltype, Ferrite.nodes_to_vtkorder(cell)))
    end
    coords = reshape(reinterpret(get_coordinate_eltype(dgrid), getnodes(dgrid)), (dim, getnnodes(dgrid)))
    return pvtk_grid(filename, coords, cls; part=part, nparts=nparts, compress=compress)
end

WriteVTK.vtk_grid(filename::AbstractString, dh::NODDofHandler{dim}; kwargs...) where {dim} = vtk_grid(filename, dh.grid; kwargs...)

"""
Enrich the VTK file with meta information about shared vertices.
"""
function vtk_shared_vertices(pvtk::WriteVTK.PVTKFile, dgrid::AbstractNODGrid)
    u = Vector{Float64}(undef, getnnodes(dgrid))
    my_rank = global_rank(dgrid)
    for rank ∈ 1:global_rank(dgrid)
        fill!(u, 0.0)
        for sv ∈ get_shared_vertices(dgrid)
            if haskey(sv.remote_vertices, rank)
                (cellidx, i) = sv.local_idx
                cell = getcells(dgrid, cellidx)
                u[Ferrite.vertices(cell)[i]] = my_rank
            end
        end
        vtk_point_data(pvtk.vtk, u, "shared vertices with $rank")
    end
end

"""
Enrich the VTK file with meta information about shared faces.
"""
function vtk_shared_faces(pvtk::WriteVTK.PVTKFile, dgrid::AbstractNODGrid)
    u = Vector{Float64}(undef, getnnodes(dgrid))
    my_rank = global_rank(dgrid)
    for rank ∈ 1:global_rank(dgrid)
        fill!(u, 0.0)
        for sf ∈ get_shared_faces(dgrid)
            if haskey(sf.remote_faces, rank)
                (cellidx, i) = sf.local_idx
                cell = getcells(dgrid, cellidx)
                facenodes = Ferrite.faces(cell)[i]
                u[[facenodes...]] .= my_rank
            end
        end
        vtk_point_data(pvtk.vtk, u, "shared faces with $rank")
    end
end

"""
Enrich the VTK file with meta information about shared edges.
"""
function vtk_shared_edges(pvtk::WriteVTK.PVTKFile, dgrid::AbstractNODGrid)
    u = Vector{Float64}(undef, getnnodes(dgrid))
    my_rank = global_rank(dgrid)
    for rank ∈ 1:global_rank(dgrid)
        fill!(u, 0.0)
        for se ∈ get_shared_edges(dgrid)
            if haskey(se.remote_edges, rank)
                (cellidx, i) = se.local_idx
                cell = getcells(dgrid, cellidx)
                edgenodes = Ferrite.edges(cell)[i]
                u[[edgenodes...]] .= my_rank
            end
        end
        vtk_point_data(pvtk.vtk, u, "shared edges with $rank")
    end
end

"""
Enrich the VTK file with partitioning meta information.
"""
function vtk_partitioning(pvtk::WriteVTK.PVTKFile, dgrid::AbstractNODGrid)
    u  = Vector{Float64}(undef, getncells(dgrid))
    u .= global_rank(dgrid)
    vtk_cell_data(pvtk.vtk, u, "partitioning")
end

function WriteVTK.vtk_point_data(
    pvtk::WriteVTK.PVTKFile,
    data::Vector{S},
    name::AbstractString
    ) where {O, D, T, M, S <: Union{AbstractFloat, Tensor{O, D, T, M}, SymmetricTensor{O, D, T, M}}}
    return vtk_point_data(pvtk.vtk, data, name)
end

function WriteVTK.vtk_point_data(
    pvtk::WriteVTK.PVTKFile,
    dh::Ferrite.AbstractDofHandler,
    u::Vector,
    suffix="")
    vtk_point_data(pvtk.vtk, dh, u, suffix)
end
