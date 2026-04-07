"""
    PVTKGridFile(filename::String, dh::NODDofHandler; kwargs...)
    PVTKGridFile(filename::String, dgrid::AbstractNODGrid; kwargs...)

Create a parallel VTK file for a distributed grid, analogous to Ferrite's `VTKGridFile`.
Use with [`write_solution`](@ref), [`write_cell_data`](@ref), [`write_node_data`](@ref),
[`vtk_shared_vertices`](@ref), [`vtk_shared_faces`](@ref), [`vtk_shared_edges`](@ref),
and [`vtk_partitioning`](@ref).

```julia
PVTKGridFile(filename, dh) do vtk
    write_solution(vtk, dh, u)
    vtk_partitioning(vtk, getglobalgrid(dh))
end
```
"""
struct PVTKGridFile{VTK <: WriteVTK.PVTKFile}
    vtk::VTK
end

function PVTKGridFile(filename::String, dgrid::AbstractNODGrid{dim}; compress::Bool=true) where {dim}
    part   = global_rank(dgrid)
    nparts = global_nranks(dgrid)
    cls = MeshCell[]
    for cell in getcells(dgrid)
        celltype = Ferrite.cell_to_vtkcell(typeof(cell))
        push!(cls, MeshCell(celltype, Ferrite.nodes_to_vtkorder(cell)))
    end
    coords = reshape(reinterpret(get_coordinate_eltype(dgrid), getnodes(dgrid)), (dim, getnnodes(dgrid)))
    pvtk = pvtk_grid(filename, coords, cls; part=part, nparts=nparts, compress=compress)
    return PVTKGridFile(pvtk)
end

PVTKGridFile(filename::String, dh::NODDofHandler; kwargs...) = PVTKGridFile(filename, dh.grid; kwargs...)

function PVTKGridFile(f::Function, args...; kwargs...)
    vtk = PVTKGridFile(args...; kwargs...)
    try
        f(vtk)
    finally
        close(vtk)
    end
    return vtk
end

function Base.close(vtk::PVTKGridFile)
    WriteVTK.vtk_save(vtk.vtk)
    return vtk
end

"""
    write_solution(vtk::PVTKGridFile, dh::NODDofHandler, u::AbstractVector, suffix="")

Save the values at the nodes in the degree of freedom vector `u` to `vtk`.
Each field in `dh` will be saved separately, and `suffix` can be used to append
to the fieldname.
"""
function Ferrite.write_solution(vtk::PVTKGridFile, dh::NODDofHandler, u::AbstractVector, suffix="")
    # TODO convert NODDofHandler to DofHandler before write
    for name in getfieldnames(dh)
        data = Ferrite._evaluate_at_grid_nodes(dh, u, name, Val(true))
        Ferrite._vtk_write_node_data(vtk.vtk.vtk, data, string(name, suffix))
    end
    return vtk
end

"""
    write_cell_data(vtk::PVTKGridFile, celldata::AbstractVector, name::String)

Write the `celldata` that is ordered by the cells in the grid to the vtk file.
"""
function Ferrite.write_cell_data(vtk::PVTKGridFile, celldata, name)
    vtk_cell_data(vtk.vtk.vtk, celldata, name)
    return vtk
end

"""
    write_node_data(vtk::PVTKGridFile, nodedata, name)

Write the `nodedata` that is ordered by the nodes in the grid to `vtk`.
"""
function Ferrite.write_node_data(vtk::PVTKGridFile, nodedata, name)
    Ferrite._vtk_write_node_data(vtk.vtk.vtk, nodedata, name)
    return vtk
end

"""
    vtk_shared_vertices(vtk::PVTKGridFile, dgrid::AbstractNODGrid)

Enrich the VTK file with meta information about shared vertices.
"""
function vtk_shared_vertices(vtk::PVTKGridFile, dgrid::AbstractNODGrid)
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
        vtk_point_data(vtk.vtk.vtk, u, "shared vertices with $rank")
    end
end

"""
    vtk_shared_faces(vtk::PVTKGridFile, dgrid::AbstractNODGrid)

Enrich the VTK file with meta information about shared faces.
"""
function vtk_shared_faces(vtk::PVTKGridFile, dgrid::AbstractNODGrid)
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
        vtk_point_data(vtk.vtk.vtk, u, "shared faces with $rank")
    end
end

"""
    vtk_shared_edges(vtk::PVTKGridFile, dgrid::AbstractNODGrid)

Enrich the VTK file with meta information about shared edges.
"""
function vtk_shared_edges(vtk::PVTKGridFile, dgrid::AbstractNODGrid)
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
        vtk_point_data(vtk.vtk.vtk, u, "shared edges with $rank")
    end
end

"""
    vtk_partitioning(vtk::PVTKGridFile, dgrid::AbstractNODGrid)

Enrich the VTK file with partitioning meta information.
"""
function vtk_partitioning(vtk::PVTKGridFile, dgrid::AbstractNODGrid)
    u  = Vector{Float64}(undef, getncells(dgrid))
    u .= global_rank(dgrid)
    vtk_cell_data(vtk.vtk.vtk, u, "partitioning")
end
