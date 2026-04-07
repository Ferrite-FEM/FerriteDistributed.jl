using FerriteDistributed
using PartitionedArrays, Metis
using IterativeSolvers

MPI.Init()

dgrid = generate_nod_grid(MPI.COMM_WORLD, Hexahedron, (10, 10, 10); partitioning_alg=FerriteDistributed.PartitioningAlgorithm.Metis(:RECURSIVE))

ref = RefHexahedron
ip = Lagrange{ref, 2}()
ip_geo = Lagrange{ref, 1}()
qr = QuadratureRule{ref}(3)
cellvalues = CellValues(qr, ip, ip_geo)

dh = DofHandler(dgrid)
add!(dh, :u, ip)
close!(dh)

ch = ConstraintHandler(dh);
∂Ω = union(getfaceset.((dgrid, ), ["left", "right", "top", "bottom", "front", "back"])...);
dbc = Dirichlet(:u, ∂Ω, (x, t) -> 0)
add!(ch, dbc)
close!(ch)

function assemble_element!(Ke, fe, cellvalues, cell_coords::AbstractVector{<:Vec{dim}}) where dim
    fill!(Ke, 0)
    fill!(fe, 0)

    n_basefuncs = getnbasefunctions(cellvalues)
    for q_point in 1:getnquadpoints(cellvalues)
        dΩ = getdetJdV(cellvalues, q_point)

        for i in 1:n_basefuncs
            v  = shape_value(cellvalues, q_point, i)
            ∇v = shape_gradient(cellvalues, q_point, i)
            # Manufactured solution of Π cos(xᵢπ)
            x = spatial_coordinate(cellvalues, q_point, cell_coords)
            fe[i] += (π/2)^2 * dim * prod(cos, x*π/2) * v * dΩ

            for j in 1:n_basefuncs
                ∇u = shape_gradient(cellvalues, q_point, j)
                Ke[i, j] += (∇v ⋅ ∇u) * dΩ
            end
        end
    end
end

function doassemble(cellvalues::CellValues, dh::FerriteDistributed.NODDofHandler{dim}) where {dim}
    assembler = start_assemble(dh, distribute_with_mpi(LinearIndices((MPI.Comm_size(MPI.COMM_WORLD),))))

    # For the local assembly nothing changes
    n_basefuncs = getnbasefunctions(cellvalues)
    Ke = zeros(n_basefuncs, n_basefuncs)
    fe = zeros(n_basefuncs)
    for cell in CellIterator(dh)
        fill!(Ke, 0)
        fill!(fe, 0)
        reinit!(cellvalues, cell)
        assemble_element!(Ke, fe, cellvalues, getcoordinates(cell))
        # Note that this call should be communication-free!
        Ferrite.assemble!(assembler, celldofs(cell), fe, Ke)
    end

    # Finally, for the `PartitionedArraysCOOAssembler` we have to call
    # `end_assemble` to construct the global sparse matrix and the global
    # right hand side vector.
    return end_assemble(assembler)
end

K, f = doassemble(cellvalues, dh)
apply!(K, f, ch)

u = cg(K, f)

vtk_grid("heat_equation_distributed", dh) do vtk
    vtk_point_data(vtk, dh, u)
    vtk_shared_vertices(vtk, dgrid)
    vtk_shared_faces(vtk, dgrid)
    vtk_partitioning(vtk, dgrid)
end

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
