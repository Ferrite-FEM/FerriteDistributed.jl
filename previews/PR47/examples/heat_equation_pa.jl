using FerriteDistributed
using PartitionedArrays, Metis
using IterativeSolvers

MPI.Init()

dgrid = generate_nod_grid(MPI.COMM_WORLD, Hexahedron, (10, 10, 10); partitioning_alg=FerriteDistributed.PartitioningAlgorithm.Metis(:RECURSIVE));

ref = RefHexahedron
ip = Lagrange{ref, 2}()
ip_geo = Lagrange{ref, 1}()
qr = QuadratureRule{ref}(3)
cellvalues = CellValues(qr, ip, ip_geo);

dh = DofHandler(dgrid)
add!(dh, :u, ip)
close!(dh);

ch = ConstraintHandler(dh);
∂Ω = union(getfacetset.((dgrid, ), ["left", "right", "top", "bottom", "front", "back"])...);
dbc = Dirichlet(:u, ∂Ω, (x, t) -> 0)
add!(ch, dbc);
close!(ch)
update!(ch, 0.0);

function doassemble(cellvalues::CellValues, dh::NODDofHandler{dim}) where {dim}
    n_basefuncs = getnbasefunctions(cellvalues)
    Ke = zeros(n_basefuncs, n_basefuncs)
    fe = zeros(n_basefuncs)

    # --------------------- Distributed assembly --------------------
    # The synchronization with the global sparse matrix is handled by
    # an assembler again. You can choose from different backends, which
    # are described in the docs and will be expanded over time. This call
    # may trigger a large amount of communication.
    # NOTE: At the time of writing the only backend available is a COO
    #       assembly via PartitionedArrays.jl .
    assembler = start_assemble(dh, distribute_with_mpi(LinearIndices((MPI.Comm_size(MPI.COMM_WORLD),))))

    # For the local assembly nothing changes
    for cell in CellIterator(dh)
        fill!(Ke, 0)
        fill!(fe, 0)

        reinit!(cellvalues, cell)
        coords = getcoordinates(cell)

        for q_point in 1:getnquadpoints(cellvalues)
            dΩ = getdetJdV(cellvalues, q_point)

            for i in 1:n_basefuncs
                v  = shape_value(cellvalues, q_point, i)
                ∇v = shape_gradient(cellvalues, q_point, i)
                # Manufactured solution of Π cos(xᵢπ)
                x = spatial_coordinate(cellvalues, q_point, coords)
                fe[i] += (π/2)^2 * dim * prod(cos, x*π/2) * v * dΩ

                for j in 1:n_basefuncs
                    ∇u = shape_gradient(cellvalues, q_point, j)
                    Ke[i, j] += (∇v ⋅ ∇u) * dΩ
                end
            end
        end

        # Note that this call should be communication-free!
        Ferrite.assemble!(assembler, celldofs(cell), Ke, fe)
    end

    # Finally, for the `PartitionedArraysCOOAssembler` we have to call
    # `end_assemble` to construct the global sparse matrix and the global
    # right hand side vector.
    return end_assemble(assembler)
end

K, f = doassemble(cellvalues, dh);
apply!(K, f, ch)

u = cg(K, f)

u_local = Vector{Float64}(undef, FerriteDistributed.num_local_dofs(dh))
FerriteDistributed.extract_local_part!(u_local, u, dh);

PVTKGridFile("heat_equation_distributed", dh) do vtk
    write_solution(vtk, dh, u_local)
    # For debugging purposes it can be helpful to enrich
    # the visualization with some meta  information about
    # the grid and its partitioning
    vtk_shared_vertices(vtk, dgrid)
    vtk_shared_faces(vtk, dgrid)
    vtk_partitioning(vtk, dgrid)
end

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
