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
∂Ω = union(getfaceset.((dgrid, ), ["left", "right", "top", "bottom", "front", "back"])...);
dbc = Dirichlet(:u, ∂Ω, (x, t) -> 0)
add!(ch, dbc);
close!(ch)
update!(ch, 0.0);

function doassemble(cellvalues::CellValues, dh::FerriteDistributed.NODDofHandler{dim}) where {dim}
    n_basefuncs = getnbasefunctions(cellvalues)
    Ke = zeros(n_basefuncs, n_basefuncs)
    fe = zeros(n_basefuncs)

    assembler = start_assemble(dh, distribute_with_mpi(LinearIndices((MPI.Comm_size(MPI.COMM_WORLD),))))

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

                x = spatial_coordinate(cellvalues, q_point, coords)
                fe[i] += (π/2)^2 * dim * prod(cos, x*π/2) * v * dΩ

                for j in 1:n_basefuncs
                    ∇u = shape_gradient(cellvalues, q_point, j)
                    Ke[i, j] += (∇v ⋅ ∇u) * dΩ
                end
            end
        end

        Ferrite.assemble!(assembler, celldofs(cell), fe, Ke)
    end

    return end_assemble(assembler)
end

K, f = doassemble(cellvalues, dh);
apply!(K, f, ch)

u = cg(K, f)

vtk_grid("heat_equation_distributed", dh) do vtk
    vtk_point_data(vtk, dh, u)

    vtk_shared_vertices(vtk, dgrid)
    vtk_shared_faces(vtk, dgrid)
    vtk_partitioning(vtk, dgrid)
end

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
