using FerriteDistributed
using HYPRE, Metis

import FerriteDistributed: getglobalgrid, global_comm, local_dof_range #TODO REMOVE THIS

MPI.Init()
HYPRE.Init()

dgrid = generate_nod_grid(MPI.COMM_WORLD, Hexahedron, (10, 10, 10); partitioning_alg=FerriteDistributed.PartitioningAlgorithm.Metis(:RECURSIVE));

ref = RefHexahedron
ip = Lagrange{ref, 2}()
ip_geo = Lagrange{ref, 1}()
qr = QuadratureRule{ref}(3)
cellvalues = CellValues(qr, ip, ip_geo);

dh = DofHandler(dgrid)
push!(dh, :u, 1, ip)
close!(dh);

ch = ConstraintHandler(dh);
∂Ω = union(getfaceset.((dgrid, ), ["left", "right", "top", "bottom", "front", "back"])...);
dbc = Dirichlet(:u, ∂Ω, (x, t) -> 0)
add!(ch, dbc);
close!(ch)

function doassemble(cellvalues::CellValues, dh::FerriteDistributed.NODDofHandler{dim}, ch::ConstraintHandler) where {dim}
    n_basefuncs = getnbasefunctions(cellvalues)
    Ke = zeros(n_basefuncs, n_basefuncs)
    fe = zeros(n_basefuncs)

    dgrid = getglobalgrid(dh)
    comm = global_comm(dgrid)
    ldofrange = local_dof_range(dh)
    K = HYPREMatrix(comm, first(ldofrange), last(ldofrange))
    f = HYPREVector(comm, first(ldofrange), last(ldofrange))

    assembler = start_assemble(K, f)

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

        apply_local!(Ke, fe, celldofs(cell), ch)

        assemble!(assembler, dh.ldof_to_gdof[celldofs(cell)], fe, Ke)
    end

    end_assemble(assembler)

    return K, f
end

K, f = doassemble(cellvalues, dh, ch);

precond = HYPRE.BoomerAMG()
solver = HYPRE.PCG(; Precond = precond)
uₕ = HYPRE.solve(solver, K, f)

u_local = Vector{Float64}(undef, FerriteDistributed.num_local_dofs(dh))
FerriteDistributed.extract_local_part!(u_local, uₕ, dh)

vtk_grid("heat_equation_distributed", dh) do vtk
    vtk_point_data(vtk, dh, u_local)

    vtk_shared_vertices(vtk, dgrid)
    vtk_shared_faces(vtk, dgrid)
    vtk_partitioning(vtk, dgrid)
end

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
