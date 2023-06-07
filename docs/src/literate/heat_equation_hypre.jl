# # Distributed Assembly with HYPRE.jl
#
# ## Introduction
#
# Now we want to solve the heat problem in parallel. To be specific, this example shows
# how to utilize process parallelism to assemble finite element matrices in parallel.
# This example presumes that the reader is familiar with solving the heat problem in
# serial with Ferrite.jl, as presented in [the first example](@ref heat_example).
#
#-
# ## Commented Program
#
# Now we solve the problem in Ferrite. What follows is a program spliced with comments.
#md # The full program, without comments, can be found in the next [section](@ref heat_equation-plain-program).
#
# First we load Ferrite, and some other packages we need
using FerriteDistributed
using HYPRE, Metis

import FerriteDistributed: getglobalgrid, global_comm, local_dof_range #TODO REMOVE THIS

# Launch MPI and HYPRE
MPI.Init()
HYPRE.Init()

# We start generating a simple grid with 10x10x10 hexahedral elements
# and distribute it across our processors using `generate_distributed_grid`.
dgrid = generate_nod_grid(MPI.COMM_WORLD, Hexahedron, (10, 10, 10); partitioning_alg=FerriteDistributed.PartitioningAlgorithm.Metis(:RECURSIVE));

# ### Trial and test functions
# Nothing changes here.
ref = RefHexahedron
ip = Lagrange{ref, 2}()
ip_geo = Lagrange{ref, 1}()
qr = QuadratureRule{ref}(3)
cellvalues = CellValues(qr, ip, ip_geo);

# ### Degrees of freedom
# Nothing changes here, too. The constructor takes care of creating the correct distributed dof handler.
dh = DofHandler(dgrid)
push!(dh, :u, 1, ip)
close!(dh);

# ### Boundary conditions
# Nothing has to be changed here either.
ch = ConstraintHandler(dh);
∂Ω = union(getfaceset.((dgrid, ), ["left", "right", "top", "bottom", "front", "back"])...);
dbc = Dirichlet(:u, ∂Ω, (x, t) -> 0)
dbc_val = 0                                 #src
dbc = Dirichlet(:u, ∂Ω, (x, t) -> dbc_val)  #src
add!(ch, dbc);
close!(ch)

# ### Assembling the linear system
# Assembling the system works also mostly analogue.
function doassemble(cellvalues::CellValues, dh::FerriteDistributed.NODDofHandler{dim}, ch::ConstraintHandler) where {dim}
    n_basefuncs = getnbasefunctions(cellvalues)
    Ke = zeros(n_basefuncs, n_basefuncs)
    fe = zeros(n_basefuncs)

    # --------------------- Distributed assembly --------------------
    # The synchronization with the global sparse matrix is handled by
    # an assembler again. You can choose from different backends, which
    # are described in the docs and will be expaned over time. This call
    # may trigger a large amount of communication.

    # TODO how to put this into an interface?
    dgrid = getglobalgrid(dh)
    comm = global_comm(dgrid)
    ldofrange = local_dof_range(dh)
    K = HYPREMatrix(comm, first(ldofrange), last(ldofrange))
    f = HYPREVector(comm, first(ldofrange), last(ldofrange))

    assembler = start_assemble(K, f)

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

        apply_local!(Ke, fe, celldofs(cell), ch)

        # TODO how to put this into an interface.
        assemble!(assembler, dh.ldof_to_gdof[celldofs(cell)], fe, Ke)
    end

    # Finally, for the `HYPREAssembler` we have to call
    # `end_assemble` to construct the global sparse matrix and the global
    # right hand side vector.
    end_assemble(assembler)

    return K, f
end
#md nothing # hide

# ### Solution of the system
# Again, we assemble our problem. Note that we applied the constraints locally.
K, f = doassemble(cellvalues, dh, ch);

# We use CG with AMG preconditioner to solve the system.
precond = HYPRE.BoomerAMG()
solver = HYPRE.PCG(; Precond = precond)
uₕ = HYPRE.solve(solver, K, f)

# And convert the solution from HYPRE to Ferrite
u_local = Vector{Float64}(undef, FerriteDistributed.num_local_dofs(dh))
FerriteDistributed.extract_local_part!(u_local, uₕ, dh)

# ### Exporting via PVTK
# To visualize the result we export the grid and our field `u`
# to a VTK-file, which can be viewed in e.g. [ParaView](https://www.paraview.org/).
vtk_grid("heat_equation_distributed", dh) do vtk
    vtk_point_data(vtk, dh, u_local)
    # For debugging purposes it can be helpful to enrich
    # the visualization with some meta  information about
    # the grid and its partitioning
    vtk_shared_vertices(vtk, dgrid)
    vtk_shared_faces(vtk, dgrid)
    vtk_shared_edges(vtk, dgrid) #src
    vtk_partitioning(vtk, dgrid)
end

## Test the result against the manufactured solution                    #src
using Test                                                              #src
for cell in CellIterator(dh)                                            #src
    reinit!(cellvalues, cell)                                           #src
    n_basefuncs = getnbasefunctions(cellvalues)                         #src
    coords = getcoordinates(cell)                                       #src
    uₑ = u_local[celldofs(cell)]                                        #src
    for q_point in 1:getnquadpoints(cellvalues)                         #src
        x = spatial_coordinate(cellvalues, q_point, coords)             #src
        for i in 1:n_basefuncs                                          #src
            uₐₙₐ    = prod(cos, x*π/2)+dbc_val                          #src
            uₐₚₚᵣₒₓ = function_value(cellvalues, q_point, uₑ)           #src
            @test isapprox(uₐₙₐ, uₐₚₚᵣₒₓ; atol=1e-1)                    #src
        end                                                             #src
    end                                                                 #src
end                                                                     #src

# Finally, we gracefully shutdown MPI
MPI.Finalize()

#md # ## [Plain program](@id distributed-assembly-plain-hypre)
#md #
#md # Here follows a version of the program without any comments.
#md # The file is also available here: [`distributed_assembly_hypre.jl`](distributed_assembly_hypre.jl).
#md #
#md # ```julia
#md # @__CODE__
#md # ```
