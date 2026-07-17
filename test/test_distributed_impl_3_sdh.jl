using FerriteDistributed
using Test

import FerriteDistributed: global_rank, num_global_dofs

MPI.Init()
comm = MPI.COMM_WORLD
@testset "MPI setup 3 sdh" begin
    @test MPI.Comm_size(comm) == 3
end

# Check that the distributed dof numbering is consistent with a serial reference:
# - the global number of dofs matches
# - every global dof is owned by exactly one rank
# - matching (cell, dof position) pairs between the serial and the distributed handler
#   induce a bijection between serial dofs and global dofs
function check_consistency(dh, serial_dh, partitioning)
    ngdofs = num_global_dofs(dh)
    @test ngdofs == ndofs(serial_dh)

    my_rank = global_rank(FerriteDistributed.getglobalgrid(dh))
    owned = zeros(Int, ngdofs)
    for (ldof, gdof) in enumerate(dh.ldof_to_gdof)
        @test 1 <= gdof <= ngdofs
        if dh.ldof_to_rank[ldof] == my_rank
            owned[gdof] += 1
        end
    end
    MPI.Allreduce!(owned, +, comm)
    @test all(==(1), owned)

    # (serial dof, global dof) pairs for each local cell; local cells appear in ascending
    # global cell order during distribution
    global_cells = findall(==(my_rank), partitioning)
    dof_pairs = Int[]
    for (ci, gci) in enumerate(global_cells)
        Ferrite.ndofs_per_cell(dh, ci) == 0 && continue
        serial_dofs = celldofs(serial_dh, gci)
        global_dofs = dh.ldof_to_gdof[celldofs(dh, ci)]
        @test length(serial_dofs) == length(global_dofs)
        for (s, g) in zip(serial_dofs, global_dofs)
            push!(dof_pairs, s, g)
        end
    end
    lengths = MPI.Allgather(Int32(length(dof_pairs)), comm)
    all_pairs = MPI.VBuffer(Vector{Int}(undef, sum(lengths)), lengths)
    MPI.Allgatherv!(dof_pairs, all_pairs, comm)
    serial_to_global = Dict{Int,Int}()
    global_to_serial = Dict{Int,Int}()
    consistent = true
    for i in 1:2:length(all_pairs.data)
        s, g = all_pairs.data[i], all_pairs.data[i + 1]
        consistent &= get!(serial_to_global, s, g) == g
        consistent &= get!(global_to_serial, g, s) == s
    end
    @test consistent
    @test length(serial_to_global) == ndofs(serial_dh)
end

# Higher order interpolations: interior edge dofs must be synchronized across the
# process boundary. The cubic case has two interior dofs per edge whose order depends
# on the edge orientation.
@testset "distributed dofs 2D higher order" begin
    for ip in (
            Lagrange{RefQuadrilateral, 2}(),
            Lagrange{RefQuadrilateral, 3}(),
            Lagrange{RefQuadrilateral, 3}()^2,
        )
        grid = generate_grid(Quadrilateral, (3, 3))
        partitioning = [1, 1, 2, 2, 3, 3, 1, 2, 3]
        dgrid = NODGrid(comm, grid, CoverTopology(grid), partitioning)

        dh = DofHandler(dgrid)
        add!(dh, :u, ip)
        close!(dh)

        serial_dh = DofHandler(grid)
        add!(serial_dh, :u, ip)
        close!(serial_dh)

        check_consistency(dh, serial_dh, partitioning)
    end
end

# Multiple SubDofHandlers: field :u everywhere, field :p only on the left half. The
# subdomain of :p touches the process boundary between ranks 2 and 3, and rank 3 does
# not know the field :p at all.
@testset "distributed dofs multiple SubDofHandlers" begin
    grid = generate_grid(Quadrilateral, (4, 2))
    addcellset!(grid, "left", Set([1, 2, 5, 6]))
    addcellset!(grid, "right", Set([3, 4, 7, 8]))
    partitioning = [1, 2, 2, 3, 1, 2, 3, 3]
    dgrid = NODGrid(comm, grid, CoverTopology(grid), partitioning)
    my_rank = global_rank(dgrid)

    ip = Lagrange{RefQuadrilateral, 1}()
    dh = DofHandler(dgrid)
    left = getcellset(dgrid, "left")
    right = getcellset(dgrid, "right")
    if !isempty(left)
        sdh_left = SubDofHandler(dh, left)
        add!(sdh_left, :u, ip)
        add!(sdh_left, :p, ip)
    end
    if !isempty(right)
        sdh_right = SubDofHandler(dh, right)
        add!(sdh_right, :u, ip)
    end
    close!(dh)

    if my_rank == 3
        @test Ferrite.getfieldnames(dh) == [:u]
    else
        @test Ferrite.getfieldnames(dh) == [:u, :p]
    end

    serial_dh = DofHandler(grid)
    serial_left = SubDofHandler(serial_dh, getcellset(grid, "left"))
    add!(serial_left, :u, ip)
    add!(serial_left, :p, ip)
    serial_right = SubDofHandler(serial_dh, getcellset(grid, "right"))
    add!(serial_right, :u, ip)
    close!(serial_dh)

    check_consistency(dh, serial_dh, partitioning)

    # ConstraintHandler on the wrapper uses the SubDofHandler-aware Ferrite machinery
    ch = ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, getfacetset(dgrid, "left"), x -> 0.0))
    close!(ch)
    nconstrained = MPI.Allreduce(length(ch.prescribed_dofs), MPI.SUM, comm)
    @test nconstrained == 3 # the left boundary lives on rank 1 only and has 3 nodes
end

MPI.Finalize()
