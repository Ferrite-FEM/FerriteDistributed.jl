using Ferrite, FerriteDistributed, MPI
using Test

MPI.Init()

@testset "MPI Setup" begin
    @test MPI.Comm_rank(MPI.COMM_WORLD) == 0
    @test MPI.Comm_size(MPI.COMM_WORLD) == 1
end

@testset "Grid construction nprocs=1" begin
    grid = generate_grid(Quadrilateral, (4,4))
    partitioning = create_partitioning(grid, FerriteDistributed.CoverTopology(grid), 4, FerriteDistributed.PartitioningAlgorithm.SFC())
    @test all(partitioning .== [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4])

    dgrid = NODGrid(MPI.COMM_WORLD, grid)
    @test isempty(get_shared_vertices(dgrid))
    @test isempty(get_shared_edges(dgrid))
    @test isempty(get_shared_faces(dgrid))
    @test getvertexsets(getlocalgrid(dgrid)) == getvertexsets(grid)
    @test getfacesets(getlocalgrid(dgrid)) == getfacesets(grid)
    @test getcellsets(getlocalgrid(dgrid)) == getcellsets(grid)
    @test getnnodes(getlocalgrid(dgrid) ) == getnnodes(grid)
    @test getncells(getlocalgrid(dgrid) ) == getncells(grid)

    # Tests 1-based indexing!
    @test FerriteDistributed.global_rank(dgrid) == 1

    @test FerriteDistributed.global_nranks(dgrid) == 1
end

using Metis
@testset "Metis integration" begin
    grid = generate_grid(Quadrilateral, (4,4))
    dgrid = NODGrid(MPI.COMM_WORLD, grid, FerriteDistributed.PartitioningAlgorithm.Metis(:RECURSIVE))

    @test isempty(get_shared_vertices(dgrid))
    @test isempty(get_shared_edges(dgrid))
    @test isempty(get_shared_faces(dgrid))
    @test getvertexsets(getlocalgrid(dgrid)) == getvertexsets(grid)
    @test getfacesets(getlocalgrid(dgrid)) == getfacesets(grid)
    @test getcellsets(getlocalgrid(dgrid)) == getcellsets(grid)
    @test getnnodes(getlocalgrid(dgrid) ) == getnnodes(grid)
    @test getncells(getlocalgrid(dgrid) ) == getncells(grid)

    @test FerriteDistributed.global_rank(dgrid) == 1

    @test FerriteDistributed.global_nranks(dgrid) == 1
end

#TODO MPI based test sets with more than 1 proc
