using FerriteDistributed
using Test

import FerriteDistributed: CoverTopology, global_rank

MPI.Init()
@testset "MPI setup 2" begin
    @test MPI.Comm_size(MPI.COMM_WORLD) == 2
end

@testset "distributed grid construction 2" begin
    # We do not cover subcommunicators for now.
    comm = MPI.COMM_WORLD

    global_grid = generate_grid(Hexahedron, (2, 1, 1))
    global_topology = CoverTopology(global_grid)
    dgrid = NODGrid(comm, global_grid, global_topology, [2, 1])
    my_rank = global_rank(dgrid)
    if my_rank == 1
        # Edges
        @test length(get_shared_edges(dgrid)) == 4
        function check_edge_correctly_shared_1(idx_local, idx_nonlocal)
            se = get_shared_edge(dgrid, idx_local)
            @test FerriteDistributed.remote_entities(se) == Dict(2 => [idx_nonlocal])
        end
        check_edge_correctly_shared_1(EdgeIndex(1,4), EdgeIndex(1,2))
        check_edge_correctly_shared_1(EdgeIndex(1,12), EdgeIndex(1,11))
        check_edge_correctly_shared_1(EdgeIndex(1,9), EdgeIndex(1,10))
        check_edge_correctly_shared_1(EdgeIndex(1,8), EdgeIndex(1,6))

        # Faces
        @test length(get_shared_faces(dgrid)) == 1
        sf = get_shared_face(dgrid, FaceIndex(1,5))
        @test FerriteDistributed.remote_entities(sf) == Dict(2 => [FaceIndex(1,3)])
    elseif my_rank == 2
        # Edges
        @test length(get_shared_edges(dgrid)) == 4
        function check_edge_correctly_shared_2(idx_nonlocal, idx_local)
            se = get_shared_edge(dgrid, idx_local)
            @test FerriteDistributed.remote_entities(se) == Dict(1 => [idx_nonlocal])
        end
        check_edge_correctly_shared_2(EdgeIndex(1,4), EdgeIndex(1,2))
        check_edge_correctly_shared_2(EdgeIndex(1,12), EdgeIndex(1,11))
        check_edge_correctly_shared_2(EdgeIndex(1,9), EdgeIndex(1,10))
        check_edge_correctly_shared_2(EdgeIndex(1,8), EdgeIndex(1,6))

        # Faces
        @test length(get_shared_faces(dgrid)) == 1
        sf = get_shared_face(dgrid, FaceIndex(1,3))
        @test FerriteDistributed.remote_entities(sf) == Dict(1 => [FaceIndex(1,5)])
    end
    MPI.Finalize()
end
