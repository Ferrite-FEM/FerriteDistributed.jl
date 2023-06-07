using Ferrite, FerriteDistributed, MPI
using Test

MPI.Init()

@testset "MPI setup 1" begin
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

if VERSION >= v"1.10.0-DEV.90"
    using Metis
    @testset "Metis integration" begin
        grid = generate_grid(Quadrilateral, (4,4))
        dgrid = NODGrid(MPI.COMM_WORLD, grid, FerriteDistributed.PartitioningAlgorithm.Ext{Metis}(:RECURSIVE))

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
end # VERSION >= v"1.10.0-DEV.90"

# @testset "Dof stuff" begin
#     # Consistency check for dof computation.
#     grid = generate_grid(Hexahedron, (2, 2, 2))
#     dh = DofHandler(grid)
#     add!(dh, :u, Lagrange{RefHexahedron,2}()^3)
#     add!(dh, :v, Lagrange{RefHexahedron,2}())
#     add!(dh, :w, Lagrange{RefHexahedron,1}()^3)
#     add!(dh, :x, Lagrange{RefHexahedron,2}()^3)
#     _, vertexdicts, edgedicts, facedicts = Ferrite.__close!(dh)
#     @test Ferrite.find_field(dh, :u) == (1,1)
#     @test Ferrite.find_field(dh, :v) == (1,2)
#     @test Ferrite.find_field(dh, :w) == (1,3)
#     for (ci,cell) ∈ enumerate(getcells(grid))
#         for field_idx ∈ 1:3
#             for vi ∈ 1:length(Ferrite.vertices(cell))
#                 vertex = VertexIndex(ci,vi)
#                 global_vertex = Ferrite.toglobal(grid, vertex)
#                 dofs = vertex_dofs(dh, field_idx, vertex)
#                 if !haskey(vertexdicts[field_idx], global_vertex)
#                     @test empty(dofs)
#                     @test has_vertex_dofs(dh, field_idx, vertex)
#                 else 
#                     @test vertexdicts[field_idx][global_vertex] == dofs[1]
#                     @test has_vertex_dofs(dh, field_idx, vertex)
#                 end
#             end

#             for fi ∈ 1:length(Ferrite.faces(cell))
#                 face = FaceIndex(ci,fi)
#                 global_face = Ferrite.toglobal(grid, face)
#                 dofs = face_dofs(dh, field_idx, face)
#                 if !haskey(facedicts[field_idx], global_face)
#                     @test isempty(dofs)
#                     @test has_face_dofs(dh, field_idx, face)
#                 else
#                     @test facedicts[field_idx][global_face] == dofs[1]
#                     @test has_face_dofs(dh, field_idx, face)
#                 end
#             end

#             for ei ∈ 1:length(Ferrite.edges(cell))
#                 edge = EdgeIndex(ci,ei)
#                 global_edge = Ferrite.toglobal(grid, edge)
#                 dofs = edge_dofs(dh, field_idx, edge)
#                 if !haskey(edgedicts[field_idx], global_edge) 
#                     @test isempty(dofs) 
#                     @test has_edge_dofs(dh, field_idx, edge)
#                 else
#                     @test edgedicts[field_idx][global_edge][1] == dofs[1]
#                     @test has_edge_dofs(dh, field_idx, edge)
#                 end
#             end
#         end
#     end
# end

@testset "FerriteMPI n=2" begin
    n = 2  # number of processes
    mpiexec() do exe  # MPI wrapper
        run(`$exe -n $n $(Base.julia_cmd()) test_distributed_impl_2.jl`)
    end
end

@testset "FerriteMPI n=3" begin
    n = 3  # number of processes
    mpiexec() do exe  # MPI wrapper
        run(`$exe -n $n $(Base.julia_cmd()) test_distributed_impl_3.jl`)
    end
end

@testset "FerriteMPI n=5" begin
    n = 5  # number of processes
    mpiexec() do exe  # MPI wrapper
        run(`$exe -n $n $(Base.julia_cmd()) test_distributed_impl_5.jl`)
    end
end

