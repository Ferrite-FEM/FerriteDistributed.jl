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
    @test getlocalgrid(dgrid).vertexsets == grid.vertexsets
    @test getlocalgrid(dgrid).facetsets == grid.facetsets
    @test getlocalgrid(dgrid).cellsets == grid.cellsets
    @test getnnodes(getlocalgrid(dgrid) ) == getnnodes(grid)
    @test getncells(getlocalgrid(dgrid) ) == getncells(grid)

    @test FerriteDistributed.global_rank(dgrid) == 1

    @test FerriteDistributed.global_nranks(dgrid) == 1
end

@testset "CoverTopology 2D" begin
    grid = generate_grid(Quadrilateral, (2,2))
    top = CoverTopology(grid)

    # cell_neighbor: all 4 cells are connected (vertex-connected counts)
    for i in 1:4
        cn = Ferrite.getneighborhood(top, grid, CellIndex(i))
        @test sort(collect(cn)) == setdiff(1:4, i)
    end

    # vertex_vertex_neighbor: center vertex (node 5) shared by all 4 cells
    # Cell 1 local vertex 3 = global node 5
    vn = Ferrite.getneighborhood(top, grid, VertexIndex(1, 3))
    @test length(vn) == 3
    @test Set(v[1] for v in vn) == Set([2, 3, 4])

    # Cell 1 local vertex 2 = global node 2, shared only with cell 2
    vn12 = Ferrite.getneighborhood(top, grid, VertexIndex(1, 2))
    @test length(vn12) == 1
    @test vn12[1] == VertexIndex(2, 1)

    # include_self
    vn_self = Ferrite.getneighborhood(top, grid, VertexIndex(1, 2), true)
    @test length(vn_self) == 2
    @test VertexIndex(1, 2) ∈ vn_self

    # face_face_neighbor: empty in 2D (faces = 2D cell surface, never shared)
    for ci in 1:4
        @test isempty(Ferrite.getneighborhood(top, grid, FaceIndex(ci, 1)))
    end
end

@testset "CoverTopology 3D" begin
    grid = generate_grid(Hexahedron, (2,1,1))
    top = CoverTopology(grid)

    # cell_neighbor: 2 cells, each other's neighbor
    @test collect(Ferrite.getneighborhood(top, grid, CellIndex(1))) == [2]
    @test collect(Ferrite.getneighborhood(top, grid, CellIndex(2))) == [1]

    # face_face_neighbor: cells share exactly 1 face pair
    face_pairs = Tuple{Int,Int}[]
    for ci in 1:2, lfi in 1:6
        for fn in Ferrite.getneighborhood(top, grid, FaceIndex(ci, lfi))
            push!(face_pairs, (ci, fn[1]))
        end
    end
    @test length(face_pairs) == 2
    @test (1, 2) ∈ face_pairs
    @test (2, 1) ∈ face_pairs

    # edge_edge_neighbor: cells share 4 edges (the face boundary)
    edge_pairs = Tuple{Int,Int}[]
    for ci in 1:2, lei in 1:12
        for en in Ferrite.getneighborhood(top, grid, EdgeIndex(ci, lei))
            push!(edge_pairs, (ci, en[1]))
        end
    end
    @test count(p -> p == (1, 2), edge_pairs) == 4
    @test count(p -> p == (2, 1), edge_pairs) == 4

    # vertex_vertex_neighbor: cells share 4 vertices
    vert_pairs = Tuple{Int,Int}[]
    for ci in 1:2, lvi in 1:8
        for vn in Ferrite.getneighborhood(top, grid, VertexIndex(ci, lvi))
            push!(vert_pairs, (ci, vn[1]))
        end
    end
    @test count(p -> p == (1, 2), vert_pairs) == 4
    @test count(p -> p == (2, 1), vert_pairs) == 4

    # Cover semantics: face-sharing cells also appear in edge AND vertex neighbors
    # (Unlike ExclusiveTopology which only stores the highest-dimensional connection)
    has_face = !isempty(face_pairs)
    has_edge = count(p -> p == (1, 2), edge_pairs) > 0
    has_vert = count(p -> p == (1, 2), vert_pairs) > 0
    @test has_face && has_edge && has_vert  # all three for face-sharing cells
end

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

