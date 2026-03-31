function Ferrite.write_solution(vtk::FerriteDistributed.PVTKGridFile, dh::FerriteDistributed.NODDofHandler, u::PVector, suffix="")
    u_local = FerriteDistributed.gather_dof_values(u, dh)
    Ferrite.write_solution(vtk, dh, u_local, suffix)
end
