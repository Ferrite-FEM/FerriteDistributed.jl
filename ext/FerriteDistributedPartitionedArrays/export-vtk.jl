function Ferrite.write_solution(vtk::FerriteDistributed.PVTKGridFile, dh::FerriteDistributed.NODDofHandler, u::PVector, suffix="")
    inv_perm = invperm(_nod_to_oag_perm(dh))
    map(local_values(u)) do u_local
        Ferrite.write_solution(vtk, dh, u_local[inv_perm], suffix)
    end
end
