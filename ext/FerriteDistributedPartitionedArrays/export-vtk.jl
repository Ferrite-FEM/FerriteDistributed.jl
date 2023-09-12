function WriteVTK.vtk_point_data(vtk, dh::Ferrite.AbstractDofHandler, u::PVector)
    map(local_values(u)) do u_local
        vtk_point_data(vtk, dh, u_local)
    end
end
