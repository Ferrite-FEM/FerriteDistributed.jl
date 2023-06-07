function WriteVTK.vtk_point_data(vtk, dh::Ferrite.AbstractDofHandler, u::PVector)
    map_parts(local_view(u, u.rows)) do u_local
        vtk_point_data(vtk, dh, u_local)
    end
end
