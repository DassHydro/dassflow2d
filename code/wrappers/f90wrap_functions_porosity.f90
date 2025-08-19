! Module fonctions_porosite_mod defined in file functions_porosity.f90

subroutine f90wrap_fonctions_porosite_mod__update_all_porosities(dof, mesh)
    use m_model, only: unk
    use fonctions_porosite_mod, only: update_all_porosities
    use m_mesh, only: msh
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    type(unk_ptr_type) :: dof_ptr
    integer, intent(in), dimension(2) :: dof
    type(msh_ptr_type) :: mesh_ptr
    integer, intent(in), dimension(2) :: mesh
    dof_ptr = transfer(dof, dof_ptr)
    mesh_ptr = transfer(mesh, mesh_ptr)
    call update_all_porosities(dof=dof_ptr%p, mesh=mesh_ptr%p)
end subroutine f90wrap_fonctions_porosite_mod__update_all_porosities

! End of module fonctions_porosite_mod defined in file functions_porosity.f90

