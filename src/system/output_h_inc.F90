!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch.
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!! $Id: output_h_inc.F90 15122 2016-02-25 23:32:21Z xavier $

  ! ---------------------------------------------------------
  subroutine output_hamiltonian(hm, st, der, dir, outp, geo, grp)
    type(hamiltonian_t),       intent(in)    :: hm
    type(states_t),            intent(inout) :: st
    type(derivatives_t),       intent(inout) :: der
    character(len=*),          intent(in)    :: dir
    type(output_t),            intent(in)    :: outp
    type(geometry_t),          intent(in)    :: geo
    type(mpi_grp_t), optional, intent(in)    :: grp !< the group that shares the same data, must contain the domains group

    integer :: is, err, idir, ispin
    character(len=MAX_PATH_LEN) :: fname
    type(base_potential_iterator_t)        :: iter
    type(base_potential_t),        pointer :: subsys_external
    type(base_hamiltonian_t),      pointer :: subsys_tnadd
    character(len=BASE_POTENTIAL_NAME_LEN) :: name
    FLOAT,         dimension(:),   pointer :: xpot
    FLOAT,         dimension(:,:), pointer :: tnadd_potential
    FLOAT, allocatable :: v0(:,:), nxc(:), potential(:)
    FLOAT, allocatable :: current(:, :, :)

    PUSH_SUB(output_hamiltonian)

    if(iand(outp%what, OPTION__OUTPUT__POTENTIAL) /= 0) then
      if(hm%cmplxscl%space) then
        call zio_function_output(outp%how, dir, "v0", der%mesh,&
          hm%ep%vpsl + M_zI*hm%ep%Imvpsl, units_out%energy, err, geo = geo, grp = grp)
      else  
        SAFE_ALLOCATE(v0(1:der%mesh%np, 1:hm%d%dim))
        v0(1:der%mesh%np, 1) = hm%ep%vpsl(1:der%mesh%np)
        call dio_function_output(outp%how, dir, "v0", der%mesh, v0(:, 1), units_out%energy, err, geo = geo, grp = grp)
        SAFE_DEALLOCATE_A(v0)
      end if

      if(hm%ep%classical_pot > 0) then
        call dio_function_output(outp%how, dir, "vc", der%mesh, hm%ep%Vclassical, units_out%energy, err, geo = geo, grp = grp)
      end if

      nullify(subsys_external, xpot)
      if(associated(hm%ep%subsys_external))then
        call base_potential_init(iter, hm%ep%subsys_external)
        do
          nullify(subsys_external, xpot)
          call base_potential_next(iter, name, subsys_external, err)
          if(err/=BASE_POTENTIAL_OK)exit
          ASSERT(associated(subsys_external))
          call base_potential_get(subsys_external, xpot)
          ASSERT(associated(xpot))
          write(fname, "(a,'-',a)") "v0", trim(adjustl(name))
          call dio_function_output(outp%how, dir, fname, der%mesh, &
            xpot, units_out%energy, err, geo = geo, grp = grp)
        end do
        call base_potential_end(iter)
        nullify(subsys_external, xpot)
      end if

      if(hm%theory_level /= INDEPENDENT_PARTICLES) then
        if (.not. hm%cmplxscl%space) then 
          call dio_function_output(outp%how, dir, 'vh', der%mesh, hm%vhartree, units_out%energy, err, geo = geo, grp = grp)
        else
          call zio_function_output(outp%how, dir, 'vh', der%mesh, &
            hm%vhartree(1:der%mesh%np) + M_zI*hm%Imvhartree(1:der%mesh%np), units_out%energy, err, geo = geo, grp = grp)
        end if

        nullify(subsys_tnadd, tnadd_potential)
        if(associated(hm%subsys_hm))then
          call base_hamiltonian_get(hm%subsys_hm, "tnadd", subsys_tnadd)
          ASSERT(associated(subsys_tnadd))
          call base_hamiltonian_get(subsys_tnadd, nspin=ispin)
          call base_hamiltonian_get(subsys_tnadd, tnadd_potential)
          ASSERT(associated(tnadd_potential))
          nullify(subsys_tnadd)
          do is = 1, min(ispin, 2)
            if(ispin == 1) then
              write(fname, '(a)') 'tnadd'
            else
              write(fname, '(a,i1)') 'tnadd-sp', is
            end if
            call dio_function_output(outp%how, dir, fname, der%mesh, &
              tnadd_potential(:,is), units_out%energy, err, geo = geo, grp = grp)
          end do
        end if
        
        SAFE_ALLOCATE(potential(1:der%mesh%np))
        do is = 1, min(hm%d%ispin, 2)
          if(hm%d%ispin == 1) then
            write(fname, '(a)') 'vxc'
          else
            write(fname, '(a,i1)') 'vxc-sp', is
          end if
          if(.not. hm%cmplxscl%space) then
            call dio_function_output(outp%how, dir, fname, der%mesh, hm%vxc(:, is), units_out%energy, err, geo = geo, grp = grp)
          else
            call zio_function_output(outp%how, dir, fname, der%mesh, &
              hm%vxc(:, is) + M_zI *  hm%Imvxc(:, is), units_out%energy, err, geo = geo, grp = grp)
          end if
          
          ! finally the full KS potential (without non-local PP contributions)
          if(associated(tnadd_potential))then
            potential = hm%ep%vpsl + hm%vhxc(:, is) - tnadd_potential(:, min(is,ispin))
          else
            potential = hm%ep%vpsl + hm%vhxc(:, is)
          end if
          if(hm%d%ispin == 1) then
            write(fname, '(a)') 'vks'
          else
            write(fname, '(a,i1)') 'vks-sp', is
          end if
          if (hm%ep%classical_pot > 0) then
            call dio_function_output(outp%how, dir, fname, der%mesh, &
              potential + hm%ep%Vclassical, units_out%energy, err, geo = geo, grp = grp)
          else
            if(.not. hm%cmplxscl%space) then
              call dio_function_output(outp%how, dir, fname, der%mesh, &
                potential, units_out%energy, err, geo = geo, grp = grp)
            else
              call zio_function_output(outp%how, dir, fname, der%mesh, &
                potential + M_zI * hm%ep%Imvpsl + M_zI * hm%Imvhxc(:, is), units_out%energy, &
                err, geo = geo, grp = grp)
            end if
          end if
        end do
        SAFE_DEALLOCATE_A(potential)
        nullify(tnadd_potential)
      end if

      if(hm%self_induced_magnetic) then
        ! unit of magnetic field is same as of electric field, and same as force (since e = 1)
        select case(der%mesh%sb%dim)
        case(3)
          do idir = 1, der%mesh%sb%dim
            call dio_function_output(outp%how, dir, 'Bind_'//index2axis(idir), der%mesh, hm%b_ind(:, idir), &
              units_out%force, err, geo = geo, grp = grp)
          end do
        case(2)
          call dio_function_output(outp%how, dir, 'Bind_z', der%mesh, hm%b_ind(:, 1), units_out%force, err, geo = geo, grp = grp)
        end select
      end if
    end if

    if(iand(outp%what, OPTION__OUTPUT__XC_DENSITY) /= 0 .and. hm%theory_level /= INDEPENDENT_PARTICLES) then
      SAFE_ALLOCATE(v0(1:der%mesh%np_part, 1))
      SAFE_ALLOCATE(nxc(1:der%mesh%np))

      do is = 1, min(hm%d%ispin, 2)
        if(hm%d%ispin == 1) then
          write(fname, '(a)') 'nxc'
        else
          write(fname, '(a,i1)') 'nxc-sp', is
        end if
                
        v0(1:der%mesh%np, 1) = hm%vxc(1:der%mesh%np, is)

        call dderivatives_lapl(der, v0(:, 1), nxc)

        call dio_function_output(outp%how, dir, fname, der%mesh, nxc, units_out%energy, err, geo = geo, grp = grp)
        
      end do

      SAFE_DEALLOCATE_A(v0)
      SAFE_DEALLOCATE_A(nxc)
    end if

    if(iand(outp%what, OPTION__OUTPUT__CURRENT) /= 0) then
      if(states_are_complex(st)) then
        ! calculate current first
        SAFE_ALLOCATE(current(1:der%mesh%np_part, 1:der%mesh%sb%dim, 1:hm%d%nspin))
        call current_calculate(outp%current_calculator, der, hm, geo, st, current)
        do is = 1, hm%d%nspin

          if(st%d%nspin == 1) then
            write(fname, '(2a)') 'current'
          else
            write(fname, '(a,i1)') 'current-sp', is
          end if
          
          call io_function_output_vector(outp%how, dir, fname, der%mesh, &
            current(:, :, is), der%mesh%sb%dim, (unit_one/units_out%time)*units_out%length**(1 - der%mesh%sb%dim), err, &
            geo = geo, grp = st%dom_st_kpt_mpi_grp, vector_dim_labels = (/'x', 'y', 'z'/))

        end do
        SAFE_DEALLOCATE_A(current)
      else
        message(1) = 'No current density output for real states since it is identically zero.'
        call messages_warning(1)
      end if
    end if
    
    POP_SUB(output_hamiltonian)
  end subroutine output_hamiltonian


  ! ---------------------------------------------------------
  subroutine output_scalar_pot(outp, gr, geo, hm, dir, time)
    type(grid_t),         intent(inout) :: gr
    type(geometry_t),     intent(in)    :: geo
    type(hamiltonian_t),  intent(inout) :: hm
    type(output_t),       intent(in)    :: outp
    character(len=*),     intent(in)    :: dir
    FLOAT, optional,      intent(in)    :: time

    integer :: is, err
    character(len=80) :: fname
    FLOAT, allocatable :: scalar_pot(:)

    PUSH_SUB(output_scalar_pot)

    if(iand(outp%what, OPTION__OUTPUT__EXTERNAL_TD_POTENTIAL) /= 0) then
      SAFE_ALLOCATE(scalar_pot(1:gr%mesh%np))
      do is = 1, hm%ep%no_lasers
        write(fname, '(a,i1)') 'scalar_pot-', is
        scalar_pot = M_ZERO
        call laser_potential(hm%ep%lasers(is), gr%mesh, scalar_pot, time=time)
        call dio_function_output(outp%how, dir, fname, gr%mesh, scalar_pot, units_out%energy, err, geo = geo)
      end do
      SAFE_DEALLOCATE_A(scalar_pot)
    end if

    POP_SUB(output_scalar_pot)
  end subroutine output_scalar_pot


  ! ---------------------------------------------------------
  subroutine output_kick(outp, mesh, geo, kick, dir)
    type(mesh_t),     intent(in) :: mesh
    type(geometry_t), intent(in) :: geo
    type(kick_t),     intent(in) :: kick
    type(output_t),   intent(in) :: outp
    character(len=*), intent(in) :: dir

    integer :: err
    CMPLX, allocatable :: kick_function(:)
    
    PUSH_SUB(output_kick)

    if(iand(outp%what, OPTION__OUTPUT__DELTA_PERTURBATION) /= 0) then
      SAFE_ALLOCATE(kick_function(1:mesh%np))
      call kick_function_get(mesh, kick, kick_function)
      call zio_function_output(outp%how, dir, "kick_function", mesh, kick_function(:), &
        units_out%energy, err, geo = geo)
      SAFE_DEALLOCATE_A(kick_function)
    end if
    
    POP_SUB(output_kick)
  end subroutine output_kick


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
