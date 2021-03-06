!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! $Id: forces_inc.F90 15590 2016-08-23 13:23:38Z ssato $

subroutine X(forces_gather)(geo, force)
  type(geometry_t), intent(in)    :: geo
  R_TYPE,           intent(inout) :: force(:, :)
  
  R_TYPE,  allocatable :: force_local(:, :)
  integer, allocatable :: recv_count(:), recv_displ(:)

  PUSH_SUB(X(forces_gather))

  call profiling_in(prof_comm, "FORCES_COMM")
    
  ! each node has a piece of the force array, they have to be
  ! collected by all nodes
  
  ASSERT(ubound(force, dim = 1) == geo%space%dim)

  SAFE_ALLOCATE(recv_count(1:geo%atoms_dist%mpi_grp%size))
  SAFE_ALLOCATE(recv_displ(1:geo%atoms_dist%mpi_grp%size))
  SAFE_ALLOCATE(force_local(1:geo%space%dim, 1:max(1, geo%atoms_dist%nlocal)))
  
  recv_count(1:geo%atoms_dist%mpi_grp%size) = geo%space%dim*geo%atoms_dist%num(0:geo%atoms_dist%mpi_grp%size - 1)
  recv_displ(1:geo%atoms_dist%mpi_grp%size) = geo%space%dim*(geo%atoms_dist%range(1, 0:geo%atoms_dist%mpi_grp%size - 1) - 1)
  
  if(geo%atoms_dist%nlocal > 0) then
    force_local(1:geo%space%dim, 1:geo%atoms_dist%nlocal) = force(1:geo%space%dim, geo%atoms_dist%start:geo%atoms_dist%end)
  end if

#ifdef HAVE_MPI  
  call MPI_Allgatherv(&
    force_local(1, 1), geo%space%dim*geo%atoms_dist%nlocal, R_MPITYPE, &
    force(1, 1), recv_count(1), recv_displ(1), R_MPITYPE, &
    geo%atoms_dist%mpi_grp%comm, mpi_err)
#endif
  
  SAFE_DEALLOCATE_A(recv_count)
  SAFE_DEALLOCATE_A(recv_displ)
  SAFE_DEALLOCATE_A(force_local)

  call profiling_out(prof_comm)
  POP_SUB(X(forces_gather))
end subroutine X(forces_gather)

!---------------------------------------------------------------------------
subroutine X(forces_from_local_potential)(gr, geo, ep, gdensity, force)
  type(grid_t),                   intent(inout) :: gr
  type(geometry_t),               intent(inout) :: geo
  type(epot_t),                   intent(inout) :: ep
  R_TYPE,                         intent(in)    :: gdensity(:, :)
  R_TYPE,                         intent(inout) :: force(:, :)

  FLOAT,  allocatable :: vloc(:)
  R_TYPE, pointer     :: zvloc(:)
  integer             :: ip, idir, iatom
 
  PUSH_SUB(X(forces_from_local_potential))

  SAFE_ALLOCATE(vloc(1:gr%mesh%np))
  SAFE_ALLOCATE(zvloc(1:gr%mesh%np))
  
  do iatom = geo%atoms_dist%start, geo%atoms_dist%end

    if(.not.simul_box_in_box(gr%mesh%sb, geo, geo%atom(iatom)%x) .and. ep%ignore_external_ions) then
      force(1:gr%mesh%sb%dim, iatom) = M_ZERO
      cycle
    end if
    
    vloc(1:gr%mesh%np) = M_ZERO
    call epot_local_potential(ep, gr%der, gr%dgrid, geo, iatom, vloc)

    forall(ip = 1:gr%mesh%np) zvloc(ip) = vloc(ip)

    do idir = 1, gr%mesh%sb%dim
      force(idir, iatom) = -X(mf_dotp)(gr%mesh, zvloc, gdensity(:, idir))
    end do

  end do

  if(geo%atoms_dist%parallel) call X(forces_gather)(geo, force)
  !if(geo%atoms_dist%parallel .and. geo%atoms_dist%nlocal > 0) call X(forces_gather)(geo, force)

  SAFE_DEALLOCATE_A(vloc)
  SAFE_DEALLOCATE_P(zvloc)

  POP_SUB(X(forces_from_local_potential))
end subroutine X(forces_from_local_potential)


!---------------------------------------------------------------------------
subroutine X(total_force_from_local_potential)(gr, ep, gdensity, force)
  type(grid_t),                   intent(inout) :: gr
  type(epot_t),                   intent(inout) :: ep
  R_TYPE,                         intent(in)    :: gdensity(:, :)
  R_TYPE,                         intent(inout) :: force(:)

  R_TYPE, pointer     :: zvloc(:)
  integer             :: idir
 
  PUSH_SUB(X(total_force_from_local_potential))
  SAFE_ALLOCATE(zvloc(1:gr%mesh%np))
  
  zvloc(1:gr%mesh%np) = ep%vpsl(1:gr%mesh%np)
  do idir = 1, gr%mesh%sb%dim
    force(idir) = force(idir) + X(mf_dotp)(gr%mesh, zvloc, gdensity(:, idir))
  end do

  SAFE_DEALLOCATE_P(zvloc)
  POP_SUB(X(total_force_from_local_potential))
end subroutine X(total_force_from_local_potential)


!---------------------------------------------------------------------------
!> Ref: Kikuji Hirose, Tomoya Ono, Yoshitaka Fujimoto, and Shigeru Tsukamoto,
!! First-principles calculations in real-space formalism: Electronic configurations
!! and transport properties of nanostructures, Imperial College Press (2005)
!! Section 1.6, page 12
subroutine X(forces_from_potential)(gr, geo, hm, st, force)
  type(grid_t),                   intent(inout) :: gr
  type(geometry_t),               intent(inout) :: geo
  type(hamiltonian_t),            intent(inout) :: hm
  type(states_t),                 intent(inout) :: st
  FLOAT,                          intent(out)   :: force(:, :)

  type(symmetrizer_t) :: symmetrizer
  integer :: iatom, ist, iq, idim, idir, np, np_part, ikpoint, iop, ii, iatom_symm
  integer :: ib, maxst, minst
  FLOAT :: ratom(1:MAX_DIM)
  R_TYPE, allocatable :: psi(:, :)
  R_TYPE, allocatable :: grad_psi(:, :, :)
  FLOAT,  allocatable :: grad_rho(:, :), force_loc(:, :), force_psi(:), force_tmp(:)
  FLOAT, allocatable :: symmtmp(:, :)
  type(batch_t) :: psib, grad_psib(1:MAX_DIM)

  PUSH_SUB(X(forces_from_potential))

  np = gr%mesh%np
  np_part = gr%mesh%np_part

  SAFE_ALLOCATE(grad_psi(1:np, 1:gr%mesh%sb%dim, 1:st%d%dim))
  SAFE_ALLOCATE(grad_rho(1:np, 1:gr%mesh%sb%dim))
  grad_rho = M_ZERO

  SAFE_ALLOCATE(force_psi(1:gr%mesh%sb%dim))
  SAFE_ALLOCATE(force_tmp(1:gr%mesh%sb%dim))

  force = M_ZERO

  ! even if there is no fine mesh, we need to make another copy
  SAFE_ALLOCATE(psi(1:np_part, 1:st%d%dim))

  !THE NON-LOCAL PART (parallel in states and k-points)
  do iq = st%d%kpt%start, st%d%kpt%end

    ikpoint = states_dim_get_kpoint_index(st%d, iq)

    do ib = st%group%block_start, st%group%block_end
      minst = states_block_min(st, ib)
      maxst = states_block_max(st, ib)

      call batch_copy(st%group%psib(ib, iq), psib)
      call batch_copy_data(gr%mesh%np, st%group%psib(ib, iq), psib)

      ! set the boundary conditions
      call boundaries_set(gr%der%boundaries, psib)

      ! set the phase for periodic systems
      if(simul_box_is_periodic(gr%sb) .and. .not. kpoints_point_is_gamma(gr%sb%kpoints, ikpoint)) then
        call X(hamiltonian_base_phase)(hm%hm_base, gr%der, gr%mesh%np_part, iq, .false., psib)
      end if

      ! calculate the gradient
      do idir = 1, gr%mesh%sb%dim
        call batch_copy(st%group%psib(ib, iq), grad_psib(idir))
        call X(derivatives_batch_perform)(gr%der%grad(idir), gr%der, psib, grad_psib(idir), set_bc = .false.)
      end do

      ! calculate the contribution to the density gradient
      call X(density_accumulate_grad)(gr, st, iq, psib, grad_psib, grad_rho)

      ! the non-local potential contribution
      if(hm%hm_base%apply_projector_matrices .and. .not. accel_is_enabled() .and. &
        .not. (st%symmetrize_density .and. gr%sb%kpoints%use_symmetries)) then

        call X(hamiltonian_base_nlocal_force)(hm%hm_base, gr%mesh, st, geo, iq, gr%mesh%sb%dim, psib, grad_psib, force)

      else 

        do ist = minst, maxst

          ! get the state and its gradient out of the batches (for the moment)
          do idim = 1, st%d%dim
            call batch_get_state(psib, (/ist, idim/), gr%mesh%np_part, psi(:, idim))
            do idir = 1, gr%mesh%sb%dim
              call batch_get_state(grad_psib(idir), (/ist, idim/), gr%mesh%np, grad_psi(:, idir, idim))
            end do
          end do

          call profiling_count_operations(np*st%d%dim*gr%mesh%sb%dim*(2 + R_MUL))

          if(st%symmetrize_density .and. gr%sb%kpoints%use_symmetries) then

            ! We use that
            !
            ! \int dr f(Rr) V_iatom(r) \nabla f(R(v)) = R\int dr f(r) V_iatom(R*r) f(r)
            !
            ! and that the operator R should map the position of atom
            ! iatom to the position of some other atom iatom_symm, so that
            !
            ! V_iatom(R*r) = V_iatom_symm(r)
            !
            do ii = 1, kpoints_get_num_symmetry_ops(gr%sb%kpoints, ikpoint)

              iop = kpoints_get_symmetry_ops(gr%sb%kpoints, ikpoint, ii)

              do iatom = 1, geo%natoms
                if(projector_is_null(hm%ep%proj(iatom))) cycle

                ratom = M_ZERO
                ratom(1:gr%sb%dim) = symm_op_apply_inv(gr%sb%symm%ops(iop), geo%atom(iatom)%x)

                call simul_box_periodic_atom_in_box(gr%sb, geo, ratom)

                ! find iatom_symm
                do iatom_symm = 1, geo%natoms
                  if(all(abs(ratom(1:gr%sb%dim) - geo%atom(iatom_symm)%x(1:gr%sb%dim)) < CNST(1.0e-5))) exit
                end do

                if(iatom_symm > geo%natoms) then
                  write(message(1),'(a,i6)') 'Internal error: could not find symmetric partner for atom number', iatom
                  write(message(2),'(a,i3,a)') 'with symmetry operation number ', iop, '.'
                  call messages_fatal(2)
                end if

                do idir = 1, gr%mesh%sb%dim
                  force_psi(idir) = - M_TWO * st%d%kweights(iq) * st%occ(ist, iq) * &
                    R_REAL(X(projector_matrix_element)(hm%ep%proj(iatom_symm), st%d%dim, iq, psi, grad_psi(:, idir, :)))
                end do


                force_tmp = symm_op_apply(gr%sb%symm%ops(iop), force_psi)

                force(1:gr%mesh%sb%dim, iatom) = force(1:gr%mesh%sb%dim, iatom) + &
                  force_tmp(1:gr%mesh%sb%dim)/kpoints_get_num_symmetry_ops(gr%sb%kpoints, ikpoint)

              end do

            end do

          else

            ! iterate over the projectors
            do iatom = 1, geo%natoms
              if(projector_is_null(hm%ep%proj(iatom))) cycle

              do idir = 1, gr%mesh%sb%dim
                force_psi(idir) = - M_TWO * st%d%kweights(iq) * st%occ(ist, iq) * &
                  R_REAL(X(projector_matrix_element)(hm%ep%proj(iatom), st%d%dim, iq, psi, grad_psi(:, idir, :)))
              end do

              force(1:gr%mesh%sb%dim, iatom) = force(1:gr%mesh%sb%dim, iatom) + force_psi(1:gr%mesh%sb%dim)
            end do

          end if

        end do

      end if


      call batch_end(psib)
      do idir = 1, gr%mesh%sb%dim
        call batch_end(grad_psib(idir))
      end do

    end do
  end do

  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(grad_psi)

#if defined(HAVE_MPI)
  if(st%parallel_in_states .or. st%d%kpt%parallel) then
    call profiling_in(prof_comm, "FORCES_COMM")
    call comm_allreduce(st%st_kpt_mpi_grp%comm, force)
    call comm_allreduce(st%st_kpt_mpi_grp%comm, grad_rho)
    call profiling_out(prof_comm)
  end if
#endif

  SAFE_ALLOCATE(force_loc(1:gr%mesh%sb%dim, 1:geo%natoms))

  if(st%symmetrize_density) then
    call symmetrizer_init(symmetrizer, gr%mesh)
    SAFE_ALLOCATE(symmtmp(1:gr%mesh%np, 1:3))

    call dsymmetrizer_apply(symmetrizer, field_vector = grad_rho, symmfield_vector = symmtmp, suppress_warning = .true.)
    grad_rho(1:gr%mesh%np, 1:3) = symmtmp(1:gr%mesh%np, 1:3)

    SAFE_DEALLOCATE_A(symmtmp)
    call symmetrizer_end(symmetrizer)
  end if

  call dforces_from_local_potential(gr, geo, hm%ep, grad_rho, force_loc)

  do iatom = 1, geo%natoms
    do idir = 1, gr%mesh%sb%dim
      force(idir, iatom) = force(idir, iatom) + force_loc(idir, iatom)
    end do
  end do


  ! transformation from reduced coordinates to Cartesian coordinates
  if (simul_box_is_periodic(gr%mesh%sb) .and. gr%mesh%sb%nonorthogonal ) then
     forall (iatom = 1:geo%natoms)
        force(1:gr%mesh%sb%dim,iatom) = &
          matmul(gr%mesh%sb%klattice_primitive(1:gr%mesh%sb%dim, 1:gr%mesh%sb%dim),force(1:gr%mesh%sb%dim,iatom))
     end forall
  end if

  SAFE_DEALLOCATE_A(force_tmp)
  SAFE_DEALLOCATE_A(force_psi)
  SAFE_DEALLOCATE_A(force_loc)
  SAFE_DEALLOCATE_A(grad_rho)
  
  POP_SUB(X(forces_from_potential))
end subroutine X(forces_from_potential)

!---------------------------------------------------------------------------
subroutine X(total_force_from_potential)(gr, geo, ep, st, x)
  type(grid_t),                   intent(inout) :: gr
  type(geometry_t),               intent(in)    :: geo
  type(epot_t),                   intent(inout) :: ep
  type(states_t),                 intent(inout) :: st
  FLOAT,                          intent(inout) :: x(1:MAX_DIM)
 
  integer :: iatom, ist, iq, idim, idir, np, np_part, ip, ikpoint
  FLOAT :: ff, kpoint(1:MAX_DIM)
  R_TYPE, allocatable :: psi(:, :)
  R_TYPE, allocatable :: grad_psi(:, :, :)
  FLOAT,  allocatable :: grad_rho(:, :), force(:, :)
  CMPLX :: phase

  PUSH_SUB(X(total_force_from_potential))

  np = gr%mesh%np
  np_part = gr%mesh%np_part

  SAFE_ALLOCATE(grad_psi(1:np, 1:gr%mesh%sb%dim, 1:st%d%dim))
  SAFE_ALLOCATE(grad_rho(1:np, 1:gr%mesh%sb%dim))
  grad_rho = M_ZERO
  SAFE_ALLOCATE(force(1:gr%mesh%sb%dim, 1:geo%natoms))
  force = M_ZERO

  ! even if there is no fine mesh, we need to make another copy
  SAFE_ALLOCATE(psi(1:np_part, 1:st%d%dim))

  !THE NON-LOCAL PART (parallel in states and k-points)
  do iq = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end

      call states_get_state(st, gr%mesh, ist, iq, psi)

      do idim = 1, st%d%dim
        call boundaries_set(gr%der%boundaries, psi(:, idim))

        ikpoint = states_dim_get_kpoint_index(st%d, iq)
        if(simul_box_is_periodic(gr%sb) .and. .not. kpoints_point_is_gamma(gr%sb%kpoints, ikpoint)) then

          kpoint = M_ZERO
          kpoint(1:gr%sb%dim) = kpoints_get_point(gr%sb%kpoints, ikpoint)

          do ip = 1, np_part
            phase = exp(-M_zI*sum(kpoint(1:gr%sb%dim)*gr%mesh%x(ip, 1:gr%sb%dim)))
            psi(ip, idim) = phase*psi(ip, idim)
          end do
        end if

        call X(derivatives_grad)(gr%der, psi(:, idim), grad_psi(:, :, idim), set_bc = .false.)

        ff = st%d%kweights(iq) * st%occ(ist, iq) * M_TWO
        do idir = 1, gr%mesh%sb%dim
          do ip = 1, np
            grad_rho(ip, idir) = grad_rho(ip, idir) + ff*R_REAL(R_CONJ(psi(ip, idim))*grad_psi(ip, idir, idim))
          end do
        end do

      end do

      call profiling_count_operations(np*st%d%dim*gr%mesh%sb%dim*(2 + R_MUL))

      ! iterate over the projectors
      do iatom = 1, geo%natoms
        if(projector_is_null(ep%proj(iatom))) cycle
        do idir = 1, gr%mesh%sb%dim

          force(idir, iatom) = force(idir, iatom) - M_TWO * st%d%kweights(iq) * st%occ(ist, iq) * &
            R_REAL(X(projector_matrix_element)(ep%proj(iatom), st%d%dim, iq, psi, grad_psi(:, idir, :)))

        end do
      end do

    end do
  end do

  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(grad_psi)

#if defined(HAVE_MPI)
  if(st%parallel_in_states .or. st%d%kpt%parallel) then
    call profiling_in(prof_comm, "FORCES_COMM")
    call comm_allreduce(st%st_kpt_mpi_grp%comm, force)
    call comm_allreduce(st%st_kpt_mpi_grp%comm, grad_rho)
    call profiling_out(prof_comm)
  end if
#endif

  call dtotal_force_from_local_potential(gr, ep, grad_rho, x)

  do iatom = 1, geo%natoms
    do idir = 1, gr%mesh%sb%dim
      x(idir) = x(idir) - force(idir, iatom)
    end do
  end do

  SAFE_DEALLOCATE_A(force)
  POP_SUB(X(total_force_from_potential))
end subroutine X(total_force_from_potential)


! --------------------------------------------------------------------------------
subroutine X(forces_derivative)(gr, geo, ep, st, lr, lr2, force_deriv)
  type(grid_t),                   intent(inout) :: gr
  type(geometry_t),               intent(inout) :: geo
  type(epot_t),                   intent(inout) :: ep
  type(states_t),                 intent(inout) :: st
  type(lr_t),                     intent(in)    :: lr
  type(lr_t),                     intent(in)    :: lr2
  CMPLX,                          intent(out)   :: force_deriv(:,:) !< (gr%mesh%sb%dim, geo%natoms)

  integer :: iatom, ist, iq, idim, idir, np, np_part, ip, ikpoint
  FLOAT :: ff, kpoint(1:MAX_DIM)
  R_TYPE, allocatable :: psi(:, :)
  R_TYPE, allocatable :: dl_psi(:, :)
  R_TYPE, allocatable :: dl_psi2(:, :)
  R_TYPE, allocatable :: grad_psi(:, :, :)
  R_TYPE, allocatable :: grad_dl_psi(:, :, :)
  R_TYPE, allocatable :: grad_dl_psi2(:, :, :)
  CMPLX,  allocatable :: grad_rho(:, :)
  CMPLX :: phase
  CMPLX, allocatable  :: force_local(:, :)

  PUSH_SUB(X(forces_derivative))

  np      = gr%mesh%np
  np_part = gr%mesh%np_part

  SAFE_ALLOCATE(grad_dl_psi(1:np, 1:gr%mesh%sb%dim, 1:st%d%dim))
  SAFE_ALLOCATE(grad_dl_psi2(1:np, 1:gr%mesh%sb%dim, 1:st%d%dim))
  SAFE_ALLOCATE(grad_psi(1:np, 1:gr%mesh%sb%dim, 1:st%d%dim))
  SAFE_ALLOCATE(grad_rho(1:np, 1:gr%mesh%sb%dim))
  grad_rho = M_ZERO
  force_deriv = M_ZERO

  ! even if there is no fine mesh, we need to make another copy
  SAFE_ALLOCATE(psi(1:np_part, 1:st%d%dim))
  SAFE_ALLOCATE(dl_psi(1:np_part, 1:st%d%dim))
  SAFE_ALLOCATE(dl_psi2(1:np_part, 1:st%d%dim))

  !THE NON-LOCAL PART (parallel in states and k-points)
  do iq = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end
      do idim = 1, st%d%dim

        call states_get_state(st, gr%mesh, idim, ist, iq, psi(:, idim))
        call boundaries_set(gr%der%boundaries, psi(:, idim))
        call lalg_copy(gr%mesh%np_part, lr%X(dl_psi)(:, idim, ist, iq), dl_psi(:, idim))
        call boundaries_set(gr%der%boundaries, dl_psi(:, idim))
        call lalg_copy(gr%mesh%np_part, lr2%X(dl_psi)(:, idim, ist, iq), dl_psi2(:, idim))
        call boundaries_set(gr%der%boundaries, dl_psi2(:, idim))

        ikpoint = states_dim_get_kpoint_index(st%d, iq)
        if(simul_box_is_periodic(gr%sb) .and. .not. kpoints_point_is_gamma(gr%sb%kpoints, ikpoint)) then

          kpoint = M_ZERO
          kpoint(1:gr%sb%dim) = kpoints_get_point(gr%sb%kpoints, ikpoint)

          do ip = 1, np_part
            phase = exp(-M_zI*sum(kpoint(1:gr%sb%dim)*gr%mesh%x(ip, 1:gr%sb%dim)))
            psi(ip, idim) = phase*psi(ip, idim)
            dl_psi(ip, idim) = phase*dl_psi(ip, idim)
            dl_psi2(ip, idim) = phase*dl_psi2(ip, idim)
          end do
        end if
       
        call X(derivatives_grad)(gr%der, psi(:, idim), grad_psi(:, :, idim), set_bc = .false.)
        call X(derivatives_grad)(gr%der, dl_psi(:, idim), grad_dl_psi(:, :, idim), set_bc = .false.)
        call X(derivatives_grad)(gr%der, dl_psi2(:, idim), grad_dl_psi2(:, :, idim), set_bc = .false.)

        !accumulate to calculate the gradient of the density
        ff = st%d%kweights(iq) * st%occ(ist, iq)
        do idir = 1, gr%mesh%sb%dim
          do ip = 1, np
            grad_rho(ip, idir) = grad_rho(ip, idir) + ff * &
              (R_CONJ(grad_psi(ip, idir, idim)) * dl_psi(ip, idim) + R_CONJ(psi(ip, idim)) * grad_dl_psi(ip, idir, idim) &
              + R_CONJ(dl_psi2(ip, idim)) * grad_psi(ip, idir, idim) + R_CONJ(grad_dl_psi2(ip, idir, idim)) * psi(ip, idim))
          end do
        end do
      end do

      ! iterate over the projectors
      do iatom = 1, geo%natoms
        if(projector_is_null(ep%proj(iatom))) cycle
        do idir = 1, gr%mesh%sb%dim

          force_deriv(idir, iatom) = force_deriv(idir, iatom) - st%d%kweights(iq) * st%occ(ist, iq) * &
            (X(projector_matrix_element)(ep%proj(iatom), st%d%dim, iq, grad_psi(:, idir, :), dl_psi) &
            + X(projector_matrix_element)(ep%proj(iatom), st%d%dim, iq, psi, grad_dl_psi(:, idir, :)) &
            + X(projector_matrix_element)(ep%proj(iatom), st%d%dim, iq, dl_psi2, grad_psi(:, idir, :)) &
            + X(projector_matrix_element)(ep%proj(iatom), st%d%dim, iq, grad_dl_psi2(:, idir, :), psi))
          
        end do
      end do

    end do
  end do
  
  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(dl_psi)
  SAFE_DEALLOCATE_A(dl_psi2)
  SAFE_DEALLOCATE_A(grad_psi)
  SAFE_DEALLOCATE_A(grad_dl_psi)
  SAFE_DEALLOCATE_A(grad_dl_psi2)

#if defined(HAVE_MPI)
  if(st%parallel_in_states .or. st%d%kpt%parallel) then
    call profiling_in(prof_comm, "FORCES_COMM")
    call comm_allreduce(st%st_kpt_mpi_grp%comm, force_deriv, dim = (/gr%mesh%sb%dim, geo%natoms/))
    call comm_allreduce(st%st_kpt_mpi_grp%comm, grad_rho)
    call profiling_out(prof_comm)
  end if
#endif
  
  SAFE_ALLOCATE(force_local(1:gr%sb%dim, 1:geo%natoms))
  call zforces_from_local_potential(gr, geo, ep, grad_rho, force_local)
  force_deriv(:,:) = force_deriv(:,:) + force_local(:,:)
  SAFE_DEALLOCATE_A(force_local)
  SAFE_DEALLOCATE_A(grad_rho)

  POP_SUB(X(forces_derivative)) 
end subroutine X(forces_derivative)

! --------------------------------------------------------------------------------
!> lr, lr2 are wfns from electric perturbation; lr is for +omega, lr2 is for -omega.
!! for each atom, Z*(i,j) = dF(j)/dE(i)
subroutine X(forces_born_charges)(gr, geo, ep, st, lr, lr2, born_charges)
  type(grid_t),                   intent(inout) :: gr
  type(geometry_t),               intent(inout) :: geo
  type(epot_t),                   intent(inout) :: ep
  type(states_t),                 intent(inout) :: st
  type(lr_t),                     intent(in)    :: lr(:)  !< (gr%mesh%sb%dim)
  type(lr_t),                     intent(in)    :: lr2(:) !< (gr%mesh%sb%dim)
  type(born_charges_t),           intent(inout) :: born_charges

  integer :: iatom, idir
  CMPLX, allocatable :: force_deriv(:, :)

  PUSH_SUB(X(forces_born_charges))

  SAFE_ALLOCATE(force_deriv(1:gr%mesh%sb%dim, 1:geo%natoms))

  do idir = 1, gr%mesh%sb%dim
    call X(forces_derivative)(gr, geo, ep, st, lr(idir), lr2(idir), force_deriv)
    do iatom = 1, geo%natoms
      born_charges%charge(:, idir, iatom) = force_deriv(:, iatom)
      born_charges%charge(idir, idir, iatom) = born_charges%charge(idir, idir, iatom) + species_zval(geo%atom(iatom)%species)
    end do
  end do

  SAFE_DEALLOCATE_A(force_deriv)

  do iatom = 1, geo%natoms
    call zsymmetrize_tensor(gr%sb%symm, born_charges%charge(:, :, iatom))
  end do

  POP_SUB(X(forces_born_charges))
end subroutine X(forces_born_charges)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
