!! Copyright (C) 2015 P. Wopperer and U. De Giovannini
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


  ! ---------------------------------------------------------
  subroutine pes_flux_output(this, mesh, sb, st, dt)
    type(pes_flux_t), intent(inout)    :: this
    type(mesh_t),        intent(in)    :: mesh
    type(simul_box_t),   intent(in)    :: sb
    type(states_t),      intent(in)    :: st
    FLOAT,               intent(in)    :: dt

    integer            :: stst, stend, kptst, kptend, sdim, mdim
    integer            :: ist, ik, isdim
    integer            :: ikp, iomk, ikp_save, iomk_save
    integer            :: ikk, ith, iph, iphi
    FLOAT              :: phik, thetak, kact

    integer            :: iunitone, iunittwo
    FLOAT, allocatable :: spctrout_cub(:), spctrout_sph(:,:)
    FLOAT              :: weight, spctrsum

    PUSH_SUB(pes_flux_output)

    stst   = st%st_start
    stend  = st%st_end
    kptst  = st%d%kpt%start
    kptend = st%d%kpt%end
    sdim   = st%d%dim
    mdim   = mesh%sb%dim

    if(this%shape == M_SPHERICAL) then
      SAFE_ALLOCATE(spctrout_sph(1:this%nk, 1:this%nstepsomegak))
      spctrout_sph = M_ZERO
    else
      SAFE_ALLOCATE(spctrout_cub(1:this%nkpnts))
      spctrout_cub = M_ZERO
    end if

    ! calculate the total spectrum
    do ik = kptst, kptend
      do ist = stst, stend
        do isdim = 1, sdim
          if(this%shape == M_SPHERICAL) then
            spctrout_sph(1:this%nk, 1:this%nstepsomegak) = spctrout_sph(1:this%nk, 1:this%nstepsomegak) + &
              abs(this%spctramp_sph(ist, isdim, ik, 1:this%nk, 1:this%nstepsomegak))**M_TWO * (dt * this%tdstepsinterval)**M_TWO
          else
            spctrout_cub(1:this%nkpnts) = spctrout_cub(1:this%nkpnts) + &
              abs(this%spctramp_cub(ist, isdim, ik, 1:this%nkpnts))**M_TWO * (dt * this%tdstepsinterval)**M_TWO
          end if
        end do
      end do
    end do

    if(st%parallel_in_states .or. st%d%kpt%parallel) then
#if defined(HAVE_MPI)
      if(this%shape == M_SPHERICAL) then
        call comm_allreduce(st%st_kpt_mpi_grp%comm, spctrout_sph)
      else
        call comm_allreduce(st%st_kpt_mpi_grp%comm, spctrout_cub)
      end if
#endif
    end if

    if(mpi_grp_is_root(mpi_world)) then
      iunittwo = io_open('td.general/PES_flux.distribution.out', action='write', position='rewind')
      iunitone = io_open('td.general/'//'PES_flux.power.sum', action='write', position='rewind')
      write(iunitone, '(a19)') '# E, total spectrum'

      if(this%shape == M_SPHERICAL) then
        write(iunittwo, '(a29)') '# k, theta, phi, distribution'
        do ikk = 1, this%nk 
          kact = ikk * this%dk
          iomk = 0
          spctrsum = M_ZERO

          do ith = 0, this%nstepsthetak
            thetak = ith * M_PI / this%nstepsthetak 

            if(ith == 0 .or. ith == this%nstepsthetak) then
              weight = (M_ONE - cos(M_PI / this%nstepsthetak / M_TWO)) * M_TWO * M_PI
            else
              weight = abs(cos(thetak - M_PI / this%nstepsthetak / M_TWO) - cos(thetak + M_PI / this%nstepsthetak / M_TWO)) &
                * M_TWO * M_PI / this%nstepsphik
            end if

            do iph = 0, this%nstepsphik - 1
              iomk = iomk + 1
              spctrsum = spctrsum + spctrout_sph(ikk, iomk) * weight 
              phik = iph * M_TWO * M_PI / this%nstepsphik
              if(iph == 0) iomk_save = iomk
              write(iunittwo,'(5(1x,e18.10E3))') kact, thetak, phik, spctrout_sph(ikk, iomk)

              ! just repeat the result for output
              if(this%nstepsphik > 1 .and. iph == (this%nstepsphik - 1)) &
                write(iunittwo,'(5(1x,e18.10E3))') kact, thetak, M_TWO * M_PI, spctrout_sph(ikk, iomk_save)

              ! just repeat the result for output and exit
              if(ith == 0 .or. ith == this%nstepsthetak) then
                if(this%nstepsphik > 1) then
                  do iphi = 1, this%nstepsphik
                    phik = iphi * M_TWO * M_PI / this%nstepsphik
                    write(iunittwo,'(5(1x,e18.10E3))') kact, thetak, phik, spctrout_sph(ikk, iomk)
                  end do
                end if
                exit
              end if
            end do

            if(this%nstepsphik > 1 .or. ith == this%nstepsthetak) write(iunittwo, '(1x)', advance='yes')
          end do
          write(iunitone, '(2(1x,e18.10E3))') kact**M_TWO / M_TWO, spctrsum * kact
        end do

      else ! this%shape == M_CUBIC
        select case(mdim)
        case(1)
          write(iunittwo, '(a17)') '# k, distribution'
          do ikp = 1, this%nkpnts
            write(iunittwo, '(5(1x,e18.10E3))') this%kcoords_cub(1, ikp), spctrout_cub(ikp)
          end do

          do ikp = this%nk + 1, this%nkpnts
            kact = this%kcoords_cub(1, ikp)
            write(iunitone, '(2(1x,e18.10E3))') kact**M_TWO / M_TWO, &
              (spctrout_cub(ikp) + spctrout_cub(this%nkpnts + 1 - ikp)) / M_TWO * kact
          end do

        case(2)
          write(iunittwo, '(a29)') '# k, phi, distribution'
          ikp = 0
          do ikk = 1, this%nk
            kact = ikk * this%dk
            
            spctrsum = M_ZERO
            do iph = 0, this%nstepsphik - 1
              ikp = ikp + 1
              if(iph == 0) ikp_save = ikp
              spctrsum = spctrsum + spctrout_cub(ikp) * M_TWO * M_PI / this%nstepsphik
              phik = iph * M_TWO * M_PI / this%nstepsphik
              write(iunittwo,'(5(1x,e18.10E3))') kact, phik, spctrout_cub(ikp)
            end do
            ! just repeat the result for output
            write(iunittwo,'(5(1x,e18.10E3))') kact, M_TWO * M_PI, spctrout_cub(ikp_save)
            write(iunittwo,'(1x)', advance = 'yes')
            write(iunitone, '(2(1x,e18.10E3))') kact**M_TWO / M_TWO, spctrsum * kact
          end do

        case(3)
          write(iunittwo, '(a29)') '# k, theta, phi, distribution'
          ikp    = 0
          do ikk = 1, this%nk
            kact = ikk * this%dk
            spctrsum = M_ZERO

            do ith = 0, this%nstepsthetak
              thetak = ith * M_PI / this%nstepsthetak 

              if(ith == 0 .or. ith == this%nstepsthetak) then
                weight = (M_ONE - cos(M_PI / this%nstepsthetak / M_TWO)) * M_TWO * M_PI
              else
                weight = abs(cos(thetak - M_PI / this%nstepsthetak / M_TWO) - cos(thetak + M_PI / this%nstepsthetak / M_TWO)) &
                  * M_TWO * M_PI / this%nstepsphik
              end if

              do iph = 0, this%nstepsphik - 1
                ikp = ikp + 1
                spctrsum = spctrsum + spctrout_cub(ikp) * weight

                phik = iph * M_TWO * M_PI / this%nstepsphik
                if(iph == 0) ikp_save = ikp
                write(iunittwo,'(5(1x,e18.10E3))') kact, thetak, phik, spctrout_cub(ikp)

                ! just repeat the result for output
                if(iph == (this%nstepsphik - 1)) &
                  write(iunittwo,'(5(1x,e18.10E3))') kact, thetak, M_TWO * M_PI, spctrout_cub(ikp_save)

                ! just repeat the result for output and exit
                if(ith == 0 .or. ith == this%nstepsthetak) then
                  do iphi = 1, this%nstepsphik
                    phik = iphi * M_TWO * M_PI / this%nstepsphik
                    write(iunittwo,'(5(1x,e18.10E3))') kact, thetak, phik, spctrout_cub(ikp)
                  end do
                  exit
                end if
              end do

              write(iunittwo, '(1x)', advance='yes')
            end do
            write(iunitone, '(2(1x,e18.10E3))') kact**M_TWO / M_TWO, spctrsum * kact
          end do
        end select
      end if

      call io_close(iunittwo)
      call io_close(iunitone)
    end if

    SAFE_DEALLOCATE_A(spctrout_cub)
    SAFE_DEALLOCATE_A(spctrout_sph)

    POP_SUB(pes_flux_output)
  end subroutine pes_flux_output

  ! ---------------------------------------------------------
  subroutine pes_flux_dump(restart, this, mesh, st, ierr)
    type(restart_t),  intent(in)  :: restart
    type(pes_flux_t), intent(in)  :: this
    type(mesh_t),     intent(in)  :: mesh
    type(states_t),   intent(in)  :: st
    integer,          intent(out) :: ierr

    integer          :: stst, stend, kptst, kptend, sdim, mdim
    integer          :: ist, ik, isdim, itot
    integer          :: err
    character(len=128) :: filename

    PUSH_SUB(pes_flux_dump)

    stst   = st%st_start
    stend  = st%st_end
    kptst  = st%d%kpt%start
    kptend = st%d%kpt%end
    sdim   = st%d%dim
    mdim   = mesh%sb%dim

    if(restart_skip(restart)) then
      POP_SUB(pes_flux_dump)
      return
    end if

    if(debug%info) then
      message(1) = "Debug: Writing pes_flux restart."
      call messages_info(1)
    end if

    do ik = kptst, kptend
      do ist = stst, stend
        do isdim = 1, sdim
!           write(filename, '(i2.2, a, i2.2, a, i2.2)') ik, '.', ist, '.', isdim
          itot = ist + (ik-1) * st%nst+  (isdim-1) * st%nst*st%d%kpt%nglobal
          write(filename,'(i10.10)') itot

          if(mpi_grp_is_root(mesh%mpi_grp)) then

            if(this%shape == M_SPHERICAL) then
              call io_binary_write(trim(restart_dir(restart))//"/pesflux1."//trim(filename)//".obf", &
                this%nk * this%nstepsomegak, this%spctramp_sph(ist, isdim, ik, :, :), err)
            else
              call io_binary_write(trim(restart_dir(restart))//"/pesflux1."//trim(filename)//".obf", &
                this%nkpnts, this%spctramp_cub(ist, isdim, ik, :), err)
            end if

            call io_binary_write(trim(restart_dir(restart))//"/pesflux2."//trim(filename)//".obf", &
              this%nsrfcpnts * this%tdsteps, this%wf(ist, isdim, ik, :, :), err)
            call io_binary_write(trim(restart_dir(restart))//"/pesflux3."//trim(filename)//".obf", &
              this%nsrfcpnts * this%tdsteps * mdim, this%gwf(ist, isdim, ik, :, :, :), err)

          end if
          
        end do
      end do
    end do

    if(this%shape == M_SPHERICAL) then
      call zrestart_write_binary(restart, 'pesflux4', this%nk * this%nstepsomegak, this%conjgphase_prev_sph, err)
    else
      call zrestart_write_binary(restart, 'pesflux4', this%nkpnts, this%conjgphase_prev_cub, err)
    end if

    call drestart_write_binary(restart, 'pesflux5', this%tdsteps * mdim, this%veca, err) 

    if(err /= 0) ierr = ierr + 1

    if(debug%info) then
      message(1) = "Debug: Writing pes_flux restart done."
      call messages_info(1)
    end if

    POP_SUB(pes_flux_dump)
  end subroutine pes_flux_dump

  ! ---------------------------------------------------------
  subroutine pes_flux_load(restart, this, mesh, st, ierr)
    type(restart_t),     intent(in)    :: restart
    type(pes_flux_t),    intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    type(states_t),      intent(in)    :: st
    integer,             intent(out)   :: ierr

    integer          :: stst, stend, kptst, kptend, sdim, mdim
    integer          :: ist, ik, isdim, itot
    integer          :: err
    character(len=128) :: filename

    PUSH_SUB(pes_flux_load)

    stst   = st%st_start
    stend  = st%st_end
    kptst  = st%d%kpt%start
    kptend = st%d%kpt%end
    sdim   = st%d%dim
    mdim   = mesh%sb%dim

    if(restart_skip(restart)) then
      ierr = -1
      POP_SUB(pes_flux_load)
      return
    end if

    if(debug%info) then
      message(1) = "Debug: Reading pes_flux restart."
      call messages_info(1)
    end if

    do ik = kptst, kptend
      do ist = stst, stend
        do isdim = 1, sdim
!           write(filename, '(i2.2, a, i2.2, a, i2.2)') ik, '.', ist, '.', isdim
          itot = ist + (ik-1) * st%nst+  (isdim-1) * st%nst*st%d%kpt%nglobal
          write(filename,'(i10.10)') itot
  
          if(this%shape == M_SPHERICAL) then
            call io_binary_read(trim(restart_dir(restart))//"/pesflux1."//trim(filename)//".obf", &
              this%nk * this%nstepsomegak, this%spctramp_sph(ist, isdim, ik, :, :), err)
          else
            call io_binary_read(trim(restart_dir(restart))//"/pesflux1."//trim(filename)//".obf", &
              this%nkpnts, this%spctramp_cub(ist, isdim, ik, :), err)
          end if

          call io_binary_read(trim(restart_dir(restart))//"/pesflux2."//trim(filename)//".obf", &
            this%nsrfcpnts * this%tdsteps, this%wf(ist, isdim, ik, :, :), err)
          call io_binary_read(trim(restart_dir(restart))//"/pesflux3."//trim(filename)//".obf", &
            this%nsrfcpnts * this%tdsteps * mdim, this%gwf(ist, isdim, ik, :, :, :), err)
        end do
      end do
    end do

    if(this%shape == M_SPHERICAL) then
      call zrestart_read_binary(restart, 'pesflux4', this%nk * this%nstepsomegak, this%conjgphase_prev_sph, err)
    else
      call zrestart_read_binary(restart, 'pesflux4', this%nkpnts, this%conjgphase_prev_cub, err)
    end if

    call drestart_read_binary(restart, 'pesflux5', this%tdsteps * mdim, this%veca, err) 

    if(err /= 0) ierr = ierr + 1
   
    if(debug%info) then
      message(1) = "Debug: Reading pes_flux restart done."
      call messages_info(1)
    end if

    POP_SUB(pes_flux_load)
  end subroutine pes_flux_load
