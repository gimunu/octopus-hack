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
!! $Id: td.F90 15582 2016-08-14 10:27:12Z philipp $

#include "global.h"

module td_oct_m
  use boundary_op_oct_m
  use calc_mode_par_oct_m
  use density_oct_m
  use energy_calc_oct_m
  use forces_oct_m
  use gauge_field_oct_m
  use geometry_oct_m
  use grid_oct_m
  use ground_state_oct_m
  use hamiltonian_oct_m
  use io_oct_m
  use ion_dynamics_oct_m
  use kick_oct_m
  use lasers_oct_m
  use loct_oct_m
  use loct_math_oct_m
  use modelmb_exchange_syms_oct_m
  use mpi_oct_m
  use parser_oct_m
  use pes_oct_m
  use poisson_oct_m
  use potential_interpolation_oct_m
  use profiling_oct_m
  use restart_oct_m
  use scdm_oct_m
  use scf_oct_m
  use simul_box_oct_m
  use states_oct_m
  use states_calc_oct_m
  use states_restart_oct_m
  use system_oct_m
  use propagator_oct_m
  use propagator_base_oct_m
  use td_write_oct_m
  use types_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use v_ks_oct_m
  use varinfo_oct_m

  implicit none

  private
  public ::               &
    td_t,                 &
    td_run,               &
    td_run_init,          &
    td_init,              &
    td_end,               &
    transform_states

  !> Parameters.
  integer, parameter :: &
    EHRENFEST = 1,   &
    BO        = 2

  type td_t
    type(propagator_t)   :: tr             !< contains the details of the time-evolution
    type(scf_t)          :: scf
    type(ion_dynamics_t) :: ions
    FLOAT                :: dt             !< time step
    integer              :: max_iter       !< maximum number of iterations to perform
    integer              :: iter           !< the actual iteration
    logical              :: recalculate_gs !< Recalculate ground-state along the evolution.

    type(pes_t)          :: pesv
    type(gauge_force_t)  :: gauge_force

    FLOAT                :: mu
    integer              :: dynamics
    integer              :: energy_update_iter
  end type td_t


contains

  subroutine td_run_init()

    PUSH_SUB(td_run_init)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .true.)

    POP_SUB(td_run_init)
  end subroutine td_run_init

  ! ---------------------------------------------------------

  subroutine td_init(td, sys, hm)
    type(td_t),            intent(inout) :: td
    type(system_t),        intent(inout) :: sys
    type(hamiltonian_t),   intent(inout) :: hm

    integer :: default
    FLOAT   :: spacing, default_dt, propagation_time

    PUSH_SUB(td_init)

    call ion_dynamics_init(td%ions, sys%geo)

    td%iter = 0

    !%Variable TDIonicTimeScale
    !%Type float
    !%Default 1.0
    !%Section Time-Dependent::Propagation
    !%Description
    !% This variable defines the factor between the timescale of ionic
    !% and electronic movement. It allows reasonably fast
    !% Born-Oppenheimer molecular-dynamics simulations based on
    !% Ehrenfest dynamics. The value of this variable is equivalent to
    !% the role of <math>\mu</math> in Car-Parrinello. Increasing it
    !% linearly accelerates the time step of the ion
    !% dynamics, but also increases the deviation of the system from the
    !% Born-Oppenheimer surface. The default is 1, which means that both
    !% timescales are the same. Note that a value different than 1
    !% implies that the electrons will not follow physical behaviour.
    !%
    !% According to our tests, values around 10 are reasonable, but it
    !% will depend on your system, mainly on the width of the gap.
    !%
    !% Important: The electronic time step will be the value of
    !% <tt>TDTimeStep</tt> divided by this variable, so if you have determined an
    !% optimal electronic time step (that we can call <i>dte</i>), it is
    !% recommended that you define your time step as:
    !%
    !% <tt>TDTimeStep</tt> = <i>dte</i> * <tt>TDIonicTimeScale</tt>
    !%
    !% so you will always use the optimal electronic time step
    !% (<a href=http://arxiv.org/abs/0710.3321>more details</a>).
    !%End
    call parse_variable('TDIonicTimeScale', CNST(1.0), td%mu)

    if (td%mu <= M_ZERO) then
      write(message(1),'(a)') 'Input: TDIonicTimeScale must be positive.'
      call messages_fatal(1)
    end if

    call messages_print_var_value(stdout, 'TDIonicTimeScale', td%mu)

    !%Variable TDTimeStep
    !%Type float
    !%Section Time-Dependent::Propagation
    !%Description
    !% The time-step for the time propagation. For most propagators you
    !% want to use the largest value that is possible without the
    !% evolution becoming unstable.
    !%
    !% The default value is the maximum value that we have found
    !% empirically that is stable for the spacing <math>h</math>:
    !% <math>dt = 0.0426 - 0.207 h + 0.808 h^2</math>
    !% (from parabolic fit to Fig. 4 of http://dx.doi.org/10.1021/ct800518j,
    !% probably valid for 3D systems only).
    !% However, you might need to adjust this value.
    !%End

    spacing = minval(sys%gr%mesh%spacing(1:sys%gr%sb%dim))
    default_dt = CNST(0.0426) - CNST(0.207)*spacing + CNST(0.808)*spacing**2
    default_dt = default_dt*td%mu

    call parse_variable('TDTimeStep', default_dt, td%dt, unit = units_inp%time)

    if (td%dt <= M_ZERO) then
      write(message(1),'(a)') 'Input: TDTimeStep must be positive.'
      call messages_fatal(1)
    end if

    call messages_print_var_value(stdout, 'TDTimeStep', td%dt, unit = units_out%time)

    td%dt = td%dt/td%mu
    
    if(parse_is_defined('TDMaxSteps') .and. parse_is_defined('TDPropagationTime')) then
      call messages_write('You cannot set TDMaxSteps and TDPropagationTime at the same time')
      call messages_fatal()
    end if

    !%Variable TDPropagationTime
    !%Type float
    !%Section Time-Dependent::Propagation
    !%Description
    !% The length of the time propagation. You cannot set this variable
    !% at the same time as <tt>TDMaxSteps</tt>. By default this variable will
    !% not be used.
    !%
    !% The units for this variable are <math>\hbar</math>/Hartree (or <math>\hbar</math>/eV if you
    !% selected <tt>ev_angstrom</tt> as input units). The approximate conversions to
    !% femtoseconds are 1 fs = 41.34 <math>\hbar</math>/Hartree = 1.52 <math>\hbar</math>/eV.
    !%End
    call parse_variable('TDPropagationTime', CNST(-1.0), propagation_time, unit = units_inp%time)

    call messages_obsolete_variable('TDMaximumIter', 'TDMaxSteps')

    !%Variable TDMaxSteps
    !%Type integer
    !%Default 1500
    !%Section Time-Dependent::Propagation
    !%Description
    !% Number of time-propagation steps that will be performed. You
    !% cannot use this variable together with <tt>TDPropagationTime</tt>.
    !%End
    default = 1500
    if(propagation_time > CNST(0.0)) default = nint(propagation_time/td%dt)
    call parse_variable('TDMaxSteps', default, td%max_iter)

    if(propagation_time <= CNST(0.0)) propagation_time = td%dt*td%max_iter

    call messages_print_var_value(stdout, 'TDPropagationTime', propagation_time, unit = units_out%time)
    call messages_print_var_value(stdout, 'TDMaxSteps', td%max_iter)

    if(td%max_iter < 1) then
      write(message(1), '(a,i6,a)') "Input: '", td%max_iter, "' is not a valid value for TDMaxSteps."
      message(2) = '(TDMaxSteps <= 1)'
      call messages_fatal(2)
    end if

    ! now the photoelectron stuff
    call pes_init(td%pesv, sys%gr%mesh, sys%gr%sb, sys%st, sys%outp%restart_write_interval, hm, td%max_iter, td%dt)

    !%Variable TDDynamics
    !%Type integer
    !%Default ehrenfest
    !%Section Time-Dependent::Propagation
    !%Description
    !% Type of dynamics to follow during a time propagation.
    !% For BO, you must set <tt>MoveIons = yes</tt>.
    !%Option ehrenfest 1
    !% Ehrenfest dynamics.
    !%Option bo 2
    !% Born-Oppenheimer (Experimental).
    !%End

    call parse_variable('TDDynamics', EHRENFEST, td%dynamics)
    if(.not.varinfo_valid_option('TDDynamics', td%dynamics)) call messages_input_error('TDDynamics')
    call messages_print_var_option(stdout, 'TDDynamics', td%dynamics)
    if(td%dynamics .ne. EHRENFEST) then
      if(.not.ion_dynamics_ions_move(td%ions)) call messages_input_error('TDDynamics')
    end if

    !%Variable RecalculateGSDuringEvolution
    !%Type logical
    !%Default no
    !%Section Time-Dependent::Propagation
    !%Description
    !% In order to calculate some information about the system along the
    !% evolution (e.g. projection onto the ground-state KS determinant,
    !% projection of the TDKS spin-orbitals onto the ground-state KS
    !% spin-orbitals), the ground-state KS orbitals are needed. If the
    !% ionic potential changes -- that is, the ions move -- one may want
    !% to recalculate the ground state. You may do this by setting this
    !% variable.
    !%
    !% The recalculation is not done every time step, but only every
    !% <tt>RestartWriteInterval</tt> time steps.
    !%End
    call parse_variable('RecalculateGSDuringEvolution', .false., td%recalculate_gs)

    call messages_obsolete_variable('TDScissor')

    call propagator_init(sys%gr, sys%st, td%tr, &
      ion_dynamics_ions_move(td%ions) .or. gauge_field_is_applied(hm%ep%gfield))
    if(hm%ep%no_lasers>0.and.mpi_grp_is_root(mpi_world)) then
      call messages_print_stress(stdout, "Time-dependent external fields")
      call laser_write_info(hm%ep%lasers, stdout, td%dt, td%max_iter)
      call messages_print_stress(stdout)
    end if

    !%Variable TDEnergyUpdateIter
    !%Type integer
    !%Section Time-Dependent::Propagation
    !%Description
    !% This variable controls after how many iterations Octopus
    !% updates the total energy during a time-propagation run. For
    !% iterations where the energy is not updated, the last calculated
    !% value is reported. If you set this variable to 1, the energy
    !% will be calculated in each step. The default value is 10,
    !% unless the ions move, in which case the default is 1.
    !%End

    default = 10
    if(ion_dynamics_ions_move(td%ions)) default = 1
    call parse_variable('TDEnergyUpdateIter', default, td%energy_update_iter)

    if(ion_dynamics_ions_move(td%ions) .and. td%energy_update_iter /= 1) then
      call messages_experimental('TDEnergyUpdateIter /= 1 when moving ions')
    end if
    
    POP_SUB(td_init)
  end subroutine td_init

  ! ---------------------------------------------------------
  
  subroutine td_end(td)
    type(td_t), intent(inout) :: td

    PUSH_SUB(td_end)

    call pes_end(td%pesv)
    call propagator_end(td%tr)  ! clean the evolution method
    call ion_dynamics_end(td%ions)

    if(td%dynamics == BO) call scf_end(td%scf)

    POP_SUB(td_end)
  end subroutine td_end

  ! ---------------------------------------------------------
  
  subroutine td_run(sys, hm, fromScratch)
    type(system_t), target, intent(inout) :: sys
    type(hamiltonian_t),    intent(inout) :: hm
    logical,                intent(inout) :: fromScratch

    type(td_t)                :: td
    type(td_write_t)          :: write_handler
    type(grid_t),     pointer :: gr   ! some shortcuts
    type(states_t),   pointer :: st
    type(geometry_t), pointer :: geo
    logical                   :: stopping, cmplxscl
    integer                   :: iter, ierr, scsteps
    real(8)                   :: etime
    type(profile_t),     save :: prof
    type(restart_t)           :: restart_load, restart_dump

    PUSH_SUB(td_run)

    cmplxscl = hm%cmplxscl%space

    ! some shortcuts
    gr  => sys%gr
    geo => sys%geo
    st  => sys%st

    if(simul_box_is_periodic(gr%mesh%sb)) call messages_experimental('Time propagation for periodic systems')

    call td_init(td, sys, hm)

    ! Allocate wavefunctions during time-propagation
    if(td%dynamics == EHRENFEST) then
      !complex wfs are required for Ehrenfest
      call states_allocate_wfns(st, gr%mesh, TYPE_CMPLX, alloc_Left = cmplxscl)
    else
      call states_allocate_wfns(st, gr%mesh, alloc_Left = cmplxscl)
      call scf_init(td%scf, sys%gr, sys%geo, sys%st, hm)
    end if

    if(hm%scdm_EXX) then
      call scdm_init(st,gr%der,psolver%cube, hm%scdm,operate_on_scdm=.true.)
      ! make sure scdm is constructed as soon as it is needed
      scdm_is_local = .false.
    end if
    
    if (gauge_field_is_applied(hm%ep%gfield)) then
      !if the gauge field is applied, we need to tell v_ks to calculate the current
      call v_ks_calculate_current(sys%ks, .true.)

      ! initialize the vector field and update the hamiltonian
      call gauge_field_init_vec_pot(hm%ep%gfield, gr%sb, st)
      call hamiltonian_update(hm, gr%mesh, time = td%dt*td%iter)
    end if

    call init_wfs()

    if(td%iter >= td%max_iter) then
      call end_()
      POP_SUB(td_run)
      return
    end if

    ! Calculate initial forces and kinetic energy
    if(ion_dynamics_ions_move(td%ions)) then
      if(td%iter > 0) then
        call td_read_coordinates()
        call hamiltonian_epot_generate(hm, gr, geo, st, time = td%iter*td%dt)
      end if

      call forces_calculate(gr, geo, hm, st, td%iter*td%dt, td%dt)

      geo%kinetic_energy = ion_dynamics_kinetic_energy(geo)
    end if

    call td_write_init(write_handler, gr, st, hm, geo, sys%ks, &
      ion_dynamics_ions_move(td%ions), gauge_field_is_applied(hm%ep%gfield), hm%ep%kick, td%iter, td%max_iter, td%dt)

    if(td%iter == 0) call td_run_zero_iter()

    if (gauge_field_is_applied(hm%ep%gfield)) call gauge_field_get_force(hm%ep%gfield, gr, st, td%gauge_force)

    !call td_check_trotter(td, sys, h)
    td%iter = td%iter + 1

    call restart_init(restart_dump, RESTART_TD, RESTART_TYPE_DUMP, st%dom_st_kpt_mpi_grp, ierr, mesh=gr%mesh)
    if (ion_dynamics_ions_move(td%ions) .and. td%recalculate_gs) then
      ! We will also use the TD restart directory as temporary storage during the time propagation
      call restart_init(restart_load, RESTART_TD, RESTART_TYPE_DUMP, st%dom_st_kpt_mpi_grp, &
        ierr, mesh=gr%mesh)
    end if

    call messages_print_stress(stdout, "Time-Dependent Simulation")
    call print_header()

    if(td%pesv%calc_spm .or. td%pesv%calc_mask .and. fromScratch) then
      call pes_init_write(td%pesv,gr%mesh,st)
    end if

    if(st%d%pack_states .and. hamiltonian_apply_packed(hm, gr%mesh)) call states_pack(st)

    etime = loct_clock()
    ! This is the time-propagation loop. It starts at t=0 and finishes
    ! at td%max_iter*dt. The index i runs from 1 to td%max_iter, and
    ! step "iter" means propagation from (iter-1)*dt to iter*dt.
    propagation: do iter = td%iter, td%max_iter

      stopping = clean_stop(sys%mc%master_comm)
      call profiling_in(prof, "TIME_STEP")

      if(iter > 1) then
        if( ((iter-1)*td%dt <= hm%ep%kick%time) .and. (iter*td%dt > hm%ep%kick%time) ) then
          if(.not. cmplxscl) then
            call kick_apply(gr%mesh, st, td%ions, geo, hm%ep%kick)
          else
            call kick_apply(gr%mesh, st, td%ions, geo, hm%ep%kick, hm%cmplxscl%theta)
          end if
          call td_write_kick(gr%mesh, hm%ep%kick, sys%outp, geo, iter)
        end if
      end if

      ! in case use scdm localized states for exact exchange and request a new localization             
      if(hm%scdm_EXX) scdm_is_local = .false.

      ! time iterate the system, one time step.
      select case(td%dynamics)
      case(EHRENFEST)
        call propagator_dt(sys%ks, hm, gr, st, td%tr, iter*td%dt, td%dt, td%energy_update_iter*td%mu, iter, td%ions, geo, &
          gauge_force = td%gauge_force, scsteps = scsteps, &
          update_energy = (mod(iter, td%energy_update_iter) == 0) .or. (iter == td%max_iter) )
      case(BO)
        call propagator_dt_bo(td%scf, gr, sys%ks, st, hm, td%gauge_force, geo, sys%mc, sys%outp, iter, td%dt, td%ions, scsteps)
      end select

      !Apply mask absorbing boundaries
      if(hm%bc%abtype == MASK_ABSORBING) call zvmask(gr, hm, st) 

      !Photoelectron stuff 
      if(td%pesv%calc_spm .or. td%pesv%calc_mask .or. td%pesv%calc_flux) &
        call pes_calc(td%pesv, gr%mesh, st, td%dt, iter, gr, hm)

      call td_write_iter(write_handler, gr, st, hm, geo, hm%ep%kick, td%dt, sys%ks, iter)

      ! write down data
      call check_point()

      ! check if debug mode should be enabled or disabled on the fly
      call io_debug_on_the_fly()

      call profiling_out(prof)
      if (stopping) exit

    end do propagation

    if(st%d%pack_states .and. hamiltonian_apply_packed(hm, gr%mesh)) call states_unpack(st)

    call restart_end(restart_dump)
    if (ion_dynamics_ions_move(td%ions) .and. td%recalculate_gs) call restart_end(restart_load)
    call td_write_end(write_handler)
    call end_()

#ifdef HAVE_MPI
    ! wait for all processors to finish
    if(st%parallel_in_states) then
      call MPI_Barrier(st%mpi_grp%comm, mpi_err)
    end if
#endif

    POP_SUB(td_run)

  contains

    subroutine print_header()

      if(.not.cmplxscl) then
        write(message(1), '(a7,1x,a14,a14,a10,a17)') 'Iter ', 'Time ', 'Energy ', 'SC Steps', 'Elapsed Time '
      else
        write(message(1), '(a7,1x,a14,a14,a14,a10,a17)') &
          'Iter ', 'Time ', 'Re(Energy) ','Im(Energy) ', 'SC Steps', 'Elapsed Time '
      end if

      call messages_info(1)
      call messages_print_stress(stdout)

    end subroutine print_header

    ! ---------------------------------------------------------
    subroutine check_point()
      PUSH_SUB(td_run.check_point)

      if(.not. cmplxscl) then
        write(message(1), '(i7,1x,2f14.6,i10,f14.3)') iter, &
          units_from_atomic(units_out%time, iter*td%dt), &
          units_from_atomic(units_out%energy, hm%energy%total + geo%kinetic_energy), &
          scsteps, loct_clock() - etime
      else
        write(message(1), '(i7,1x,3f14.6,i10,f14.3)') iter, &
          units_from_atomic(units_out%time, iter*td%dt), &
          units_from_atomic(units_out%energy, hm%energy%total + geo%kinetic_energy), &
          units_from_atomic(units_out%energy, hm%energy%Imtotal), &
          scsteps, loct_clock() - etime
      end if

      call messages_info(1)
      etime = loct_clock()

      if((sys%outp%output_interval > 0 .and. mod(iter, sys%outp%output_interval) == 0) .or. &
        iter == td%max_iter .or. stopping) then ! output
        ! TODO this now overwrites wf inside st. If this is not wanted need to add an optional overwrite=no flag
        if (st%modelmbparticles%nparticle > 0) then
          call modelmb_sym_all_states (gr, st, geo)
        end if
        call td_write_output(write_handler, gr, st, hm, sys%ks, sys%outp, geo, iter, td%dt)
      end if

      if (mod(iter, sys%outp%restart_write_interval) == 0 .or. iter == td%max_iter .or. stopping) then ! restart
        !if(iter == td%max_iter) sys%outp%iter = ii - 1
        call td_write_data(write_handler, gr, st, hm, sys%ks, sys%outp, geo, iter, td%dt)
        call td_dump(restart_dump, gr, st, hm, td, iter, ierr)
        if (ierr /= 0) then
          message(1) = "Unable to write time-dependent restart information."
          call messages_warning(1)
        end if

        call pes_output(td%pesv, gr%mesh, st, iter, sys%outp, td%dt, gr, geo)

        if (ion_dynamics_ions_move(td%ions) .and. td%recalculate_gs) then
          call messages_print_stress(stdout, 'Recalculating the ground state.')
          fromScratch = .false.
          call ground_state_run(sys, hm, fromScratch)
          call states_load(restart_load, st, gr, ierr, iter=iter)
          if (ierr /= 0) then
            message(1) = "Unable to load TD states."
            call messages_fatal(1)
          end if
          call messages_print_stress(stdout, "Time-dependent simulation proceeds")
          call print_header()
        end if
      end if

      POP_SUB(td_run.check_point)
    end subroutine check_point

    ! ---------------------------------------------------------
    subroutine end_()
      PUSH_SUB(td_run.end_)

      ! free memory
      call states_deallocate_wfns(st)
      call ion_dynamics_end(td%ions)
      call td_end(td)
      if (ion_dynamics_ions_move(td%ions) .and. td%recalculate_gs) call restart_end(restart_load)

      POP_SUB(td_run.end_)
    end subroutine end_

    ! ---------------------------------------------------------
    subroutine init_wfs()

      integer :: ierr, freeze_orbitals
      FLOAT :: x
      logical :: freeze_hxc
      type(restart_t) :: restart

      PUSH_SUB(td_run.init_wfs)

      if (.not. fromscratch) then
        call restart_init(restart, RESTART_TD, RESTART_TYPE_LOAD, st%dom_st_kpt_mpi_grp, ierr, mesh=gr%mesh)
        if(ierr == 0) &
          call td_load(restart, gr, st, hm, td, ierr)
        if(ierr /= 0) then
          fromScratch = .true.
          td%iter = 0
          message(1) = "Unable to read time-dependent restart information: Starting from scratch"
          call messages_warning(1)
        end if
        call restart_end(restart)
      end if

      if (td%iter >= td%max_iter) then
        message(1) = "All requested iterations have already been done. Use FromScratch = yes if you want to redo them."
        call messages_info(1)
        POP_SUB(td_run.init_wfs)
        return
      end if

      if (fromScratch) then
        call restart_init(restart, RESTART_GS, RESTART_TYPE_LOAD, st%dom_st_kpt_mpi_grp, ierr, mesh=gr%mesh, exact=.true.)

        if(.not. st%only_userdef_istates) then
          if(ierr == 0) call states_load(restart, st, gr, ierr, label = ": gs")
          if (ierr /= 0) then
            message(1) = 'Unable to read ground-state wavefunctions.'
            call messages_fatal(1)
          end if
        end if

        ! check if we should deploy user-defined wavefunctions.
        ! according to the settings in the input file the routine
        ! overwrites orbitals that were read from restart/gs
        if(parse_is_defined('UserDefinedStates')) call states_read_user_def_orbitals(gr%mesh, st)

        call transform_states(st, restart, gr)
        call restart_end(restart)
      end if



      !%Variable TDFreezeOrbitals
      !%Type integer
      !%Default 0
      !%Section Time-Dependent
      !%Description
      !% (Experimental) You have the possibility of "freezing" a number of orbitals during a time-propagation.
      !% The Hartree and exchange-correlation potential due to these orbitals (which
      !% will be the lowest-energy ones) will be added during the propagation, but the orbitals
      !% will not be propagated.
      !%Option sae -1
      !% Single-active-electron approximation. This option is only valid for time-dependent
      !% calculations (<tt>CalculationMode = td</tt>). Also, the nuclei should not move.
      !% The idea is that all orbitals except the last one are frozen. The orbitals are to
      !% be read from a previous ground-state calculation. The active orbital is then treated
      !% as independent (whether it contains one electron or two) -- although it will
      !% feel the Hartree and exchange-correlation potentials from the ground-state electronic
      !% configuration.
      !%
      !% It is almost equivalent to setting <tt>TDFreezeOrbitals = N-1</tt>, where <tt>N</tt> is the number
      !% of orbitals, but not completely.
      !%End
      call parse_variable('TDFreezeOrbitals', 0, freeze_orbitals)

      if(freeze_orbitals /= 0) call messages_experimental('TDFreezeOrbitals')

      if(.not. cmplxscl) then
        call density_calc(st, gr, st%rho)
      else
        call density_calc(st, gr, st%zrho%Re, st%zrho%Im)
      end if

      if(freeze_orbitals > 0) then
        ! In this case, we first freeze the orbitals, then calculate the Hxc potential.
        call states_freeze_orbitals(st, gr, sys%mc, freeze_orbitals)
        write(message(1),'(a,i4,a,i4,a)') 'Info: The lowest', freeze_orbitals, &
          ' orbitals have been frozen.', st%nst, ' will be propagated.'
        call messages_info(1)
        call v_ks_calc(sys%ks, hm, st, sys%geo, calc_eigenval=.true., time = td%iter*td%dt)
      elseif(freeze_orbitals < 0) then
        ! This means SAE approximation. We calculate the Hxc first, then freeze all
        ! orbitals minus one.
        write(message(1),'(a)') 'Info: The single-active-electron approximation will be used.'
        call messages_info(1)
        call v_ks_calc(sys%ks, hm, st, sys%geo, calc_eigenval=.true., time = td%iter*td%dt)
        call states_freeze_orbitals(st, gr, sys%mc, n = st%nst-1)
        call v_ks_freeze_hxc(sys%ks)
        call density_calc(st, gr, st%rho)
      else
        ! Normal run.
        call v_ks_calc(sys%ks, hm, st, sys%geo, calc_eigenval=.true., time = td%iter*td%dt)
      end if

      !%Variable TDFreezeHXC
      !%Type logical
      !%Default no
      !%Section Time-Dependent
      !%Description
      !% The electrons are evolved as independent particles feeling the Hartree and 
      !% exchange-correlation potentials from the ground-state electronic configuration.
      !%End
      call parse_variable('TDFreezeHXC', .false., freeze_hxc)
      if(freeze_hxc) then 
        write(message(1),'(a)') 'Info: Freezing Hartree and exchange-correlation potentials.'
        call messages_info(1)
        call v_ks_freeze_hxc(sys%ks)
      end if

      x = minval(st%eigenval(st%st_start, :))
#ifdef HAVE_MPI
      if(st%parallel_in_states) then
        call MPI_Bcast(x, 1, MPI_FLOAT, 0, st%mpi_grp%comm, mpi_err)
      end if
#endif
      call hamiltonian_span(hm, minval(gr%mesh%spacing(1:gr%mesh%sb%dim)), x)
      ! initialize Fermi energy
      call states_fermi(st, gr%mesh)
      call energy_calc_total(hm, gr, st)

      POP_SUB(td_run.init_wfs)
    end subroutine init_wfs


    ! ---------------------------------------------------------
    subroutine td_run_zero_iter()
      PUSH_SUB(td_run.td_run_zero_iter)

      call td_write_iter(write_handler, gr, st, hm, geo, hm%ep%kick, td%dt, sys%ks, 0)

      ! I apply the delta electric field *after* td_write_iter, otherwise the
      ! dipole matrix elements in write_proj are wrong
      if(hm%ep%kick%time  ==  M_ZERO) then
        if(.not. cmplxscl) then
          call kick_apply(gr%mesh, st, td%ions, geo, hm%ep%kick)
        else
          call kick_apply(gr%mesh, st, td%ions, geo, hm%ep%kick, hm%cmplxscl%theta)
        end if
        call td_write_kick(gr%mesh, hm%ep%kick, sys%outp, geo, 0)
      end if
      call propagator_run_zero_iter(hm, gr, td%tr)
      if (sys%outp%output_interval > 0) then
        call td_write_data(write_handler, gr, st, hm, sys%ks, sys%outp, geo, 0)
        call td_write_output(write_handler, gr, st, hm, sys%ks, sys%outp, geo, 0)
      end if

      POP_SUB(td_run.td_run_zero_iter)
    end subroutine td_run_zero_iter


    ! ---------------------------------------------------------
    !> reads the pos and vel from coordinates file
    subroutine td_read_coordinates() 
      integer :: iatom, iter, iunit
      PUSH_SUB(td_run.td_read_coordinates)

      call io_assign(iunit)
      open(unit = iunit, file = io_workpath('td.general/coordinates'), action='read', status='old')

      if(iunit < 0) then
        message(1) = "Could not open file '"//trim(io_workpath('td.general/coordinates'))//"'."
        message(2) = "Starting simulation from initial geometry."
        call messages_warning(2)
        POP_SUB(td_run.td_read_coordinates)
        return
      end if

      call io_skip_header(iunit)
      do iter = 0, td%iter - 1
        read(iunit, *) ! skip previous iterations... sorry, but no portable seek in Fortran
      end do
      read(iunit, '(28x)', advance='no') ! skip the time index.

      do iatom = 1, geo%natoms
        read(iunit, '(3es20.12)', advance='no') geo%atom(iatom)%x(1:gr%mesh%sb%dim)
        geo%atom(iatom)%x(:) = units_to_atomic(units_inp%length, geo%atom(iatom)%x(:))
      end do
      do iatom = 1, geo%natoms
        read(iunit, '(3es20.12)', advance='no') geo%atom(iatom)%v(1:gr%mesh%sb%dim)
        geo%atom(iatom)%v(:) = units_to_atomic(units_inp%velocity, geo%atom(iatom)%v(:))
      end do
      do iatom = 1, geo%natoms
        read(iunit, '(3es20.12)', advance='no') geo%atom(iatom)%f(1:gr%mesh%sb%dim)
        geo%atom(iatom)%f(:) = units_to_atomic(units_inp%force, geo%atom(iatom)%f(:))
      end do

      call io_close(iunit)
      POP_SUB(td_run.td_read_coordinates)
    end subroutine td_read_coordinates

  end subroutine td_run


  ! ---------------------------------------------------------
  subroutine transform_states(st, restart, gr, prefix)
    type(states_t),             intent(inout) :: st
    type(restart_t),            intent(inout) :: restart
    type(grid_t),               intent(in)    :: gr
    character(len=*), optional, intent(in)    :: prefix

    type(states_t) :: stin
    type(block_t) :: blk
    CMPLX, allocatable :: rotation_matrix(:,:), psi(:, :)
    integer :: ist, jst, ncols, iqn
    character(len=256) :: block_name
    
    PUSH_SUB(transform_states)

    block_name = trim(optional_default(prefix, "")) // "TransformStates"
    
    !%Variable TransformStates
    !%Type block
    !%Default no
    !%Section States
    !%Description
    !% Before starting the <tt>td</tt> calculation, the initial states (that are
    !% read from the <tt>restart/gs</tt> directory, which should have been
    !% generated in a previous ground-state calculation) can be "transformed"
    !% among themselves. The block <tt>TransformStates</tt> gives the transformation matrix
    !% to be used. The number of rows and columns of the matrix should equal the number
    !% of the states present in the time-dependent calculation (the independent
    !% spin and <i>k</i>-point subspaces are all transformed equally); the number of
    !% columns should be equal to the number of states present in the
    !% <tt>restart/gs</tt> directory. This number may be different: for example,
    !% one could have run previously in <tt>unocc</tt> mode in order to obtain unoccupied
    !% Kohn-Sham states, and therefore <tt>restart/gs</tt> will contain more states.
    !% These states can be used in the transformation. 
    !%
    !% Note that the code will not check the orthonormality of the new states!
    !%
    !% Each line provides the coefficients of the new states, in terms of
    !% the old ones. The coefficients are complex, but the imaginary part will be
    !% ignored for real wavefunctions.
    !% Note: This variable cannot be used when parallel in states.
    !%End
    if(parse_is_defined(trim(block_name))) then
      if(parse_block(trim(block_name), blk) == 0) then
        if(st%parallel_in_states) &
          call messages_not_implemented(trim(block_name) // " parallel in states")
        if(parse_block_n(blk) /= st%nst) then
          message(1) = "Number of rows in block " // trim(block_name) // " must equal number of states in this calculation."
          call messages_fatal(1)
        end if
        call states_copy(stin, st, exclude_wfns = .true.)
        call states_look_and_load(restart, stin, gr)

        ! FIXME: rotation matrix should be R_TYPE
        SAFE_ALLOCATE(rotation_matrix(1:stin%nst, 1:stin%nst))
        SAFE_ALLOCATE(psi(1:gr%mesh%np, 1:st%d%dim))
        
        rotation_matrix = M_z0
        forall(ist = 1:stin%nst) rotation_matrix(ist, ist) = CNST(1.0)
        
        do ist = 1, st%nst
          ncols = parse_block_cols(blk, ist-1)
          if(ncols /= stin%nst) then            
            write(message(1),'(a,i6,a,i6,3a,i6,a)') "Number of columns (", ncols, ") in row ", ist, " of block ", &
              trim(block_name), " must equal number of states (", stin%nst, ") read from gs restart."
            call messages_fatal(1)
          end if
          do jst = 1, stin%nst
            call parse_block_cmplx(blk, ist - 1, jst - 1, rotation_matrix(jst, ist))
          end do
        end do

        call parse_block_end(blk)

        do iqn = st%d%kpt%start, st%d%kpt%end
          if(states_are_real(st)) then
            call states_rotate(gr%mesh, stin, real(rotation_matrix, REAL_PRECISION), iqn)
          else
            call states_rotate(gr%mesh, stin, rotation_matrix, iqn)
          end if

          do ist = st%st_start, st%st_end 
            call states_get_state(stin, gr%mesh, ist, iqn, psi)
            call states_set_state(st, gr%mesh, ist, iqn, psi)
          end do

        end do

        SAFE_DEALLOCATE_A(rotation_matrix)
        SAFE_DEALLOCATE_A(psi)

        call states_end(stin)

      else
        call messages_input_error(trim(block_name), '"' // trim(block_name) // '" has to be specified as block.')
      end if
      
    end if

    POP_SUB(transform_states)
  end subroutine transform_states

  ! ---------------------------------------------------------
  subroutine td_dump(restart, gr, st, hm, td, iter, ierr)
    type(restart_t),     intent(in)  :: restart
    type(grid_t),        intent(in)  :: gr
    type(states_t),      intent(in)  :: st
    type(hamiltonian_t), intent(in)  :: hm
    type(td_t),          intent(in)  :: td
    integer,             intent(in)  :: iter
    integer,             intent(out) :: ierr

    logical :: cmplxscl
    integer :: err, err2

    PUSH_SUB(td_dump)

    ierr = 0

    if (restart_skip(restart)) then
      POP_SUB(td_dump)
      return
    end if

    if (debug%info) then
      message(1) = "Debug: Writing td restart."
      call messages_info(1)
    end if

    ! first write resume file
    call states_dump(restart, st, gr, err, iter=iter)
    if (err /= 0) ierr = ierr + 1

    cmplxscl = st%cmplxscl%space      
    call potential_interpolation_dump(td%tr%vksold, restart, gr, cmplxscl, st%d%nspin, err2)
    if (err2 /= 0) ierr = ierr + 2

    call pes_dump(restart, td%pesv, st, gr%mesh, err)
    if (err /= 0) ierr = ierr + 4

    ! Gauge field restart
    if (gauge_field_is_applied(hm%ep%gfield)) then
      call gauge_field_dump(restart, hm%ep%gfield, ierr)
    end if

    if (debug%info) then
      message(1) = "Debug: Writing td restart done."
      call messages_info(1)
    end if

    POP_SUB(td_dump)
  end subroutine td_dump

  ! ---------------------------------------------------------
  subroutine td_load(restart, gr, st, hm, td, ierr)
    type(restart_t),     intent(in)    :: restart
    type(grid_t),        intent(in)    :: gr
    type(states_t),      intent(inout) :: st
    type(hamiltonian_t), intent(inout) :: hm
    type(td_t),          intent(inout) :: td
    integer,             intent(out)   :: ierr

    logical :: cmplxscl
    integer :: err, err2
    PUSH_SUB(td_load)

    ierr = 0

    if (restart_skip(restart)) then
      ierr = -1
      POP_SUB(td_load)
      return
    end if

    if (debug%info) then
      message(1) = "Debug: Reading td restart."
      call messages_info(1)
    end if

    ! Read states
    call states_load(restart, st, gr, err, iter=td%iter, read_left = st%have_left_states, label = ": td")
    if (err /= 0) then
      ierr = ierr + 1
    end if

    ! read potential from previous interactions
    cmplxscl = st%cmplxscl%space
    call potential_interpolation_load(td%tr%vksold, restart, gr, cmplxscl, st%d%nspin, err2)

    if (err2 /= 0) ierr = ierr + 2

    ! read PES restart
    if (td%pesv%calc_spm .or. td%pesv%calc_mask .or. td%pesv%calc_flux) then
      call pes_load(restart, td%pesv, st, gr%mesh, err)
      if (err /= 0) ierr = ierr + 4
    end if

    ! Gauge field restart
    if (gauge_field_is_applied(hm%ep%gfield)) then
      call gauge_field_load(restart, hm%ep%gfield, err)
      if (err /= 0) then
        ierr = ierr + 8
      else
        call hamiltonian_update(hm, gr%mesh, time = td%dt*td%iter)
      end if
    end if

    if (debug%info) then
      message(1) = "Debug: Reading td restart done."
      call messages_info(1)
    end if

    POP_SUB(td_load)
  end subroutine td_load

end module td_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
