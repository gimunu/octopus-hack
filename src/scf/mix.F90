!! Copyright (C) 2002-2014 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira
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
!! $Id: mix.F90 15523 2016-07-24 23:22:19Z xavier $

#include "global.h"

module mix_oct_m
  use derivatives_oct_m
  use global_oct_m
  use io_oct_m
  use io_function_oct_m
  use lalg_adv_oct_m
  use lalg_basic_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use nl_operator_oct_m
  use parser_oct_m
  use profiling_oct_m
  use restart_oct_m
  use stencil_cube_oct_m
  use types_oct_m
  use unit_system_oct_m
  use varinfo_oct_m

  implicit none

  private
  public ::                     &
    mix_set_mixing,             &
    mix_t,                      &
    mix_init,                   &
    mix_clear,                  &
    mix_end,                    &
    mix_dump,                   &
    mix_load,                   &
    dmixing,                    &
    zmixing,                    &
    mix_coefficient

  type mix_t
    private
    integer :: scheme           !< the mixing scheme used (linear, broyden, etc)

    FLOAT :: coeff              !< the mixing coefficient (in linear mixing: vnew = (1-coeff)*vin + coeff*vout)

    integer :: iter             !< number of SCF iterations already done. In case of restart, this number must
                                !< include the iterations done in previous calculations.

    type(type_t) :: func_type   !< type of the functions to be mixed
    integer :: ns               !< number of steps used to extrapolate the new vector

    integer :: d1, d2, d3, d4   !< the dimensions of the arrays that store the information from the previous iterations

    integer :: last_ipos        !< where is the information about the last iteration stored in arrays df and dv

    integer :: interval

    FLOAT, pointer :: ddf(:, :, :, :)
    FLOAT, pointer :: ddv(:, :, :, :)
    FLOAT, pointer :: df_old(:, :, :)
    FLOAT, pointer :: dvin_old(:, :, :)

    CMPLX, pointer :: zdf(:, :, :, :)
    CMPLX, pointer :: zdv(:, :, :, :)
    CMPLX, pointer :: zf_old(:, :, :)
    CMPLX, pointer :: zvin_old(:, :, :)

    type(derivatives_t), pointer :: der
    logical                      :: precondition
    type(nl_operator_t)          :: preconditioner

    FLOAT :: residual_coeff
    
  end type mix_t

contains

  ! ---------------------------------------------------------
  subroutine mix_init(smix, der, d1, d2, d3, def_, func_type, prefix_)
    type(mix_t),                   intent(out) :: smix
    type(derivatives_t), target,   intent(in)  :: der
    integer,                       intent(in)  :: d1, d2, d3
    integer,             optional, intent(in)  :: def_
    type(type_t),        optional, intent(in)  :: func_type
    character(len=*),    optional, intent(in)  :: prefix_

    integer :: def
    character(len=32) :: prefix

    PUSH_SUB(mix_init)

    smix%der => der

    def = OPTION__MIXINGSCHEME__BROYDEN
    if(present(def_)) def = def_
    if(present(func_type)) then 
      smix%func_type = func_type
    else 
      smix%func_type = TYPE_FLOAT
    end if
    prefix = ""
    if(present(prefix_)) prefix = prefix_

    call messages_obsolete_variable('TypeOfMixing', 'MixingScheme')
    
    !%Variable MixingScheme
    !%Type integer
    !%Default broyden
    !%Section SCF::Mixing
    !%Description
    !% The scheme used to produce, at each iteration in the self-consistent cycle
    !% that attempts to solve the Kohn-Sham equations, the input density from the value
    !% of the input and output densities of previous iterations.
    !%Option linear 0
    !% Simple linear mixing.
    !%Option broyden 2
    !% Broyden scheme [C. G Broyden, <i>Math. Comp.</i> <b>19</b>, 577 (1965); 
    !% D. D. Johnson, <i>Phys. Rev. B</i> <b>38</b>, 12807 (1988)].
    !% For complex functions (e.g. Sternheimer with <tt>EMEta</tt> > 0), we use the generalization
    !% with a complex dot product.
    !%Option diis 9
    !% Direct inversion in the iterative subspace (diis)
    !% scheme [P. Pulay, <i>Chem. Phys. Lett.</i>, <b>73</b>, 393
    !% (1980)] as described in [G. Kresse, and J. Hurthmueller,
    !% <i>Phys. Rev. B</i> <b>54</b>, 11169 (1996)].
    !%Option bowler_gillan 1
    !% The Guaranteed-reduction modification of the Pulay scheme by
    !% Bowler and Gillan [D. R. Bowler and M. J. Gillan,
    !% <i>Chem. Phys.  Lett.</i> <b>325</b>, 473 (2000)].
    !%End
    call parse_variable(trim(prefix)//'MixingScheme', def, smix%scheme)
    if(.not.varinfo_valid_option('MixingScheme', smix%scheme)) call messages_input_error('MixingScheme', 'invalid option')
    call messages_print_var_option(stdout, "MixingScheme", smix%scheme)

    if(smix%scheme == OPTION__MIXINGSCHEME__DIIS) call messages_experimental('MixingScheme = diis')

    !%Variable MixingPreconditioner
    !%Type logical
    !%Default false
    !%Section SCF::Mixing
    !%Description
    !% (Experimental) If set to yes, Octopus will use a preconditioner
    !% for the mixing operator.
    !%End
    call parse_variable(trim(prefix)+'MixingPreconditioner', .false., smix%precondition)
    if(smix%precondition) call messages_experimental('MixingPreconditioner')
    
    !%Variable Mixing
    !%Type float
    !%Default 0.3
    !%Section SCF::Mixing
    !%Description
    !% The linear, Broyden and DIIS scheme depend on a "mixing parameter", set by this variable. 
    !% Must be 0 < <tt>Mixing</tt> <= 1.
    !%End
    call parse_variable(trim(prefix)+'Mixing', CNST(0.3), smix%coeff)
    if(smix%coeff <= M_ZERO .or. smix%coeff > M_ONE) then
      call messages_input_error('Mixing', 'Value should be positive and smaller than one.')
    end if
    
    !%Variable MixingResidual
    !%Type float
    !%Default 0.05
    !%Section SCF::Mixing
    !%Description
    !% In the DIIS mixing it is benefitial to include a bit of
    !% residual into the mixing. This parameter controls this amount.
    !%End
    call parse_variable(trim(prefix)+'MixingResidual', CNST(0.05), smix%residual_coeff)
    if(smix%residual_coeff <= M_ZERO .or. smix%residual_coeff > M_ONE) then
      call messages_input_error('MixingResidual', 'Value should be positive and smaller than one.')
    end if
    
    !%Variable MixNumberSteps
    !%Type integer
    !%Default 3
    !%Section SCF::Mixing
    !%Description
    !% In the Broyden and Bowler_Gillan schemes, the new input density or potential is constructed
    !% from the values of the densities/potentials of a given number of previous iterations.
    !% This number is set by this variable. Must be greater than 1.
    !%End
    if (smix%scheme /= OPTION__MIXINGSCHEME__LINEAR) then
      call parse_variable(trim(prefix)//'MixNumberSteps', 3, smix%ns)
      if(smix%ns <= 1) call messages_input_error('MixNumberSteps')
    else
      smix%ns = 0
    end if
    
    !%Variable MixInterval
    !%Type integer
    !%Default 1
    !%Section SCF::Mixing
    !%Description
    !% When this variable is set to a value different than 1 (the
    !% defaul) a combined mixing scheme will be used, with MixInterval
    !% - 1 steps of linear mixing followed by 1 step of the selected
    !% mixing. For the moment this variable only works with DIIS mixing.
    !%End
    call parse_variable(trim(prefix)//'MixInterval', 1, smix%interval)
    if(smix%interval < 1) call messages_input_error('MixInterval', 'MixInterval must be larger or equal than 1')
    
    smix%iter = 0

    nullify(smix%ddf)
    nullify(smix%ddv)
    nullify(smix%df_old)
    nullify(smix%dvin_old)

    nullify(smix%zdf)
    nullify(smix%zdv)
    nullify(smix%zf_old)
    nullify(smix%zvin_old)

    smix%d1 = d1
    smix%d2 = d2
    smix%d3 = d3
    select case (smix%scheme)
    case (OPTION__MIXINGSCHEME__LINEAR, OPTION__MIXINGSCHEME__BROYDEN, OPTION__MIXINGSCHEME__DIIS)
      smix%d4 = smix%ns
    case (OPTION__MIXINGSCHEME__BOWLER_GILLAN)
      smix%d4 = smix%ns + 1
    end select

    if (smix%scheme /= OPTION__MIXINGSCHEME__LINEAR) then
      if(smix%func_type == TYPE_FLOAT) then 
        SAFE_ALLOCATE(     smix%ddf(1:d1, 1:d2, 1:d3, 1:smix%d4))
        SAFE_ALLOCATE(smix%dvin_old(1:d1, 1:d2, 1:d3))
        SAFE_ALLOCATE(     smix%ddv(1:d1, 1:d2, 1:d3, 1:smix%d4))
        SAFE_ALLOCATE(  smix%df_old(1:d1, 1:d2, 1:d3))
      else
        SAFE_ALLOCATE(     smix%zdf(1:d1, 1:d2, 1:d3, 1:smix%d4))
        SAFE_ALLOCATE(smix%zvin_old(1:d1, 1:d2, 1:d3))
        SAFE_ALLOCATE(     smix%zdv(1:d1, 1:d2, 1:d3, 1:smix%d4))
        SAFE_ALLOCATE(  smix%zf_old(1:d1, 1:d2, 1:d3))
      end if
    end if

    call mix_clear(smix)

    if(smix%precondition) call init_preconditioner()

    POP_SUB(mix_init)

  contains

    subroutine init_preconditioner()

      integer :: ns, maxp, ip, is
      FLOAT, parameter :: weight = 50.0
      
      ! This the mixing preconditioner from GPAW:
      !
      !   https://wiki.fysik.dtu.dk/gpaw/documentation/densitymix/densitymix.html
      !

      ASSERT(.not. der%mesh%use_curvilinear)
      ASSERT(der%mesh%sb%dim == 3)
      
      call nl_operator_init(smix%preconditioner, "Mixing preconditioner")
      call stencil_cube_get_lapl(smix%preconditioner%stencil, der%mesh%sb%dim, 1)
      call nl_operator_build(der%mesh, smix%preconditioner, der%mesh%np, const_w = .not. der%mesh%use_curvilinear)
      
      ns = smix%preconditioner%stencil%size

      if (smix%preconditioner%const_w) then
        maxp = 1
      else
        maxp = der%mesh%np
      end if

      do ip = 1, maxp

        do is = 1, ns
          select case(sum(abs(smix%preconditioner%stencil%points(1:der%mesh%sb%dim, is))))
          case(0)
            smix%preconditioner%w_re(is, ip) = CNST(1.0) + weight/CNST(8.0)
          case(1)
            smix%preconditioner%w_re(is, ip) = weight/CNST(16.0)
          case(2)
            smix%preconditioner%w_re(is, ip) = weight/CNST(32.0)
          case(3)
            smix%preconditioner%w_re(is, ip) = weight/CNST(64.0)
          case default
            ASSERT(.false.)
          end select

        end do
      end do
      
      call nl_operator_update_weights(smix%preconditioner)

    end subroutine init_preconditioner

  end subroutine mix_init


  ! ---------------------------------------------------------
  subroutine mix_clear(smix)
    type(mix_t),             intent(inout) :: smix
    
    PUSH_SUB(mix_clear)

    if (smix%scheme /= OPTION__MIXINGSCHEME__LINEAR) then
      if(smix%func_type == TYPE_FLOAT) then 
        smix%ddf = M_ZERO
        smix%ddv = M_ZERO
        smix%dvin_old = M_ZERO
        smix%df_old = M_ZERO
      else
        smix%zdf = M_z0
        smix%zdv = M_z0
        smix%zvin_old = M_z0
        smix%zf_old = M_z0
      end if
    end if

    smix%iter = 0
    smix%last_ipos = 0

    POP_SUB(mix_clear)
  end subroutine mix_clear


  ! ---------------------------------------------------------
  subroutine mix_end(smix)
    type(mix_t), intent(inout) :: smix

    PUSH_SUB(mix_end)

    if(smix%precondition) call nl_operator_end(smix%preconditioner)

    ! Arrays got allocated for all mixing schemes, except linear mixing
    if (smix%scheme /= OPTION__MIXINGSCHEME__LINEAR) then
      SAFE_DEALLOCATE_P(smix%ddf)
      SAFE_DEALLOCATE_P(smix%ddv)
      SAFE_DEALLOCATE_P(smix%dvin_old)
      SAFE_DEALLOCATE_P(smix%df_old)

      SAFE_DEALLOCATE_P(smix%zdf)
      SAFE_DEALLOCATE_P(smix%zdv)
      SAFE_DEALLOCATE_P(smix%zvin_old)
      SAFE_DEALLOCATE_P(smix%zf_old)
    end if

    POP_SUB(mix_end)
  end subroutine mix_end


  ! ---------------------------------------------------------
  subroutine mix_set_mixing(smix, newmixing)
    type(mix_t), intent(inout) :: smix
    FLOAT, intent(in):: newmixing

    PUSH_SUB(mix_set_mixing)
    
    if(smix%scheme == OPTION__MIXINGSCHEME__LINEAR) then
      smix%coeff = newmixing
    else
    !  message(1) = "Mixing can only be adjusted in linear mixing scheme."
    !  call messages_fatal(1)
    end if
    
    POP_SUB(mix_set_mixing)
  end subroutine mix_set_mixing


  ! ---------------------------------------------------------
  subroutine mix_dump(restart, smix, mesh, ierr)
    type(restart_t), intent(in)  :: restart
    type(mix_t),     intent(in)  :: smix
    type(mesh_t),    intent(in)  :: mesh
    integer,         intent(out) :: ierr

    integer :: iunit, id2, id3, id4, err, err2(4)
    character(len=40) :: lines(8)
    character(len=80) :: filename

    PUSH_SUB(mix_dump)

    ierr = 0

    if (restart_skip(restart)) then
      POP_SUB(mix_dump)
      return
    end if

    if (debug%info) then
      message(1) = "Debug: Writing mixing restart."
      call messages_info(1)
    end if

    ! functions to be written need to be compatible with the mesh
    ASSERT(mesh%np == smix%d1)

    ! First we write some information about the mixing
    iunit = restart_open(restart, 'mixing')
    write(lines(1), '(a11,i1)')  'scheme=    ', smix%scheme
    ! Number of global mesh points have to be written, not only smix%d1
    write(lines(2), '(a11,i10)') 'd1=        ', mesh%np_global
    write(lines(3), '(a11,i10)') 'd2=        ', smix%d2
    write(lines(4), '(a11,i10)') 'd3=        ', smix%d3
    write(lines(5), '(a11,i10)') 'd4=        ', smix%d4
    write(lines(6), '(a11,i10)') 'iter=      ', smix%iter
    write(lines(7), '(a11,i10)') 'ns=        ', smix%ns
    write(lines(8), '(a11,i10)') 'last_ipos= ', smix%last_ipos
    call restart_write(restart, iunit, lines, 8, err)
    if (err /= 0) ierr = ierr + 1
    call restart_close(restart, iunit)

    ! Now we write the different functions. 
    ! These are not needed when using linear mixing, so we will make sure we skip this step in that case.
    err2 = 0
    if (smix%scheme /= OPTION__MIXINGSCHEME__LINEAR) then
      do id2 = 1, smix%d2
        do id3 = 1, smix%d3
          do id4 = 1, smix%d4

            write(filename,'(a3,i2.2,i2.2,i2.2)') 'df_', id2, id3, id4
            if (smix%func_type == TYPE_FLOAT) then
              call drestart_write_mesh_function(restart, filename, mesh, smix%ddf(1:mesh%np, id2, id3, id4), err)
            else
              call zrestart_write_mesh_function(restart, filename, mesh, smix%zdf(1:mesh%np, id2, id3, id4), err)
            end if
            if (err /= 0) err2(1) = err2(1) + 1

            write(filename,'(a3,i2.2,i2.2,i2.2)') 'dv_', id2, id3, id4
            if (smix%func_type == TYPE_FLOAT) then
              call drestart_write_mesh_function(restart, filename, mesh, smix%ddv(1:mesh%np, id2, id3, id4), err)
            else
              call zrestart_write_mesh_function(restart, filename, mesh, smix%zdv(1:mesh%np, id2, id3, id4), err)
            end if
            if (err /= 0) err2(2) = err2(2) + 1
              
          end do

          write(filename,'(a6,i2.2,i2.2)') 'f_old_', id2, id3
          if (smix%func_type == TYPE_FLOAT) then
            call drestart_write_mesh_function(restart, filename, mesh, smix%df_old(1:mesh%np, id2, id3), err)
          else
            call zrestart_write_mesh_function(restart, filename, mesh, smix%zf_old(1:mesh%np, id2, id3), err)
          end if
          if (err /= 0) err2(3) = err2(3) + 1

          write(filename,'(a8,i2.2,i2.2)') 'vin_old_', id2, id3
          if (smix%func_type == TYPE_FLOAT) then
            call drestart_write_mesh_function(restart, filename, mesh, smix%dvin_old(1:mesh%np, id2, id3), err)
          else
            call zrestart_write_mesh_function(restart, filename, mesh, smix%zvin_old(1:mesh%np, id2, id3), err)
          end if
          if (err /= 0) err2(4) = err2(4) + 1
          
        end do
      end do

      if (err2(1) /= 0) ierr = ierr + 2
      if (err2(2) /= 0) ierr = ierr + 4
      if (err2(3) /= 0) ierr = ierr + 8
      if (err2(4) /= 0) ierr = ierr + 16
    end if

    if (debug%info) then
      message(1) = "Debug: Writing mixing restart done."
      call messages_info(1)
    end if

    POP_SUB(mix_dump)
  end subroutine mix_dump


  !---------------------------------------------------------
  subroutine mix_load(restart, smix, mesh, ierr)
    type(restart_t), intent(in)    :: restart
    type(mix_t),     intent(inout) :: smix
    type(mesh_t),    intent(in)    :: mesh
    integer,         intent(out)   :: ierr

    integer :: iunit, err, err2(4)
    integer :: scheme, d1, d2, d3, d4, ns
    integer :: id2, id3, id4
    character(len=11)  :: str
    character(len=80)  :: filename
    character(len=256) :: lines(8)

    PUSH_SUB(mix_load)

    ierr = 0

    if (restart_skip(restart)) then
      ierr = -1
      POP_SUB(mix_load)
      return
    end if

    if (debug%info) then
      message(1) = "Debug: Reading mixing restart."
      call messages_info(1)
    end if

    ! First we read some information about the mixing
    iunit = restart_open(restart, 'mixing')
    call restart_read(restart, iunit, lines, 8, err)
    if (err /= 0) then
      ierr = ierr + 1
    else
      read(lines(1), *) str, scheme
      read(lines(2), *) str, d1
      read(lines(3), *) str, d2
      read(lines(4), *) str, d3
      read(lines(5), *) str, d4
      read(lines(6), *) str, smix%iter
      read(lines(7), *) str, ns
      read(lines(8), *) str, smix%last_ipos
    end if
    call restart_close(restart, iunit)


    if (ierr == 0) then
      ! We can only use the restart information if the mixing scheme and the number of steps used remained the same
      if (scheme /= smix%scheme .or. ns /= smix%ns) then
        message(1) = "The mixing scheme from the restart data is not the same as the one used in the current calculation."
        call messages_warning(1)
        ierr = ierr + 2
      end if

      ! Check the dimensions of the arrays to be read
      if (mesh%np_global /= d1 .or. mesh%np /= smix%d1 .or. d2 /= smix%d2 .or. d3 /= smix%d3 ) then
        message(1) = "The dimensions of the arrays from the mixing restart data"
        message(2) = "are not the same as the ones used in this calculation."
        call messages_warning(2)
        ierr = ierr + 4
      end if
    end if


    ! Now we read the different functions.
    ! Note that we may have more or less functions than the ones needed (d4 /= smix%d4)
    if (ierr == 0) then
      if (smix%scheme /= OPTION__MIXINGSCHEME__LINEAR) then
        err2 = 0
        do id2 = 1, smix%d2
          do id3 = 1, smix%d3
            do id4 = 1, smix%d4

              write(filename,'(a3,i2.2,i2.2,i2.2)') 'df_', id2, id3, id4
              if (smix%func_type == TYPE_FLOAT) then
                call drestart_read_mesh_function(restart, trim(filename), mesh, smix%ddf(1:mesh%np, id2, id3, id4), err)
              else
                call zrestart_read_mesh_function(restart, trim(filename), mesh, smix%zdf(1:mesh%np, id2, id3, id4), err)
              end if
              if (err /= 0) err2(1) = err2(1) + 1
          
              write(filename,'(a3,i2.2,i2.2,i2.2)') 'dv_', id2, id3, id4
              if (smix%func_type == TYPE_FLOAT) then
                call drestart_read_mesh_function(restart, trim(filename), mesh, smix%ddv(1:mesh%np, id2, id3, id4), err)
              else
                call zrestart_read_mesh_function(restart, trim(filename), mesh, smix%zdv(1:mesh%np, id2, id3, id4), err)
              end if
              if (err /= 0) err2(2) = err2(2) + 1

            end do

            write(filename,'(a6,i2.2,i2.2)') 'f_old_', id2, id3
            if (smix%func_type == TYPE_FLOAT) then
              call drestart_read_mesh_function(restart, trim(filename), mesh, smix%df_old(1:mesh%np, id2, id3), err)
            else
              call zrestart_read_mesh_function(restart, trim(filename), mesh, smix%zf_old(1:mesh%np, id2, id3), err)
            end if
            if (err /= 0) err2(3) = err2(3) + 1

            write(filename,'(a8,i2.2,i2.2)') 'vin_old_', id2, id3
            if (smix%func_type == TYPE_FLOAT) then
              call drestart_read_mesh_function(restart, trim(filename), mesh, smix%dvin_old(1:mesh%np, id2, id3), err)
            else
              call zrestart_read_mesh_function(restart, trim(filename), mesh, smix%zvin_old(1:mesh%np, id2, id3), err)
            end if
            if (err /= 0) err2(4) = err2(4) + 1

          end do
        end do

        if (err2(1) /= 0) ierr = ierr + 8
        if (err2(2) /= 0) ierr = ierr + 16
        if (err2(3) /= 0) ierr = ierr + 32
        if (err2(4) /= 0) ierr = ierr + 64
      end if
    end if

    if (ierr /= 0) then
      ! Something went wront, so make sure we start from scratch
      call mix_clear(smix)
    end if

    if (debug%info) then
      message(1) = "Debug: Reading mixing restart done."
      call messages_info(1)
    end if

    POP_SUB(mix_load)
  end subroutine mix_load

  FLOAT pure function mix_coefficient(this) result(coefficient)
    type(mix_t), intent(in) :: this
    
    coefficient = this%coeff
  end function mix_coefficient
  
  
#include "undef.F90"
#include "real.F90"

#include "mix_inc.F90"

#include "undef.F90"
#include "complex.F90"

#include "mix_inc.F90"

end module mix_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
