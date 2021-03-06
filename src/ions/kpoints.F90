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
!! $Id: kpoints.F90 15370 2016-05-18 10:43:13Z nicolastd $

#include "global.h"
  
module kpoints_oct_m
  use geometry_oct_m
  use global_oct_m
  use loct_oct_m
  use messages_oct_m
  use parser_oct_m
  use profiling_oct_m
  use sort_oct_m
  use symmetries_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use utils_oct_m
  
  implicit none
  
  private
  
  public ::                       &
    kpoints_grid_t,               &
    kpoints_t,                    &
    kpoints_nullify,              &
    kpoints_init,                 &
    kpoints_end,                  &
    kpoints_copy,                 &
    kpoints_number,               &
    kpoints_get_weight,           &
    kpoints_get_point,            &
    kpoints_write_info,           &
    kpoints_point_is_gamma,       &
    kpoints_get_symmetry_ops,     &
    kpoints_get_num_symmetry_ops, &
    kpoints_kweight_denominator,  &
    kpoints_grid_generate,        &
    kpoints_have_zero_weight_path,&
    kpoints_to_absolute

  type kpoints_grid_t
    FLOAT, pointer :: point(:, :)
    FLOAT, pointer :: red_point(:, :)
    FLOAT, pointer :: weight(:)
    integer          :: nshifts            !< number of shifts
    FLOAT, pointer   :: shifts(:,:)
    integer        :: npoints
    integer        :: dim
  end type kpoints_grid_t

  type kpoints_t
    type(kpoints_grid_t) :: full
    type(kpoints_grid_t) :: reduced

    integer        :: method

    logical        :: use_symmetries
    logical        :: use_time_reversal
    integer        :: nik_skip=0 !< number of user defined points with zero weight

    !> For the modified Monkhorst-Pack scheme
    integer        :: nik_axis(MAX_DIM)    !< number of MP divisions
    integer, pointer :: symmetry_ops(:, :)  !< (reduced%npoints, nops)
    integer, pointer :: num_symmetry_ops(:) !< (reduced%npoints)

    FLOAT, pointer :: klattice(:, :)
  end type kpoints_t

  integer, parameter ::                &
    KPOINTS_GAMMA       =  1,          &
    KPOINTS_MONKH_PACK  =  2,          &
    KPOINTS_USER        =  3

contains

  elemental subroutine  kpoints_grid_nullify(this)
    type(kpoints_grid_t), intent(out) :: this

    nullify(this%point, this%red_point, this%weight, this%shifts)
    this%npoints = 0
    this%dim = 0
    this%nshifts = 1
  
  end subroutine kpoints_grid_nullify

  subroutine kpoints_grid_init(dim, this, npoints, nshifts)
    integer,              intent(in)  :: dim
    type(kpoints_grid_t), intent(out) :: this
    integer,              intent(in)  :: npoints
    integer,              intent(in)  :: nshifts

    PUSH_SUB(kpoints_grid_init)

    this%dim = dim
    this%npoints = npoints
    this%nshifts = nshifts
    SAFE_ALLOCATE(this%red_point(1:dim, 1:npoints))
    SAFE_ALLOCATE(this%point(1:dim, 1:npoints))
    SAFE_ALLOCATE(this%weight(1:npoints))
    SAFE_ALLOCATE(this%shifts(1:dim,1:nshifts)) 

    POP_SUB(kpoints_grid_init)
  end subroutine kpoints_grid_init

  ! ---------------------------------------------------------
  subroutine kpoints_grid_end(this)
    type(kpoints_grid_t), intent(inout) :: this

    PUSH_SUB(kpoints_grid_end)

    SAFE_DEALLOCATE_P(this%red_point)
    SAFE_DEALLOCATE_P(this%point)
    SAFE_DEALLOCATE_P(this%weight)
    SAFE_DEALLOCATE_P(this%shifts)

    POP_SUB(kpoints_grid_end)
  end subroutine kpoints_grid_end


  ! ---------------------------------------------------------
  subroutine kpoints_grid_copy(bb, aa)
    type(kpoints_grid_t), intent(in)  :: bb
    type(kpoints_grid_t), intent(out) :: aa

    PUSH_SUB(kpoints_grid_copy)
    
    call kpoints_grid_init(bb%dim, aa, bb%npoints, bb%nshifts)
    
    aa%weight = bb%weight
    aa%point  = bb%point
    aa%red_point = bb%red_point
    aa%shifts = bb%shifts

    POP_SUB(kpoints_grid_copy)
  end subroutine kpoints_grid_copy


  elemental subroutine kpoints_nullify(this)
    type(kpoints_t), intent(out) :: this

    call kpoints_grid_nullify(this%full)
    call kpoints_grid_nullify(this%reduced)
    this%method = 0
    this%use_symmetries = .false.
    this%use_time_reversal = .false.
    this%nik_skip = 0
    this%nik_axis = 0
    nullify(this%symmetry_ops, this%num_symmetry_ops)
    nullify(this%klattice)

  end subroutine kpoints_nullify

  ! ---------------------------------------------------------
  subroutine kpoints_init(this, symm, dim, rlattice, klattice, only_gamma)
    type(kpoints_t),    intent(out) :: this
    type(symmetries_t), intent(in)  :: symm
    integer,            intent(in)  :: dim
    FLOAT,              intent(in)  :: rlattice(:,:), klattice(:,:)
    logical,            intent(in)  :: only_gamma

    integer :: ik, idir, is
    character(len=100) :: str_tmp
    FLOAT :: weight_sum

    PUSH_SUB(kpoints_init)

    ASSERT(dim <= MAX_DIM)

    nullify(this%symmetry_ops)
    nullify(this%num_symmetry_ops)
    nullify(this%klattice)

    !%Variable KPointsUseSymmetries
    !%Type logical
    !%Default no
    !%Section Mesh::KPoints
    !%Description
    !% This variable defines whether symmetries are taken into account
    !% or not for the choice of <i>k</i>-points. If it is set to no, the <i>k</i>-point
    !% sampling will range over the full Brillouin zone.
    !%
    !% When a perturbation is applied to the system, the full
    !% symmetries of the system cannot be used. In this case you must
    !% not use symmetries or use the <tt>SymmetryBreakDir</tt> to tell
    !% Octopus the direction of the perturbation (for the moment this
    !% has to be done by hand by the user, in the future it will be
    !% automatic).
    !%
    !%End
    call parse_variable('KPointsUseSymmetries', .false., this%use_symmetries)

    !%Variable KPointsUseTimeReversal
    !%Type logical
    !%Section Mesh::KPoints
    !%Description
    !% If symmetries are used to reduce the number of <i>k</i>-points,
    !% this variable defines whether time-reversal symmetry is taken
    !% into account or not. If it is set to no, the <i>k</i>-point
    !% sampling will not be reduced according to time-reversal
    !% symmetry.
    !%
    !% The default is yes, unless symmetries are broken in one
    !% direction by the <tt>SymmetryBreakDir</tt> block.
    !% 
    !% Warning: For time propagation runs with an external field,
    !% time-reversal symmetry should not be used.
    !%
    !%End
    call parse_variable('KPointsUseTimeReversal', .not. symmetries_have_break_dir(symm), this%use_time_reversal)

    if(only_gamma) then
      this%method = KPOINTS_GAMMA
      call read_MP(gamma_only = .true.)
    else 
      ! do this always to allow setting KPointsGrid for BerkeleyGW output even if KPoints(Reduced) is also set
      call read_MP(gamma_only = .false.)

      if(read_user_kpoints()) then
        this%method = KPOINTS_USER
      else
        this%method = KPOINTS_MONKH_PACK
 
        write(message(1),'(a)') ' '
        write(message(2),'(1x,i3,a)') this%reduced%npoints, ' k-points generated from parameters :'
        write(message(3),'(1x,a)') '---------------------------------------------------'
        write(message(4),'(4x,a)') 'n ='
        do idir = 1, dim
          write(str_tmp,'(i5)') this%nik_axis(idir)
          message(4) = trim(message(4)) // trim(str_tmp)
        end do
        call messages_info(4)
        
        do is = 1, this%reduced%nshifts
          write(message(1),'(a)') ' '
          write(message(2),'(4x,a,i1,a)') 's', is, '  ='
          do idir = 1, dim
            write(str_tmp,'(f6.2)') this%reduced%shifts(idir,is)
            message(2) = trim(message(2)) // trim(str_tmp)
          end do
          call messages_info(2)
        enddo

      end if

      write(message(1),'(a)') ' '
      write(message(2),'(a)') ' index |    weight    |             coordinates              |'
      call messages_info(2)

      do ik = 1, this%reduced%npoints
        write(str_tmp,'(i6,a,f12.6,a)') ik, " | ", this%reduced%weight(ik), " |"
        message(1) =  str_tmp
        do idir = 1, dim
          write(str_tmp,'(f12.6)') this%reduced%red_point(idir, ik)
          message(1) = trim(message(1)) // trim(str_tmp)
        end do
        write(str_tmp,'(a)') "  |"
        message(1) = trim(message(1)) // trim(str_tmp)
        call messages_info(1)
      end do

      write(message(1),'(a)') ' '
      call messages_info(1)

    end if

    SAFE_ALLOCATE(this%klattice(1:dim, 1:dim))

    this%klattice(1:dim, 1:dim) = klattice(1:dim, 1:dim)

    POP_SUB(kpoints_init)

  contains

    ! ---------------------------------------------------------
    subroutine read_MP(gamma_only)
      logical, intent(in) :: gamma_only

      logical       :: gamma_only_
      integer       :: ii, is, ncols, nshifts
      type(block_t) :: blk
      integer, allocatable :: symm_ops(:, :), num_symm_ops(:)
      FLOAT, allocatable :: shifts(:,:)
      

      PUSH_SUB(kpoints_init.read_MP)

      call messages_obsolete_variable('KPointsMonkhorstPack', 'KPointsGrid')

      !%Variable KPointsGrid
      !%Type block
      !%Default <math>\Gamma</math>-point only
      !%Section Mesh::KPoints
      !%Description
      !% When this block is given (and the <tt>KPoints</tt> block is not present),
      !% <i>k</i>-points are distributed in a uniform grid, according to a modified
      !% version of the Monkhorst-Pack scheme. For the original MP scheme, see
      !% James D. Pack and Hendrik J. Monkhorst,
      !% <i>Phys. Rev. B</i> <b>13</b>, 5188 (1976) and <i>Phys. Rev. B</i> <b>16</b>, 1748 (1977).
      !%
      !% The first row of the block is a set of integers defining
      !% the number of <i>k</i>-points to be used along each direction
      !% in reciprocal space. The numbers refer to the whole Brillouin
      !% zone, and the actual number of <i>k</i>-points is usually
      !% reduced exploiting the symmetries of the system.  By default
      !% the grid will always include the <math>\Gamma</math>-point. Optional
      !% rows can be added to specify multiple shifts in the <i>k</i>-points (between 0.0 and 1.0),
      !% in units of the Brillouin zone divided by the number in the first row.
      !% The number of columns should be equal to <tt>Dimensions</tt>,
      !% but the grid and shift numbers should be 1 and zero in finite directions.
      !%
      !% For example, the following input samples the BZ with 100 points in the 
      !% <i>xy</i>-plane of reciprocal space:
      !%
      !% <tt>%KPointsGrid
      !% <br>&nbsp;&nbsp;10 | 10 | 1
      !% <br>%</tt>
      !%
      !%End

      gamma_only_ = gamma_only
      if(.not. gamma_only_) &
        gamma_only_ = (parse_block('KPointsGrid', blk) /= 0)

      this%nik_axis(1:MAX_DIM) = 1

      if(.not. gamma_only_) then
        nshifts = max(parse_block_n(blk)-1,1) 
      else
        nshifts = 1
      end if
      SAFE_ALLOCATE(shifts(1:MAX_DIM,1:nshifts))
      shifts(1:MAX_DIM,1:nshifts) = M_ZERO

      if(.not. gamma_only_) then
        ncols = parse_block_cols(blk, 0)
        if(ncols /= dim) then
          write(message(1),'(a,i3,a,i3)') 'KPointsGrid first row has ', ncols, ' columns but must have ', dim
          call messages_fatal(1)
        end if
        do ii = 1, dim
          call parse_block_integer(blk, 0, ii - 1, this%nik_axis(ii))
        end do

        if (any(this%nik_axis(1:dim) < 1)) then
          message(1) = 'Input: KPointsGrid is not valid.'
          call messages_fatal(1)
        end if

        if(parse_block_n(blk) > 1) then ! we have a shift, or even more
          ncols = parse_block_cols(blk, 1)
          if(ncols /= dim) then
            write(message(1),'(a,i3,a,i3)') 'KPointsGrid second row has ', ncols, ' columns but must have ', dim
            call messages_fatal(1)
          end if
          do is = 1, nshifts
            do ii = 1, dim
              call parse_block_float(blk, is, ii - 1, shifts(ii,is))
            end do
          end do
        end if

        call parse_block_end(blk)
      end if

      call kpoints_grid_init(dim, this%full, product(this%nik_axis(1:dim))*nshifts, nshifts)

      !We move the k-points into this%shifts
      do is = 1, nshifts
        do ii = 1, dim
          this%full%shifts(ii,is) = shifts(ii,is)
        end do
      end do
      SAFE_DEALLOCATE_A(shifts)

      call kpoints_grid_generate(dim, this%nik_axis(1:dim), this%full%nshifts, &
               this%full%shifts(1:dim,1:this%full%nshifts), this%full%red_point)

      this%full%weight = M_ONE / this%full%npoints

      call kpoints_grid_copy(this%full, this%reduced)

      SAFE_ALLOCATE(this%num_symmetry_ops(1:this%reduced%npoints))

      if(this%use_symmetries) then

        SAFE_ALLOCATE(num_symm_ops(1:this%full%npoints))

        if(this%use_time_reversal) then
          SAFE_ALLOCATE(symm_ops(1:this%full%npoints, 1:2*symmetries_number(symm)))
        else
          SAFE_ALLOCATE(symm_ops(1:this%full%npoints, 1:symmetries_number(symm)))
        end if
        
        call kpoints_grid_reduce(symm, this%use_time_reversal, &
          this%reduced%npoints, dim, this%reduced%red_point, this%reduced%weight, symm_ops, num_symm_ops)
        
        ! sanity checks
        ASSERT(maxval(num_symm_ops) >= 1)
        if(this%use_time_reversal) then
          ASSERT(maxval(num_symm_ops) <= 2*symmetries_number(symm))
        else
          ASSERT(maxval(num_symm_ops) <= symmetries_number(symm))
        end if
        ! the total number of symmetry operations in the list has to be equal to the number of k-points
        ASSERT(sum(num_symm_ops(1:this%reduced%npoints)) == this%full%npoints)

        do ik = 1, this%reduced%npoints
          ASSERT(all(symm_ops(ik, 1:num_symm_ops(ik)) <= symm%nops))
        end do

        SAFE_ALLOCATE(this%symmetry_ops(1:this%reduced%npoints, 1:maxval(num_symm_ops)))
        
        this%num_symmetry_ops(1:this%reduced%npoints) = num_symm_ops(1:this%reduced%npoints)
        this%symmetry_ops(1:this%reduced%npoints, 1:maxval(num_symm_ops)) = &
          symm_ops(1:this%reduced%npoints, 1:maxval(num_symm_ops))
        
        SAFE_DEALLOCATE_A(num_symm_ops)
        SAFE_DEALLOCATE_A(symm_ops)

      else

        SAFE_ALLOCATE(this%symmetry_ops(1:this%reduced%npoints, 1:1))
        this%num_symmetry_ops(1:this%reduced%npoints) = 1
        this%symmetry_ops(1:this%reduced%npoints, 1) = 1
        
      end if
      
      do ik = 1, this%full%npoints
        call kpoints_to_absolute(klattice, this%full%red_point(:, ik), this%full%point(:, ik), dim)
      end do

      do ik = 1, this%reduced%npoints
        call kpoints_to_absolute(klattice, this%reduced%red_point(:, ik), this%reduced%point(:, ik), dim)
      end do

      POP_SUB(kpoints_init.read_MP)
    end subroutine read_MP


    ! ---------------------------------------------------------
    logical function read_user_kpoints()
      type(block_t) :: blk
      logical :: reduced
      integer :: ik, idir

      PUSH_SUB(kpoints_init.read_user_kpoints)

      !%Variable KPoints
      !%Type block
      !%Section Mesh::KPoints
      !%Description
      !% This block defines an explicit set of <i>k</i>-points and their weights for
      !% a periodic-system calculation. The first column is the weight
      !% of each <i>k</i>-point and the following are the components of the <i>k</i>-point
      !% vector. You only need to specify the components for the
      !% periodic directions. Note that the <i>k</i>-points should be given in
      !% Cartesian coordinates (not in reduced coordinates), in the units of inverse length.
      !% The weights will be renormalized so they sum to 1 (and must be rational numbers).
      !%
      !% For example, if you want to include only the Gamma point, you can
      !% use:
      !%
      !% <tt>%KPoints
      !% <br>&nbsp;&nbsp;1.0 | 0 | 0 | 0
      !% <br>%</tt>
      !%
      !%End

      !%Variable KPointsReduced
      !%Type block
      !%Section Mesh::KPoints
      !%Description
      !% Same as the block <tt>KPoints</tt> but this time the input is given in reduced 
      !% coordinates, <i>i.e.</i>
      !% what <tt>Octopus</tt> writes in a line in the ground-state standard output as
      !% 
      !% <tt>#k =   1, k = (    0.154000,    0.154000,    0.154000)</tt>.
      !%End

      reduced = .false.
      if(parse_block('KPoints', blk) /= 0) then
        if(parse_block('KPointsReduced', blk) == 0) then
          reduced = .true.
        else
          read_user_kpoints = .false.
          POP_SUB(kpoints_init.read_user_kpoints)
          return
        end if
      end if
      read_user_kpoints = .true.

      ! end the one initialized by KPointsGrid already
      call kpoints_end(this)

      this%use_symmetries = .false.

      call kpoints_grid_init(dim, this%full, parse_block_n(blk), 1)

      this%full%red_point = M_ZERO
      this%full%point = M_ZERO
      this%full%weight = M_ZERO
      this%full%shifts = M_ZERO

      if(reduced) then
        do ik = 1, this%full%npoints
          call parse_block_float(blk, ik - 1, 0, this%full%weight(ik))
          do idir = 1, dim
            call parse_block_float(blk, ik - 1, idir, this%full%red_point(idir, ik))
          end do
          ! generate also the absolute coordinates
          call kpoints_to_absolute(klattice, this%full%red_point(:, ik), this%full%point(:, ik), dim)
        end do
      else
        do ik = 1, this%full%npoints
          call parse_block_float(blk, ik - 1, 0, this%full%weight(ik))
          do idir = 1, dim
            call parse_block_float(blk, ik - 1, idir, this%full%point(idir, ik), unit_one/units_inp%length)
          end do
          ! generate also the reduced coordinates
          call kpoints_to_reduced(rlattice, this%full%point(:, ik), this%full%red_point(:, ik), dim)
        end do
      end if
      call parse_block_end(blk)

      this%nik_skip = 0
      if(any(this%full%weight(:) < M_EPSILON)) then
        call messages_experimental('K-points with zero weight')
        message(1) = "Found k-points with zero weight. They are excluded from density calculation"
        call messages_warning(1)
        ! count k-points with zero weight and  make sure the points are given in
        ! a block after all regular k-points. This is for convenience, so they can be skipped
        ! easily and not a big restraint for the user who has to provide the k-points
        ! explicitly anyway.
        do ik=1,this%full%npoints
          if(this%full%weight(ik) < M_EPSILON) then
            ! check there are no points with positive weight following a zero weighted one
            if(ik < this%full%npoints) then
              if(this%full%weight(ik+1) > M_EPSILON) then
                message(1) = "K-points with zero weight must follow all regular k-points in a block"
                call messages_fatal(1)
              endif
            end if
            this%nik_skip = this%nik_skip + 1
            ! set to machine zero
            this%full%weight(ik) = M_ZERO
          end if
        end do
      end if
      ! renormalize weights
      weight_sum = sum(this%full%weight(1:this%full%npoints))
      if(weight_sum < M_EPSILON) then
        message(1) = "k-point weights must sum to a positive number."
        call messages_fatal(1)
      end if
      this%full%weight = this%full%weight / weight_sum

      ! for the moment we do not apply symmetries to user kpoints
      call kpoints_grid_copy(this%full, this%reduced)

      write(message(1), '(a,i4,a)') 'Input: ', this%full%npoints, ' k-points were read from the input file'
      call messages_info(1)

      POP_SUB(kpoints_init.read_user_kpoints)
    end function read_user_kpoints

  end subroutine kpoints_init

  ! ---------------------------------------------------------
  subroutine kpoints_end(this)
    type(kpoints_t), intent(inout) :: this

    PUSH_SUB(kpoints_end)

    call kpoints_grid_end(this%full)
    call kpoints_grid_end(this%reduced)

    SAFE_DEALLOCATE_P(this%klattice)
    SAFE_DEALLOCATE_P(this%symmetry_ops)
    SAFE_DEALLOCATE_P(this%num_symmetry_ops)
   

    POP_SUB(kpoints_end)
  end subroutine kpoints_end


  ! ---------------------------------------------------------
  subroutine kpoints_to_absolute(klattice, kin, kout, dim)
    FLOAT,   intent(in)  :: klattice(:,:)
    FLOAT,   intent(in)  :: kin(:)
    FLOAT,   intent(out) :: kout(:)
    integer, intent(in)  :: dim

    integer :: ii
    
    PUSH_SUB(kpoints_to_absolute)

    kout(1:dim) = M_ZERO
    do ii = 1, dim
      kout(1:dim) = kout(1:dim) + kin(ii) * klattice(1:dim, ii)
    end do

    POP_SUB(kpoints_to_absolute)
  end subroutine kpoints_to_absolute


  ! ---------------------------------------------------------
  subroutine kpoints_to_reduced(rlattice, kin, kout, dim)
    FLOAT,   intent(in)  :: rlattice(:,:)
    FLOAT,   intent(in)  :: kin(:)
    FLOAT,   intent(out) :: kout(:)
    integer, intent(in)  :: dim

    integer :: ii

    PUSH_SUB(kpoints_to_reduced)

    kout(1:dim) = M_ZERO
    do ii = 1, dim
      kout(1:dim) = kout(1:dim) + kin(ii) * rlattice(ii, 1:dim)
    end do
    kout(1:dim) = kout(1:dim) / (M_TWO * M_PI)

    POP_SUB(kpoints_to_reduced)
  end subroutine kpoints_to_reduced


  ! ---------------------------------------------------------
  subroutine kpoints_copy(kin, kout)
    type(kpoints_t), intent(in)  :: kin
    type(kpoints_t), intent(out) :: kout

    PUSH_SUB(kpoints_copy)

    kout%method = kin%method

    call kpoints_grid_copy(kin%full, kout%full)
    call kpoints_grid_copy(kin%reduced, kout%reduced)

    kout%use_symmetries = kin%use_symmetries
    kout%use_time_reversal = kin%use_time_reversal

    kout%nik_axis(1:kin%full%dim) = kin%nik_axis(1:kin%full%dim)

    SAFE_ALLOCATE(kout%klattice(1:kin%full%dim, kin%full%dim))
    kout%klattice(1:kin%full%dim, kin%full%dim) = kin%klattice(1:kin%full%dim, kin%full%dim)

    POP_SUB(kpoints_copy)
  end subroutine kpoints_copy


  ! ----------------------------------------------------------
  integer pure function kpoints_number(this) result(number)
    type(kpoints_t), intent(in) :: this

    number = this%reduced%npoints
    
  end function kpoints_number


  ! ----------------------------------------------------------
  pure function kpoints_get_point(this, ik, absolute_coordinates) result(point)
    type(kpoints_t),   intent(in) :: this
    integer,           intent(in) :: ik
    logical, optional, intent(in) :: absolute_coordinates    !< .true. by default
    FLOAT                         :: point(1:this%full%dim)

    if(optional_default(absolute_coordinates, .true.)) then
      point(1:this%full%dim) = this%reduced%point(1:this%full%dim, ik)
    else
      point(1:this%full%dim) = this%reduced%red_point(1:this%full%dim, ik)
    end if

  end function kpoints_get_point


  ! ----------------------------------------------------------
  FLOAT pure function kpoints_get_weight(this, ik) result(weight)
    type(kpoints_t), intent(in) :: this
    integer,         intent(in) :: ik

    weight = this%reduced%weight(ik)

  end function kpoints_get_weight


  ! ----------------------------------------------------------
  !> Generates the k-points grid.
  !! Sets up a uniform array of k-points. Use a modification of the normal Monkhorst-Pack scheme, 
  !! which is equivalent to the normal MP scheme in the case of even number of kpoints (i.e. naxis (i) even)  
  !! used with a shift of (1/2, 1/2, 1/2).
  !! naxis(i) are the number of points in the three directions determined by the lattice vectors.
  !! shift(i) and sz shift the grid of integration points from the origin.
  subroutine kpoints_grid_generate(dim, naxis, nshifts, shift, kpoints, lk123)  
    integer,           intent(in)  :: dim
    integer,           intent(in)  :: naxis(:)
    integer,           intent(in)  :: nshifts
    FLOAT,             intent(in)  :: shift(:,:)
    FLOAT,             intent(out) :: kpoints(:, :)
    integer, optional, intent(out) :: lk123(:,:)      !< lk123(1:nkpt,1:3): maps ik to a triplet of indices on a cube
                                                      !< running from 0 to naxis(1:3).
  
    FLOAT :: dx(1:MAX_DIM), maxcoord
    integer :: ii, jj, divisor, ik, idir, npoints, is
    integer, allocatable :: ix(:), lk123_(:,:),idx(:)
    FLOAT, allocatable :: nrm(:), shell(:), coords(:, :)
    logical, allocatable :: move_to_minus_half(:)

    PUSH_SUB(kpoints_grid_generate)
   
    dx(1:dim) = M_ONE/(M_TWO*naxis(1:dim))

    npoints = product(naxis(1:dim))

    SAFE_ALLOCATE(move_to_minus_half(1:dim))
    SAFE_ALLOCATE(ix(1:dim))
    
    if (present(lk123)) then
      SAFE_ALLOCATE(lk123_(1:npoints*nshifts,1:dim))
      SAFE_ALLOCATE(idx(1:npoints*nshifts))
    end if

    move_to_minus_half(1:dim) = .true.
    
    do is = 1, nshifts
      do ii = 0, npoints - 1
        ik = npoints*is - ii
        jj = ii
        divisor = npoints

        do idir = 1, dim
          divisor = divisor / naxis(idir)
          ix(idir) = jj / divisor + 1
          jj = mod(jj, divisor)

          kpoints(idir, ik) = (M_TWO*ix(idir) - M_ONE*naxis(idir) + M_TWO*shift(idir,is))*dx(idir)

          if(mod(naxis(idir), 2) /= 0) then    
            kpoints(idir, ik) = kpoints(idir, ik) - dx(idir)
          end if
  
          !bring back point to first Brillouin zone, except for points at 1/2
          if ( abs(kpoints(idir, ik) - CNST(0.5)) > CNST(1e-14) )  then
            kpoints(idir, ik) = mod(kpoints(idir, ik) + M_HALF, M_ONE) - M_HALF
          else
            ! alternate the assignation of points at 1/2 and -1/2 such that the total sum of k-points is zero.
            if(move_to_minus_half(idir)) kpoints(idir,ik) = -CNST(0.5)
            move_to_minus_half(idir) = .not. move_to_minus_half(idir)
          end if

        end do
        if (present(lk123)) then
          lk123_(ik, 1:dim) = ix(1:dim)
          idx(ik) = ik
        end if
      end do
    end do

    SAFE_DEALLOCATE_A(ix)    

    ! sort the k-points

    SAFE_ALLOCATE(nrm(1:npoints*nshifts))
    SAFE_ALLOCATE(shell(1:npoints*nshifts))
    SAFE_ALLOCATE(coords(1:dim, 1:npoints*nshifts))
    
    do ik = 1, npoints*nshifts
      shell(ik) = sum((kpoints(1:dim, ik)/dx(1:dim))**2)
      do idir = 1, dim
        coords(idir, ik) = kpoints(idir, ik)
        if(coords(idir, ik) < CNST(0.0)) coords(idir, ik) = coords(idir, ik) + CNST(1.0)
        coords(idir, ik) = coords(idir, ik)/(dx(idir)*CNST(2.0))
      end do
    end do

    nrm(1:npoints*nshifts) = M_ZERO

    maxcoord = CNST(1.0)
    do idir = dim, 1, -1
      do ik = 1, npoints*nshifts
        nrm(ik) = nrm(ik) + coords(idir, ik)*maxcoord
      end do
      maxcoord = maxcoord*max(CNST(1.0), maxval(coords(idir, 1:npoints*nshifts)))
    end do

    do ik = 1, npoints*nshifts
      nrm(ik) = nrm(ik) + shell(ik)*maxcoord
    end do

    if (present(lk123)) then
      call sort(nrm, idx)      
      do ik = 1, npoints*nshifts
        lk123(ik,1:dim) = lk123_(idx(ik),1:dim)
      end do
      SAFE_DEALLOCATE_A(lk123_)
      SAFE_DEALLOCATE_A(idx)
    end if

    call sort(nrm, kpoints)      


    SAFE_DEALLOCATE_A(nrm)
    SAFE_DEALLOCATE_A(shell)
    SAFE_DEALLOCATE_A(coords)
    SAFE_DEALLOCATE_A(move_to_minus_half)

    POP_SUB(kpoints_grid_generate)
  end subroutine kpoints_grid_generate

  
  ! --------------------------------------------------------------------------------------------
  subroutine kpoints_grid_reduce(symm, time_reversal, nkpoints, dim, kpoints, weights, symm_ops, num_symm_ops)
    type(symmetries_t), intent(in)    :: symm
    logical,            intent(in)    :: time_reversal
    integer,            intent(inout) :: nkpoints
    integer,            intent(in)    :: dim
    FLOAT,              intent(inout) :: kpoints(:, :)
    FLOAT,              intent(out)   :: weights(:)
    integer,            intent(out)   :: symm_ops(:, :)
    integer,            intent(out)   :: num_symm_ops(:)

    integer :: nreduced
    FLOAT, allocatable :: reduced(:, :)
    
    FLOAT :: dw
    integer ik, iop, ik2
    FLOAT :: tran(MAX_DIM), tran_inv(MAX_DIM)
    integer, allocatable :: kmap(:)

    PUSH_SUB(kpoints_grid_reduce)

    ! reduce to irreducible zone

    ! kmap is used to mark reducible k-points and also to
    ! map reducible to irreducible k-points

    SAFE_ALLOCATE(kmap(1:nkpoints))
    SAFE_ALLOCATE(reduced(1:dim, 1:nkpoints))

    forall(ik = 1:nkpoints) kmap(ik) = ik
    
    dw = M_ONE / nkpoints

    nreduced = 0

    num_symm_ops = 1
    symm_ops(:, 1) = symmetries_identity_index(symm)

    do ik = 1, nkpoints
      if (kmap(ik) /= ik) cycle
      
      ! new irreducible point
      ! mark reduced with negative kmap
      
      nreduced = nreduced + 1
      reduced(1:dim, nreduced) = kpoints(1:dim, ik)
      
      kmap(ik) = -nreduced
      
      weights(nreduced) = dw

      if (ik == nkpoints) cycle
      
      ! operate with the symmetry operations
      
      do iop = 1, symmetries_number(symm)
        if(iop == symmetries_identity_index(symm)) cycle ! no need to check for the identity

        call symmetries_apply_kpoint(symm, iop, reduced(1:dim, nreduced), tran)
        tran_inv(1:dim) = -tran(1:dim)
           
        ! remove (mark) k-points related to irreducible reduced by symmetry
        do ik2 = ik + 1, nkpoints
          if (kmap(ik2) /= ik2) cycle
          
          ! both the transformed rk ...
          if (all( abs(tran(1:dim) - kpoints(1:dim, ik2)) <= CNST(1.0e-5))) then
            weights(nreduced) = weights(nreduced) + dw
            kmap(ik2) = nreduced
            INCR(num_symm_ops(nreduced), 1)
            symm_ops(nreduced, num_symm_ops(nreduced)) = iop
            cycle
          end if

          ! and its inverse
          if (time_reversal .and. all(abs(tran_inv(1:dim) - kpoints(1:dim, ik2)) <= CNST(1.0e-5)) ) then
            weights(nreduced) = weights(nreduced) + dw
            kmap(ik2) = nreduced
            INCR(num_symm_ops(nreduced), 1)
            symm_ops(nreduced, num_symm_ops(nreduced)) = iop
          end if
          
        end do
      end do
    end do
    
    nkpoints = nreduced
    do ik = 1, nreduced
      kpoints(1:dim, ik) = reduced(1:dim, ik)
    end do

    POP_SUB(kpoints_grid_reduce)
  end subroutine kpoints_grid_reduce


  ! ---------------------------------------------------------
  subroutine kpoints_write_info(this, iunit, absolute_coordinates)
    type(kpoints_t),    intent(in) :: this
    integer,            intent(in) :: iunit
    logical, optional,  intent(in) :: absolute_coordinates
    
    integer :: ik, idir
    character(len=100) :: str_tmp
    character :: index
    
    PUSH_SUB(kpoints_write_info)
    
    call messages_print_stress(iunit, 'Brillouin zone sampling')

    if(this%method == KPOINTS_MONKH_PACK) then

      call messages_write('Dimensions of the k-point grid      =')
      do idir = 1, this%full%dim
        call messages_write(this%nik_axis(idir), fmt = '(i3,1x)')
      end do
      call messages_new_line()
      
      call messages_write('Total number of k-points            =')
      call messages_write(this%full%npoints)
      call messages_new_line()

      call messages_write('Number of symmetry-reduced k-points =')
      call messages_write(this%reduced%npoints)

      call messages_info(iunit = iunit)
      
    else

      call messages_write('Total number of k-points            =')
      call messages_write(this%full%npoints)
      call messages_new_line()
      call messages_info(iunit = iunit)

    end if

    call messages_new_line()
    call messages_write('List of k-points:')
    call messages_info(iunit = iunit)

    write(message(1), '(6x,a)') 'ik'
    do idir = 1, this%full%dim
      index = index2axis(idir)
      write(str_tmp, '(9x,2a)') 'k_', index
      message(1) = trim(message(1)) // trim(str_tmp)
    end do
    write(str_tmp, '(6x,a)') 'Weight'
    message(1) = trim(message(1)) // trim(str_tmp)
    message(2) = '---------------------------------------------------------'
    call messages_info(2, iunit)
    
    do ik = 1, kpoints_number(this)
      write(message(1),'(i8,1x)') ik
      do idir = 1, this%full%dim
        if(optional_default(absolute_coordinates, .false.)) then
          write(str_tmp,'(f12.4)') this%reduced%point(idir, ik)
        else  
          write(str_tmp,'(f12.4)') this%reduced%red_point(idir, ik)
        end if
        message(1) = trim(message(1)) // trim(str_tmp)
      end do
      write(str_tmp,'(f12.4)') kpoints_get_weight(this, ik)
      message(1) = trim(message(1)) // trim(str_tmp)
      call messages_info(1, iunit)
    end do

    call messages_info(iunit = iunit)

    call messages_print_stress(iunit)

    POP_SUB(kpoints_write_info)
  end subroutine kpoints_write_info
  

  ! ---------------------------------------------------------
  logical pure function kpoints_point_is_gamma(this, ik) result(is_gamma)
    type(kpoints_t),    intent(in) :: this
    integer,            intent(in) :: ik
    
    is_gamma = (maxval(abs(kpoints_get_point(this, ik))) < M_EPSILON)
    
  end function kpoints_point_is_gamma

  !--------------------------------------------------------

  integer pure function kpoints_get_num_symmetry_ops(this, ik) result(num)
    type(kpoints_t),    intent(in) :: this
    integer,            intent(in) :: ik

   if(this%use_symmetries) then
     num = this%num_symmetry_ops(ik)
   else
     num = 1
   end if

  end function kpoints_get_num_symmetry_ops
  !--------------------------------------------------------

  integer pure function kpoints_get_symmetry_ops(this, ik, index) result(iop)
    type(kpoints_t),    intent(in) :: this
    integer,            intent(in) :: ik
    integer,            intent(in) :: index

    if(this%use_symmetries) then
      iop = this%symmetry_ops(ik, index)
    else
      iop = 1
    end if

  end function kpoints_get_symmetry_ops

  !--------------------------------------------------------
  integer function kpoints_kweight_denominator(this)
    type(kpoints_t),    intent(in) :: this

    integer :: denom, nik, nik_skip

    PUSH_SUB(kpoints_kweight_denominator)

    nik = this%full%npoints
    nik_skip = this%nik_skip
    
    if(this%method == KPOINTS_MONKH_PACK) then
      kpoints_kweight_denominator = this%full%npoints
    else
      kpoints_kweight_denominator = 0
      ! NB largest reasonable value is: # k-points x 48. from space-group symmetries
      do denom = 1, 100000
        if(all(abs(int(this%full%weight(1:nik-nik_skip)*denom + CNST(10)*M_EPSILON) - &
          this%full%weight(1:nik-nik_skip)*denom) < CNST(100)*M_EPSILON)) then
          kpoints_kweight_denominator = denom
          exit
        end if
      end do
    end if

    POP_SUB(kpoints_kweight_denominator)
  end function kpoints_kweight_denominator

  !--------------------------------------------------------
  logical  pure function kpoints_have_zero_weight_path(this) result(have_zerow)
    type(kpoints_t),    intent(in) :: this
    
    if (this%nik_skip > 0) then
      have_zerow = .true.
    else 
      have_zerow = .false.
    end if

  end function kpoints_have_zero_weight_path

end module kpoints_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
