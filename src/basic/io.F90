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
!! $Id: io.F90 15203 2016-03-19 13:15:05Z xavier $

#include "global.h"

module io_oct_m
  use debug_oct_m
  use global_oct_m
  use loct_oct_m
  use messages_oct_m
  use mpi_oct_m
  use parser_oct_m

  implicit none

  private
  public ::              &
    io_workpath,         &
    io_open,             &
    io_mkdir,            &
    io_init,             &
    io_end,              &
    io_status,           &
    io_dump_file,        &
    io_free,             &
    io_close,            &
    io_assign,           &
    io_get_extension,    &
    io_switch_status,    &
    io_debug_on_the_fly, &
    iopar_read,          &
    iopar_backspace,     &
    iopar_find_line,     &
    io_skip_header,      &
    io_file_exists

  integer, parameter :: min_lun=10, max_lun=99
  logical            :: lun_is_free(min_lun:max_lun)
  character(len=MAX_PATH_LEN) :: work_dir    !< name of the output directory

contains

  ! ---------------------------------------------------------
  !> If the argument defaults is present and set to true, then the routine
  !! will not try to read anything from the inp file, but set everything
  !! to the default values.
  subroutine io_init(defaults)
    logical, optional, intent(in) :: defaults

    character(len=MAX_PATH_LEN) :: filename
    character(len=256) :: node_hook
    logical :: file_exists, mpi_debug_hook
    integer :: sec, usec

    ! cannot use push/pop before initializing io

    if(present(defaults)) then
      if(defaults) then
        lun_is_free(min_lun:max_lun)=.true.
        stdin  = 5
        stdout = 6
        stderr = 0
        work_dir = '.'
        flush_messages = .false.
        return
      end if
    end if

    lun_is_free(min_lun:max_lun)=.true.
    stdin = 5

    !%Variable stdout
    !%Type string
    !%Default "-"
    !%Section Execution::IO
    !%Description
    !% The standard output by default goes to, well, to standard output. This can
    !% be changed by setting this variable: if you give it a name (other than "-")
    !% the output stream is printed in that file instead.
    !%End
    call parse_variable('stdout', '-', filename)
    stdout = 6
    if(trim(filename) /= '-') then
      close(stdout)
      open(stdout, file=filename, status='unknown')
    end if

    !%Variable stderr
    !%Type string
    !%Default "-"
    !%Section Execution::IO
    !%Description
    !% The standard error by default goes to, well, to standard error. This can
    !% be changed by setting this variable: if you give it a name (other than "-")
    !% the output stream is printed in that file instead.
    !%End
    call parse_variable('stderr', '-', filename)
    stderr = 0
    if(trim(filename) /= '-') then
      close(stderr)
      open(stderr, file=filename, status='unknown')
    end if

    !%Variable WorkDir
    !%Type string
    !%Default "."
    !%Section Execution::IO
    !%Description
    !% By default, all files are written and read from the working directory,
    !% <i>i.e.</i> the directory from which the executable was launched. This behavior can
    !% be changed by setting this variable: if you give it a name (other than ".")
    !% the files are written and read in that directory.
    !%End
    call parse_variable('WorkDir', '.', work_dir)
    ! ... and if necessary create workdir (will not harm if work_dir is already there)
    if (work_dir /= '.') call loct_mkdir(trim(work_dir))

    !%Variable FlushMessages
    !%Type logical
    !%Default no
    !%Section Execution::IO
    !%Description
    !% In addition to writing to stdout and stderr, the code messages may also be
    !% flushed to <tt>messages.stdout</tt> and <tt>messages.stderr</tt>, if this variable is
    !% set to yes.
    !%End
    call parse_variable('FlushMessages', .false., flush_messages)

    ! delete files so that we start writing to empty ones
    if(flush_messages) then
      call loct_rm('messages.stdout')
      call loct_rm('messages.stderr')
    end if

    if(debug%info) then
      call io_mkdir('debug')
    end if

    if(debug%trace_file) then
      !wipe out debug trace files from previous runs to start fresh rather than appending
      call delete_debug_trace()
    end if

    ! create static directory
    call io_mkdir(STATIC_DIR)

    if(debug%info) then
      !%Variable MPIDebugHook
      !%Type logical
      !%Default no
      !%Section Execution::Debug
      !%Description
      !% When debugging the code in parallel it is usually difficult to find the origin
      !% of race conditions that appear in MPI communications. This variable introduces 
      !% a facility to control separate MPI processes. If set to yes, all nodes will 
      !% start up, but will get trapped in an endless loop. In every cycle of the loop 
      !% each node is sleeping for one second and is then checking if a file with the 
      !% name <tt>node_hook.xxx</tt> (where <tt>xxx</tt> denotes the node number) exists. A given node can 
      !% only be released from the loop if the corresponding file is created. This allows 
      !% to selectively run, <i>e.g.</i>, a compute node first followed by the master node. Or, by
      !% reversing the file creation of the node hooks, to run the master first followed
      !% by a compute node.
      !%End
      call parse_variable('MPIDebugHook', .false., mpi_debug_hook)
      if (mpi_debug_hook) then
        call loct_gettimeofday(sec, usec)
        call epoch_time_diff(sec,usec)
        write(message(1),'(a,i6,a,i6.6,20x,a)') '* I ',sec,'.',usec,' | MPI debug hook'
        call messages_debug(1)

        write(stdout,'(a,i3,a)') 'node:', mpi_world%rank, ' In debug hook'
        write(node_hook,'(i3.3)') mpi_world%rank
        file_exists = .false.

        do while (.not.file_exists)
          inquire(file='node_hook.'//node_hook, exist=file_exists)
          call loct_nanosleep(1,0)
          write(stdout,'(a,i3,a)') 'node:', mpi_world%rank, &
            ' - still sleeping. To release me touch: node_hook.'//trim(node_hook)
        end do

        write(stdout,'(a,i3,a)') 'node:', mpi_world%rank, ' Leaving debug hook'
        ! remove possible debug hooks
        call loct_rm( 'node_hook.'//trim(node_hook) )

        call loct_gettimeofday(sec, usec)
        call epoch_time_diff(sec,usec)
        write(message(1),'(a,i6,a,i6.6,20x,a)') '* O ', sec, '.', usec,' | MPI debug hook'
        call messages_debug(1)
      end if
    end if

  end subroutine io_init

  ! ---------------------------------------------------------
  subroutine io_end()
    
    ! no PUSH/POP, because the POP would write to stderr after it was closed.

    if(stderr /= 0) call io_close(stderr)
    if(stdin  /= 5) call io_close(stdin)
    if(stdout /= 6) call io_close(stdout)

  end subroutine io_end


  ! ---------------------------------------------------------
  subroutine io_assign(got_lun)
    integer, intent(out) :: got_lun

    integer :: iostat, lun
    logical :: used

    PUSH_SUB(io_assign)

    got_lun = -1

    ! Looks for a free unit and assigns it to lun
    do lun = min_lun, max_lun
      if (lun_is_free(lun)) then
        inquire(unit=lun, opened=used, iostat=iostat)

        if (iostat /= 0) used = .true.
        lun_is_free(lun) = .false.
        if (.not. used) then
          got_lun = lun
          exit
        end if
      end if
    end do

    POP_SUB(io_assign)

  end subroutine io_assign


  ! ---------------------------------------------------------
  subroutine io_free(lun)
    integer, intent(in) :: lun

    PUSH_SUB(io_free)

    if (lun >= min_lun .and. lun  <=  max_lun) &
      lun_is_free(lun) = .true.

    POP_SUB(io_free)

  end subroutine io_free


  ! ---------------------------------------------------------
  character(len=MAX_PATH_LEN) function io_workpath(path) result(wpath)
    character(len=*),  intent(in) :: path

    PUSH_SUB(io_workpath)

    if(path(1:1)  ==  '/') then
      ! we do not change absolute path names
      wpath = trim(path)
    else
      write(wpath, '(3a)') trim(work_dir), "/", trim(path)
    end if

    POP_SUB(io_workpath)

  end function io_workpath


  ! ---------------------------------------------------------
  subroutine io_mkdir(fname, parents)
    character(len=*),  intent(in) :: fname
    logical, optional, intent(in) :: parents

    logical :: parents_
    integer :: last_slash, pos, length

    PUSH_SUB(io_mkdir)

    parents_ = .false.
    if (present(parents)) parents_ = parents

    if (.not. parents_) then
      call loct_mkdir(trim(io_workpath(fname)))
    else
      last_slash = max(index(fname, "/", .true.), len_trim(fname))
      pos = 1
      length = index(fname, '/') - 1
      do while (pos < last_slash)
        call loct_mkdir(trim(io_workpath(fname(1:pos+length-1))))
        pos = pos + length + 1
        length = index(fname(pos:), "/") - 1
        if (length < 1) length = len_trim(fname(pos:))
      end do

    end if

    POP_SUB(io_mkdir)
  end subroutine io_mkdir


  ! ---------------------------------------------------------
  integer function io_open(file, action, status, form, position, die, recl, grp) result(iunit)
    character(len=*), intent(in) :: file, action
    character(len=*), intent(in), optional :: status, form, position
    logical,          intent(in), optional :: die
    integer,          intent(in), optional :: recl
    type(mpi_grp_t),  intent(in), optional :: grp

    character(len=20)  :: status_, form_, position_
    character(len=MAX_PATH_LEN) :: file_
    logical            :: die_
    integer            :: iostat
    type(mpi_grp_t)    :: grp_

    PUSH_SUB(io_open)

    if(present(grp)) then
      grp_%comm = grp%comm
      grp_%rank = grp%rank
      grp_%size = grp%size
    else
      call mpi_grp_init(grp_, -1)
    end if


    if(mpi_grp_is_root(grp_)) then

      status_ = 'unknown'
      if(present(status  )) status_   = status
      form_   = 'formatted'
      if(present(form    )) form_     = form
      position_ = 'asis'
      if(present(position)) position_ = position
      die_    = .true.
      if(present(die     )) die_      = die

      call io_assign(iunit)
      if(iunit<0) then
        if(die_) then
          write(message(1), '(a)') '*** IO Error: Too many files open.'
          call messages_fatal(1)
        end if
        POP_SUB(io_open)
        return
      end if

      file_ = io_workpath(file)

      if(present(recl)) then
        open(unit=iunit, file=trim(file_), status=trim(status_), form=trim(form_), &
          recl=recl, action=trim(action), position=trim(position_), iostat=iostat)
      else
        open(unit=iunit, file=trim(file_), status=trim(status_), form=trim(form_), &
          action=trim(action), position=trim(position_), iostat=iostat)
      end if
      
      if(iostat /= 0) then
        call io_free(iunit)
        iunit = -1
        if(die_) then
          write(message(1), '(5a,i6)') '*** IO Error: Could not open file "', trim(file_), &
            '" for action="', trim(action), '". Error code = ', iostat
          call messages_fatal(1)
        end if
      end if

    end if

#if defined(HAVE_MPI)
    if(grp_%size > 1) then
      call MPI_Bcast(iunit, 1, MPI_INTEGER, 0, grp_%comm, mpi_err)
    end if
#endif

    POP_SUB(io_open)
    
  end function io_open


  ! ---------------------------------------------------------
  subroutine io_close(iunit, grp)
    integer, intent(inout) :: iunit
    type(mpi_grp_t),  intent(in), optional :: grp

    type(mpi_grp_t)    :: grp_

    PUSH_SUB(io_close)

    if(present(grp)) then
      grp_%comm = grp%comm
      grp_%rank = grp%rank
      grp_%size = grp%size
    else
      call mpi_grp_init(grp_, -1)
    end if

    if(mpi_grp_is_root(grp_)) then
      close(iunit)
      call io_free(iunit)
      iunit = -1
    end if

#if defined(HAVE_MPI)
    if(grp_%size > 1) then
      call MPI_Bcast(iunit, 1, MPI_INTEGER, 0, grp_%comm, mpi_err)
    end if
#endif

    POP_SUB(io_close)

  end subroutine io_close


  ! ---------------------------------------------------------
  !> Prints a list of the connected logical units and the names of
  !! the associated files
  ! ---------------------------------------------------------
  subroutine io_status(iunit)
    integer, intent(in) :: iunit

    integer :: ii, iostat
    logical :: opened, named
    character(len=MAX_PATH_LEN) :: filename
    character(len=11) :: form

    PUSH_SUB(io_status)

    write(iunit, '(a)') '******** io_status ********'
    do ii = 0, max_lun
      inquire(ii, opened=opened, named=named, name=filename, form=form, iostat=iostat)
      if (iostat  ==  0) then
        if (opened) then
          if(.not. named) filename = 'No name available'
          write(iunit, '(i4,5x,a,5x,a)') ii, form, filename
        end if
      else
        write(iunit, '(i4,5x,a)') ii, 'Iostat error'
      end if
    end do
    write(iunit,'(a)') '********           ********'

    POP_SUB(io_status)

  end subroutine io_status


  ! ---------------------------------------------------------
  subroutine io_dump_file(ounit, filename)
    integer,          intent(in) :: ounit
    character(len=*), intent(in) :: filename

    integer :: iunit, err
    character(len=80) :: line

    if(.not. mpi_grp_is_root(mpi_world)) return

    PUSH_SUB(io_dump_file)

    call io_assign(iunit)
    open(unit=iunit, file=filename, iostat=err, action='read', status='old')

    if(flush_messages) then
      open(unit=iunit_out, file='messages.stdout', &
        action='write', position='append')
    end if

    do while(err == 0)
      read(iunit, fmt='(a80)', iostat=err) line
      if(err==0) then
        write(ounit, '(a)') trim(line)
        if(flush_messages) then
          write(iunit_out, '(a)') trim(line)
        end if
      end if
    end do
    
    if(flush_messages) then
      close(iunit_out)
    end if

    call io_close(iunit)
    POP_SUB(io_dump_file)

  end subroutine io_dump_file


  ! ---------------------------------------------------------
  !> Given a path, it returns the extension (if it exists) of the file
  !! (that is, the part of the name that comes after its last point).
  !! If the filename does not have an extension, it returns the empty string.
  character(len=8) function io_get_extension(path) result(ext)
    character(len = * ), intent(in)  :: path
    integer :: i, j

    PUSH_SUB(io_get_extension)

    i = index(path, ".", back = .true.)
    j = index(path(i+1:), "/")
    if(i == 0 .or. j /= 0) then
      ext = ""
    else
      ext = path(i+1:)
    end if

    POP_SUB(io_get_extension)

  end function io_get_extension


  ! ---------------------------------------------------------
  !> create status file for asynchronous communication with GUI
  subroutine io_switch_status(status)
    character(len=*), intent(in) :: status
    
    integer :: iunit

    PUSH_SUB(io_switch_status)

    ! only root node is taking care of file I/O
    if(mpi_grp_is_root(mpi_world)) then 
      ! remove possible leftovers first before we switch to new status
      call loct_rm_status_files()

      ! create empty status file 
      iunit = io_open('exec/oct-status-'//trim(status), &
        action='write', status='new')
      call io_close(iunit)
    end if

    POP_SUB(io_switch_status)

  end subroutine io_switch_status


  ! ---------------------------------------------------------
  !> check if debug mode or message flushing should be enabled or 
  !! disabled on the fly
  subroutine io_debug_on_the_fly()

    PUSH_SUB(io_debug_on_the_fly)

    ! only root node performs the check
    if(mpi_grp_is_root(mpi_world)) then
      if(io_file_exists('enable_debug_mode', msg='Enabling DebugMode')) then
        call debug_enable(debug)
        ! this call does not hurt if the directory is already there
        ! but is otherwise required
        call io_mkdir('debug')
        ! we have been notified by the user, so we can cleanup the file
        call loct_rm('enable_debug_mode')
        ! artificially increase sub stack to avoid underflow
        no_sub_stack = no_sub_stack + 8
      end if

      if(io_file_exists('enable_flush_messages', msg='Enabling flushing of messages')) then
        flush_messages   = .true.
        ! we have been notified by the user, so we can cleanup the file
        call loct_rm('enable_flush_messages')
      end if

      if(io_file_exists('disable_debug_mode', msg='Disabling DebugMode')) then
        call debug_disable(debug)
        ! we have been notified by the user, so we can cleanup the file
        call loct_rm('disable_debug_mode')
      end if

      if(io_file_exists('disable_flush_messages', msg='Disabling flushing of messages')) then
        flush_messages   = .false.
        ! we have been notified by the user, so we can cleanup the file
        call loct_rm('disable_flush_messages')
      end if
    end if

    POP_SUB(io_debug_on_the_fly)

  end subroutine io_debug_on_the_fly

  
  !> Returns true if a file with name 'filename' exists
  !! and issues a reminder.
  ! ---------------------------------------------------------
  logical function io_file_exists(filename, msg) result(file_exists)
    character(len=*),           intent(in)  :: filename
    character(len=*), optional, intent(in)  :: msg

    PUSH_SUB(io_file_exists)

    file_exists = .false.
    inquire(file=trim(filename), exist=file_exists)
    if(file_exists .and. present(msg)) then
      message(1) = trim(msg)
      call messages_warning(1)
    end if

    POP_SUB(io_file_exists)
  end function io_file_exists


  ! ---------------------------------------------------------
  subroutine iopar_read(grp, iunit, lines, n_lines, ierr)
    type(mpi_grp_t),  intent(in)  :: grp
    integer,          intent(in)  :: iunit
    character(len=*), intent(out) :: lines(:)
    integer,          intent(in)  :: n_lines
    integer,          intent(out) :: ierr

    integer :: il

    PUSH_SUB(iopar_read)

    ASSERT(n_lines <= size(lines))

    if(mpi_grp_is_root(grp)) then
      do il = 1, n_lines
        read(iunit, '(a)', iostat=ierr) lines(il)
        if (ierr /= 0) exit
      end do
    end if

#if defined(HAVE_MPI)
    if(grp%size > 1) then
      call MPI_Bcast(ierr, 1, MPI_INTEGER, 0, grp%comm, mpi_err)
      call MPI_Bcast(lines(1), len(lines(1))*n_lines, MPI_CHARACTER, 0, grp%comm, mpi_err)
      ! MPI_Bcast is not a synchronization point, this can cause
      ! problems when iopar_read is called several times (in a loop,
      ! for example). So we add a barrier to make sure the calls are properly synchronized.
      call MPI_Barrier(grp%comm, mpi_err)
    end if
#endif

    POP_SUB(iopar_read)
  end subroutine iopar_read

  ! ---------------------------------------------------------
  subroutine iopar_backspace(grp, iunit)
    type(mpi_grp_t), intent(in)  :: grp
    integer,         intent(in)  :: iunit

    PUSH_SUB(iopar_backspace)

    if(mpi_grp_is_root(grp)) then
      backspace(iunit)
    end if

    POP_SUB(iopar_backspace)

  end subroutine iopar_backspace


  ! ---------------------------------------------------------
  subroutine iopar_find_line(grp, iunit, line, ierr)
    type(mpi_grp_t),  intent(in)  :: grp
    integer,          intent(in)  :: iunit
    character(len=*), intent(in)  :: line
    integer,          intent(out) :: ierr

    character(len=80) :: read_line

    PUSH_SUB(io_find_line)

    if(mpi_grp_is_root(grp)) then
      rewind(iunit)
      do
        read(iunit, '(a)', iostat=ierr) read_line
        if (ierr /= 0 .or. trim(line)  ==  trim(read_line)) exit
      end do
    end if

#if defined(HAVE_MPI)
    if(grp%size > 1) then
      call MPI_Bcast(ierr, 1, MPI_INTEGER, 0, grp%comm, mpi_err)
      ! MPI_Bcast is not a synchronization point, so we add a barrier 
      ! to make sure the calls are properly synchronized.
      call MPI_Barrier(grp%comm, mpi_err)
    end if
#endif

    POP_SUB(io_find_line)
  end subroutine iopar_find_line


  ! ---------------------------------------------------------
  subroutine io_skip_header(iunit)
    integer, intent(in) :: iunit

    character(len=1) :: a

    PUSH_SUB(io_skip_header)

    rewind(iunit)
    read(iunit,'(a)') a
    do while(a=='#')
      read(iunit,'(a)') a
    end do
    backspace(iunit)

    POP_SUB(io_skip_header)

  end subroutine io_skip_header

end module io_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
