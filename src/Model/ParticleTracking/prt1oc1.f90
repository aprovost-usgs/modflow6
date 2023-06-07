module PrtOcModule

  use BaseDisModule, only: DisBaseType
  use KindModule, only: DP, I4B
  use ConstantsModule, only: LENMODELNAME, MNORMAL
  use OutputControlModule, only: OutputControlType
  use OutputControlDataModule, only: OutputControlDataType, ocd_cr
  use SimVariablesModule, only: errmsg, warnmsg

  implicit none
  private
  public PrtOcType, oc_cr

  !> @ brief Output control for PRT
  !!
  !!  Concrete implementation of OutputControlType for the
  !!  PRT Model
  !<
  type, extends(OutputControlType) :: PrtOcType

  contains
    procedure :: oc_ar
    procedure :: read_options => prt_oc_read_options
  end type PrtOcType

contains

  !> @ brief Create GwtOcType
  !!
  !!  Create by allocating a new PrtOcType object and initializing
  !!  member variables.
  !!
  !<
  subroutine oc_cr(ocobj, name_model, inunit, iout)
    ! -- dummy
    type(PrtOcType), pointer :: ocobj !< PrtOcType object
    character(len=*), intent(in) :: name_model !< name of the model
    integer(I4B), intent(in) :: inunit !< unit number for input
    integer(I4B), intent(in) :: iout !< unit number for output
    !
    ! -- Create the object
    allocate (ocobj)
    !
    ! -- Allocate scalars
    call ocobj%allocate_scalars(name_model)
    !
    ! -- Save unit numbers
    ocobj%inunit = inunit
    ocobj%iout = iout
    !
    ! -- Initialize block parser
    call ocobj%parser%Initialize(inunit, iout)
    !
    ! -- Return
    return
  end subroutine oc_cr

  !> @ brief Allocate and read PrtOcType
  !!
  !!  Setup concentration, budget, and particle tracks as output control variables.
  !!  todo: how to pass track data to OC? each column as an array via init_dbl etc?
  !<
  subroutine oc_ar(this, conc, dis, dnodata)
    ! -- dummy
    class(PrtOcType) :: this !< PrtOcType object
    real(DP), dimension(:), pointer, contiguous, intent(in) :: conc !< model concentration
    class(DisBaseType), pointer, intent(in) :: dis !< model discretization package
    real(DP), intent(in) :: dnodata !< no data value
    ! -- local
    integer(I4B) :: i, nocdobj, inodata
    type(OutputControlDataType), pointer :: ocdobjptr
    real(DP), dimension(:), pointer, contiguous :: nullvec => null()
    !
    ! -- Initialize variables
    inodata = 0
    nocdobj = 2
    allocate (this%ocdobj(nocdobj))
    do i = 1, nocdobj
      call ocd_cr(ocdobjptr)
      select case (i)
      case (1)
        call ocdobjptr%init_dbl('BUDGET', nullvec, dis, 'PRINT LAST ', &
                                'COLUMNS 10 WIDTH 11 DIGITS 4 GENERAL ', &
                                this%iout, dnodata)
      case (2)
        call ocdobjptr%init_dbl('CONCENTRATION', conc, dis, 'PRINT LAST ', &
                                'COLUMNS 10 WIDTH 11 DIGITS 4 GENERAL ', &
                                this%iout, dnodata)
      end select
      this%ocdobj(i) = ocdobjptr
      deallocate (ocdobjptr)
    end do
    !
    ! -- Read options or set defaults if this package not on
    if (this%inunit > 0) then
      call this%read_options()
    end if
    !
    ! -- Return
    return
  end subroutine oc_ar

  subroutine prt_oc_read_options(this)
    ! -- modules
    use ConstantsModule, only: LINELENGTH
    use TrackDataModule, only: TRACKHEADERS, TRACKTYPES
    use SimModule, only: store_error, store_error_unit
    use InputOutputModule, only: openfile, getunit
    ! -- dummy
    class(PrtOcType) :: this
    ! -- local
    character(len=LINELENGTH) :: keyword
    character(len=LINELENGTH) :: keyword2
    character(len=LINELENGTH) :: fname
    character(len=:), allocatable :: line
    integer(I4B) :: ierr
    integer(I4B) :: ipos
    logical :: isfound, found, endOfBlock
    type(OutputControlDataType), pointer :: ocdobjptr
    ! -- formats
    character(len=*), parameter :: fmttrkbin = &
      "(4x, 'PARTICLE TRACKS WILL BE SAVED TO BINARY FILE: ', a, /4x, &
    &'OPENED ON UNIT: ', I0)"
    character(len=*), parameter :: fmttrkcsv = &
      "(4x, 'PARTICLE TRACKS WILL BE SAVED TO CSV FILE: ', a, /4x, &
    &'OPENED ON UNIT: ', I0)"
    !
    ! -- get options block
    call this%parser%GetBlock('OPTIONS', isfound, ierr, &
                              supportOpenClose=.true., blockRequired=.false.)
    !
    ! -- parse options block if detected
    if (isfound) then
      write (this%iout, '(/,1x,a,/)') 'PROCESSING OC OPTIONS'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call this%parser%GetStringCaps(keyword)
        found = .false.
        select case (keyword)
        case ('BUDGETCSV')
          call this%parser%GetStringCaps(keyword2)
          if (keyword2 /= 'FILEOUT') then
            errmsg = "BUDGETCSV must be followed by FILEOUT and then budget &
              &csv file name.  Found '"//trim(keyword2)//"'."
            call store_error(errmsg)
            call this%parser%StoreErrorUnit()
          end if
          call this%parser%GetString(fname)
          this%ibudcsv = GetUnit()
          call openfile(this%ibudcsv, this%iout, fname, 'CSV', &
                        filstat_opt='REPLACE')
          found = .true.
          ! case ('TRACK')
          !   call this%parser%GetStringCaps(keyword)
          !   if (keyword == 'FILEOUT') then
          !     call this%parser%GetString(fname)
          !     this%itrkout = getunit()
          !     call openfile(this%itrkout, this%iout, fname, 'DATA(BINARY)', &
          !                   form, access, filstat_opt='REPLACE', &
          !                   mode_opt=MNORMAL)
          !     write (this%iout, fmttrkbin) trim(adjustl(fname)), this%itrkout
          !     ! open and write ascii header file
          !     this%itrkhdr = getunit()
          !     fname = trim(fname)//'.hdr'
          !     call openfile(this%itrkhdr, this%iout, fname, 'CSV', &
          !                   filstat_opt='REPLACE', mode_opt=MNORMAL)
          !     write (this%itrkhdr, '(a,/,a)') &
          !       TRACKHEADERS, &
          !       TRACKTYPES
          !   else
          !     call store_error('OPTIONAL TRACK KEYWORD MUST BE '// &
          !                     'FOLLOWED BY FILEOUT')
          !   end if
          !   found = .true.
          ! case ('TRACKCSV')
          !   call this%parser%GetStringCaps(keyword)
          !   if (keyword == 'FILEOUT') then
          !     ! parse filename
          !     call this%parser%GetString(fname)
          !     ! open CSV file and write headers
          !     this%itrkcsv = getunit()
          !     call openfile(this%itrkcsv, this%iout, fname, 'CSV', &
          !                   filstat_opt='REPLACE')
          !     write (this%iout, fmttrkcsv) trim(adjustl(fname)), this%itrkcsv
          !     write (this%itrkcsv, '(a)') TRACKHEADERS
          !   else
          !     call store_error('OPTIONAL TRACKCSV KEYWORD MUST BE &
          !       &FOLLOWED BY FILEOUT')
          !   end if
          !   found = .true.
        case default
          found = .false.
        end select

        if (.not. found) then
          do ipos = 1, size(this%ocdobj)
            ocdobjptr => this%ocdobj(ipos)
            if (keyword == trim(ocdobjptr%cname)) then
              found = .true.
              exit
            end if
          end do
          if (.not. found) then
            errmsg = "UNKNOWN OC OPTION '"//trim(keyword)//"'."
            call store_error(errmsg)
            call this%parser%StoreErrorUnit()
          end if
          call this%parser%GetRemainingLine(line)
          call ocdobjptr%set_option(line, this%parser%iuactive, this%iout)
        end if
      end do
      write (this%iout, '(1x,a)') 'END OF OC OPTIONS'
    end if
    !
    ! -- return
    return
  end subroutine prt_oc_read_options

end module PrtOcModule
