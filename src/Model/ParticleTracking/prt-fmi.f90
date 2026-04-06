module PrtFmiModule

  use KindModule, only: DP, I4B, LGP
  use ErrorUtilModule, only: pstop
  use ConstantsModule, only: DZERO, LENAUXNAME, LENPACKAGENAME, LENVARNAME
  use SimModule, only: store_error
  use SimVariablesModule, only: errmsg
  use FlowModelInterfaceModule, only: FlowModelInterfaceType
  use BaseDisModule, only: DisBaseType
  use BudgetObjectModule, only: BudgetObjectType
  use BffModule

  implicit none
  private
  public :: PrtFmiType
  public :: fmi_cr

  character(len=LENPACKAGENAME) :: text = '    PRTFMI'

  type, extends(FlowModelInterfaceType) :: PrtFmiType
    private
    integer(I4B), public :: max_faces !< maximum number of 3d cell faces
    real(DP), allocatable, public :: StorageFlows(:) ! cell storage flows array
    type(BffType), pointer, public :: bff => NULL() ! boundary-face flows object

  contains

    procedure :: fmi_ad
    procedure :: fmi_df => prtfmi_df
    procedure :: fmi_ar => prtfmi_ar
    procedure :: fmi_da => prtfmi_da
    procedure, private :: accumulate_flows

  end type PrtFmiType

contains

  !> @brief Create a new PrtFmi object
  subroutine fmi_cr(fmiobj, name_model, input_mempath, inunit, iout)
    ! dummy
    type(PrtFmiType), pointer :: fmiobj
    character(len=*), intent(in) :: name_model
    character(len=*), intent(in) :: input_mempath
    integer(I4B), intent(inout) :: inunit
    integer(I4B), intent(in) :: iout

    ! Create the object
    allocate (fmiobj)

    ! create name and memory path
    call fmiobj%set_names(1, name_model, 'FMI', 'FMI', input_mempath)
    fmiobj%text = text

    ! Allocate scalars
    call fmiobj%allocate_scalars()

    ! Set variables
    fmiobj%inunit = inunit
    fmiobj%iout = iout

    ! Assign dependent variable label
    fmiobj%depvartype = 'TRACKS          '

  end subroutine fmi_cr

  !> @brief Time step advance
  subroutine fmi_ad(this)
    ! modules
    use ConstantsModule, only: DHDRY
    ! dummy
    class(PrtFmiType) :: this
    ! local
    integer(I4B) :: n
    character(len=15) :: nodestr
    character(len=*), parameter :: fmtdry = &
     &"(/1X,'WARNING: DRY CELL ENCOUNTERED AT ',a,';  RESET AS INACTIVE')"
    character(len=*), parameter :: fmtrewet = &
     &"(/1X,'DRY CELL REACTIVATED AT ', a)"

    ! Set flag to indicated that flows are being updated.  For the case where
    ! flows may be reused (only when flows are read from a file) then set
    ! the flag to zero to indicated that flows were not updated
    this%iflowsupdated = 1

    ! If reading flows from a budget file, read the next set of records
    if (this%iubud /= 0) call this%advance_bfr()

    ! If reading heads from a head file, read the next set of records
    if (this%iuhds /= 0) call this%advance_hfr()

    ! If mover flows are being read from file, read the next set of records
    if (this%iumvr /= 0) &
      call this%mvrbudobj%bfr_advance(this%dis, this%iout)

    ! Accumulate flows
    call this%accumulate_flows()

    ! if flow cell is dry, then set this%ibound = 0
    do n = 1, this%dis%nodes
      ! Calculate the ibound-like array that has 0 if saturation
      ! is zero and 1 otherwise
      if (this%gwfsat(n) > DZERO) then
        this%ibdgwfsat0(n) = 1
      else
        this%ibdgwfsat0(n) = 0
      end if

      ! Check if active model cell is inactive for flow
      if (this%ibound(n) > 0) then
        if (this%gwfhead(n) == DHDRY) then
          ! cell should be made inactive
          this%ibound(n) = 0
          call this%dis%noder_to_string(n, nodestr)
          write (this%iout, fmtdry) trim(nodestr)
        end if
      end if

      ! Convert dry model cell to active if flow has rewet
      if (this%ibound(n) == 0) then
        if (this%gwfhead(n) /= DHDRY) then
          ! cell is now wet
          this%ibound(n) = 1
          call this%dis%noder_to_string(n, nodestr)
          write (this%iout, fmtrewet) trim(nodestr)
        end if
      end if
    end do

  end subroutine fmi_ad

  !> @brief Define the flow model interface
  subroutine prtfmi_df(this, dis, idryinactive)
    class(PrtFmiType) :: this
    class(DisBaseType), pointer, intent(in) :: dis
    integer(I4B), intent(in) :: idryinactive

    call this%FlowModelInterfaceType%fmi_df(dis, idryinactive)
    allocate (this%bff)
    call this%bff%bff_df(this%dis, this%nflowpack, this%gwfpackages)

    this%max_faces = this%dis%get_max_npolyverts() + 2
    if (this%max_faces > 32) then
      write (errmsg, '(a,i0,a,i0,a)') &
        'DISV grid contains a cell with ', this%max_faces - 2, &
        ' lateral faces. Cells may have at most 30 lateral faces.'
      call store_error(errmsg)
      call this%parser%StoreErrorUnit()
      return
    end if
    
  end subroutine prtfmi_df

  !> @brief Allocate the package
  !<
  subroutine prtfmi_ar(this, ibound)
    ! -- modules
    ! -- dummy
    class(PrtFmiType) :: this
    integer(I4B), dimension(:), pointer, contiguous :: ibound
    !
    call this%FlowModelInterfaceType%fmi_ar(ibound)
    call this%bff%bff_ar(ibound)
    !
    allocate (this%StorageFlows(this%dis%nodes))

  end subroutine prtfmi_ar

  !> @brief Accumulate flows
  subroutine accumulate_flows(this)
    ! dummy
    class(PrtFmiType) :: this
    ! local

    this%StorageFlows = DZERO
    if (this%igwfstrgss /= 0) &
      this%StorageFlows = this%StorageFlows + this%gwfstrgss
    if (this%igwfstrgsy /= 0) &
      this%StorageFlows = this%StorageFlows + this%gwfstrgsy

    ! Accumulate boundary face flows
    call this%bff%accumulate_flows()
    
  end subroutine accumulate_flows

  !> @brief Deallocate variables
  !<
  subroutine prtfmi_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(PrtFmiType) :: this
    !
    call this%FlowModelInterfaceType%fmi_da()
    !
    ! -- deallocate bff object
    deallocate (this%bff)
    !
  end subroutine prtfmi_da

end module PrtFmiModule
