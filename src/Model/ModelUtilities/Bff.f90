module BffModule

  use KindModule, only: DP, I4B, LGP
  use ErrorUtilModule, only: pstop
  use ConstantsModule, only: DZERO, LENAUXNAME, LENPACKAGENAME, LENVARNAME
  use SimModule, only: store_error
  use SimVariablesModule, only: errmsg
  use BaseDisModule, only: DisBaseType
  use BlockParserModule, only: BlockParserType
  use PackageBudgetModule, only: PackageBudgetType

  implicit none
  private
  public :: BffType
  !!public :: bff_cr
  public :: IFLOWFACE_TOP, IFLOWFACE_BOTTOM

  character(len=LENPACKAGENAME) :: text = '       BFF'

  !> @brief IFLOWFACE numbers for top and bottom faces
  enum, bind(c)
    enumerator :: IFLOWFACE_TOP = -1
    enumerator :: IFLOWFACE_BOTTOM = -2
  end enum

  type :: BffType
    private
    integer(I4B), public :: max_faces !< maximum number of 3d cell faces
    real(DP), allocatable, public :: SourceFlows(:) ! cell source flows array
    real(DP), allocatable, public :: SinkFlows(:) ! cell sink flows array
    real(DP), allocatable, public :: BoundaryFlows(:, :) ! cell boundary flows array
    integer(I4B), allocatable, public :: BoundaryFaces(:) ! bitmask of assigned boundary faces
    class(DisBaseType), pointer :: dis => null() !< model discretization object
    integer(I4B), dimension(:), pointer, contiguous :: ibound => null() !< pointer to this model ibound
    integer(I4B), pointer :: nflowpack => null() !< number of GWF flow packages
    type(PackageBudgetType), dimension(:), pointer :: gwfpackages => NULL() !< used to get flows between a package and gwf
    type(BlockParserType) :: parser !< parser object for reading blocks of information

  contains

    procedure :: bff_df
    procedure :: bff_ar
    procedure :: accumulate_flows
    procedure :: mark_boundary_face
    procedure :: is_boundary_face
    procedure :: is_net_out_boundary_face
    procedure, private :: iflowface_to_icellface

  end type BffType

contains

  !> @brief Define the flow model interface
  subroutine bff_df(this, dis, nflowpack, gwfpackages)
    class(BffType) :: this
    class(DisBaseType), pointer, intent(in) :: dis
    integer(I4B), pointer :: nflowpack
    type(PackageBudgetType), dimension(:), allocatable, target :: gwfpackages

    ! -- Store pointers
    this%dis => dis
    this%nflowpack => nflowpack
    this%gwfpackages => gwfpackages

    this%max_faces = this%dis%get_max_npolyverts() + 2
    if (this%max_faces > 32) then
      write (errmsg, '(a,i0,a,i0,a)') &
        'DISV grid contains a cell with ', this%max_faces - 2, &
        ' lateral faces. Cells may have at most 30 lateral faces.'
      call store_error(errmsg)
      call this%parser%StoreErrorUnit()
      return
    end if

    allocate (this%SourceFlows(this%dis%nodes))
    allocate (this%SinkFlows(this%dis%nodes))
    allocate (this%BoundaryFlows(this%dis%nodes, this%max_faces))
    allocate (this%BoundaryFaces(this%dis%nodes))

  end subroutine bff_df
  
  !> @brief Allocate the package
  !<
  subroutine bff_ar(this, ibound)
    ! -- modules
    ! -- dummy
    class(BffType) :: this
    integer(I4B), dimension(:), pointer, contiguous :: ibound
    !
    this%ibound => ibound

  end subroutine bff_ar

  !> @brief Accumulate flows
  subroutine accumulate_flows(this)
    ! dummy
    class(BffType) :: this
    ! local
    integer(I4B) :: j, i, ip, ib
    integer(I4B) :: iflowface, iauxiflowface, icellface
    real(DP) :: qbnd
    character(len=LENAUXNAME) :: auxname
    integer(I4B) :: naux

    this%SourceFlows = DZERO
    this%SinkFlows = DZERO
    this%BoundaryFlows = DZERO
    this%BoundaryFaces = 0
    do ip = 1, this%nflowpack
      iauxiflowface = 0
      naux = this%gwfpackages(ip)%naux
      if (naux > 0) then
        do j = 1, naux
          auxname = this%gwfpackages(ip)%auxname(j)
          if (trim(adjustl(auxname)) == "IFLOWFACE") then
            iauxiflowface = j
            exit
          end if
        end do
      end if
      do ib = 1, this%gwfpackages(ip)%nbound
        i = this%gwfpackages(ip)%nodelist(ib)
        if (i <= 0) cycle
        if (this%ibound(i) <= 0) cycle
        qbnd = this%gwfpackages(ip)%get_flow(ib)
        ! todo, after initial release: default iflowface values for different packages
        iflowface = 0
        icellface = 0
        if (iauxiflowface > 0) then
          iflowface = NINT(this%gwfpackages(ip)%auxvar(iauxiflowface, ib))
          icellface = this%iflowface_to_icellface(iflowface)
        end if
        if (icellface > 0) then
          call this%mark_boundary_face(i, icellface)
          this%BoundaryFlows(i, icellface) = &
            this%BoundaryFlows(i, icellface) + qbnd
        else if (qbnd .gt. DZERO) then
          this%SourceFlows(i) = this%SourceFlows(i) + qbnd
        else if (qbnd .lt. DZERO) then
          this%SinkFlows(i) = this%SinkFlows(i) + qbnd
        end if
      end do
    end do

  end subroutine accumulate_flows

  !> @brief Mark a face as a boundary face.
  subroutine mark_boundary_face(this, ic, icellface)
    class(BffType) :: this
    integer(I4B), intent(in) :: ic !< node number (reduced)
    integer(I4B), intent(in) :: icellface !< cell face number
    ! local
    integer(I4B) :: bit_pos

    if (ic <= 0 .or. ic > this%dis%nodes) then
      print *, 'Invalid cell number: ', ic
      print *, 'Expected a value in range [1, ', this%dis%nodes, ']'
      call pstop(1)
    end if
    if (icellface <= 0) then
      print *, 'Invalid face number: ', icellface
      print *, 'Expected a value in range [1, ', this%max_faces, ']'
      call pstop(1)
    end if
    bit_pos = icellface - 1 ! bit position 0-based
    if (bit_pos < 0 .or. bit_pos > 31) then
      print *, 'Invalid bitmask position: ', bit_pos
      print *, 'Expected a value in range [0, 31]'
      call pstop(1)
    end if
    this%BoundaryFaces(ic) = ibset(this%BoundaryFaces(ic), bit_pos)
  end subroutine mark_boundary_face

  !> @brief Check if a face is assigned to a boundary package.
  function is_boundary_face(this, ic, icellface) result(is_boundary)
    class(BffType) :: this
    integer(I4B), intent(in) :: ic !< node number (reduced)
    integer(I4B), intent(in) :: icellface !< cell face number
    logical(LGP) :: is_boundary
    ! local
    integer(I4B) :: bit_pos

    is_boundary = .false.
    if (ic <= 0 .or. ic > this%dis%nodes) then
      print *, 'Invalid cell number: ', ic
      print *, 'Expected a value in range [1, ', this%dis%nodes, ']'
      call pstop(1)
    end if
    if (icellface <= 0) then
      print *, 'Invalid face number: ', icellface
      print *, 'Expected a value in range [1, ', this%max_faces, ']'
      call pstop(1)
    end if
    bit_pos = icellface - 1 ! bit position 0-based
    if (bit_pos < 0 .or. bit_pos > 31) then
      print *, 'Invalid bitmask position: ', bit_pos
      print *, 'Expected a value in range [0, 31]'
      call pstop(1)
    end if
    is_boundary = btest(this%BoundaryFaces(ic), bit_pos)
  end function is_boundary_face

  !> @brief Check if a face is an assigned boundary with net outflow.
  function is_net_out_boundary_face(this, ic, icellface) &
    result(is_net_out_boundary)
    class(BffType) :: this
    integer(I4B), intent(in) :: ic !< node number (reduced)
    integer(I4B), intent(in) :: icellface !< cell face number
    logical(LGP) :: is_net_out_boundary

    is_net_out_boundary = .false.
    if (.not. this%is_boundary_face(ic, icellface)) return
    if (this%BoundaryFlows(ic, icellface) < DZERO) &
      is_net_out_boundary = .true.
  end function is_net_out_boundary_face

  !> @brief Convert an iflowface number to a cell face number.
  !! Maps bottom (-2) -> max_faces - 1, top (-1) -> max_faces.
  function iflowface_to_icellface(this, iflowface) result(icellface)
    class(BffType), intent(inout) :: this
    integer(I4B), intent(in) :: iflowface
    integer(I4B) :: icellface

    icellface = iflowface
    if (icellface < 0) icellface = icellface + this%max_faces - IFLOWFACE_TOP
  end function iflowface_to_icellface

end module BffModule
