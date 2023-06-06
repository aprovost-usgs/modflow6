module TrackDataModule

  use KindModule, only: DP, I4B

  implicit none

  private
  public :: TrackDataType

  ! Structure of arrays to hold particle tracks. Arrays are in long format.
  ! Each particle has 1+ rows representing its track over the model domain.
  ! Particles can be uniquely identified by a combination of column values:
  !   - todo: model ID?
  !   - PRP ID
  !   - particle release location ID
  !   - particle release time (retrieved from partlist)
  type :: TrackDataType
    ! integer arrays
    integer(I4B), pointer :: nrows => null() ! total count of track data
    integer(I4B), dimension(:), pointer, contiguous :: kper ! stress period
    integer(I4B), dimension(:), pointer, contiguous :: kstp ! time step
    integer(I4B), dimension(:), pointer, contiguous :: irpt ! particle ID
    integer(I4B), dimension(:), pointer, contiguous :: iprp ! PRP ID
    integer(I4B), dimension(:), pointer, contiguous :: icell ! cell ID
    integer(I4B), dimension(:), pointer, contiguous :: izone ! todo zone number
    integer(I4B), dimension(:), pointer, contiguous :: istatus ! particle status
    integer(I4B), dimension(:), pointer, contiguous :: ireason ! reason for datum
    ! ireason can take values:
    !   0: release
    !   1: cross cell boundary
    !   2: time step selected (not implemented yet)
    !   3: termination
    !   4: inactive?

    ! double arrays
    real(DP), dimension(:), pointer, contiguous :: x ! current x coordinate
    real(DP), dimension(:), pointer, contiguous :: y ! current y coordinate
    real(DP), dimension(:), pointer, contiguous :: z ! current z coordinate
    real(DP), dimension(:), pointer, contiguous :: t ! current time

  end type TrackDataType

end module TrackDataModule
