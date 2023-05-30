module TrackDataModule

  use KindModule, only: DP, I4B

  implicit none

  private
  public :: TrackDataType

  ! Structure of arrays to hold particle tracks. Arrays are in long format.
  ! Each particle has 1+ rows representing its track over the model domain.
  ! Particles can be uniquely identified by a combination of column values:
  !   - todo: model ID
  !   - particle release package (PRP) ID
  !   - particle ID (release location ID)
  !   - todo: particle release time
  type :: TrackDataType
    ! integer arrays
    integer(I4B), pointer :: nrows => null() ! track data counter
    integer(I4B), dimension(:), pointer, contiguous :: kper ! stress period
    integer(I4B), dimension(:), pointer, contiguous :: kstp ! time step
    integer(I4B), dimension(:), pointer, contiguous :: iprpid ! PRP ID
    integer(I4B), dimension(:), pointer, contiguous :: ipartid ! particle ID
    integer(I4B), dimension(:), pointer, contiguous :: icellid ! cell ID
    integer(I4B), dimension(:), pointer, contiguous :: istatus ! particle status
    integer(I4B), dimension(:), pointer, contiguous :: ireason ! reason for datum
    ! integer(I4B), dimension(:), pointer, contiguous :: izoneno ! todo zone number
    !   0: release
    !   1: cross cell boundary
    !   2: time step selected
    !   3: termination

    ! double arrays
    real(DP), dimension(:), pointer, contiguous :: x ! current x coordinate
    real(DP), dimension(:), pointer, contiguous :: y ! current y coordinate
    real(DP), dimension(:), pointer, contiguous :: z ! current z coordinate
    real(DP), dimension(:), pointer, contiguous :: t ! current time
    ! real(DP), dimension(:), pointer, contiguous :: trelease ! todo release time

  end type TrackDataType

end module TrackDataModule
