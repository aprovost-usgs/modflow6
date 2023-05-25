module test_prt_particle
    use KindModule, only : I4B
    use testdrive, only : error_type, unittest_type, new_unittest, check
    use UtilMiscModule, only : FirstNonBlank, GetLayerRowColumn
    use ParticleModule, only: ParticleListType, resize_particle_list
    use ConstantsModule, only : LINELENGTH
    implicit none
    private
    public :: collect_prt_particle

    contains

    subroutine collect_prt_particle(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
        new_unittest("resize_particle_list", test_resize_particle_list) &
      ]
    end subroutine collect_prt_particle

    subroutine test_resize_particle_list(error)
      type(error_type), allocatable, intent(out) :: error
      type(ParticleListType), pointer :: partlist => null()
      integer(I4B) :: npartmax1 = 10
      integer(I4B) :: npartmax2 = 20

      ! allocate particle arrays
      allocate(partlist)
      allocate(partlist%x(npartmax1))
      allocate(partlist%y(npartmax1))
      allocate(partlist%z(npartmax1))
      allocate(partlist%trelease(npartmax1))
      allocate(partlist%tstop(npartmax1))
      allocate(partlist%ttrack(npartmax1))
      allocate(partlist%istopweaksink(npartmax1))
      allocate(partlist%istopzone(npartmax1))
      allocate(partlist%istatus(npartmax1))
      allocate(partlist%irpt(npartmax1))
      allocate(partlist%iTrackingDomain(npartmax1, 0:4))
      allocate(partlist%iTrackingDomainBoundary(npartmax1, 0:4))

      ! check initial array sizes
      call check(error, size(partlist%x) == npartmax1)
      call check(error, size(partlist%y) == npartmax1)
      call check(error, size(partlist%z) == npartmax1)
      call check(error, size(partlist%trelease) == npartmax1)
      call check(error, size(partlist%tstop) == npartmax1)
      call check(error, size(partlist%ttrack) == npartmax1)
      call check(error, size(partlist%istopweaksink) == npartmax1)
      call check(error, size(partlist%istopzone) == npartmax1)
      call check(error, size(partlist%istatus) == npartmax1)
      call check(error, size(partlist%irpt) == npartmax1)
      call check(error, size(partlist%iTrackingDomain, 1) == npartmax1)
      call check(error, size(partlist%iTrackingDomain, 2) == 5)
      call check(error, size(partlist%iTrackingDomainBoundary, 1) == npartmax1)
      call check(error, size(partlist%iTrackingDomainBoundary, 2) == 5)
      if (allocated(error)) return

      ! resize particle arrays
      call resize_particle_list(partlist, npartmax2)

      ! check that arrays have been resized
      call check(error, size(partlist%x) == npartmax2)
      call check(error, size(partlist%y) == npartmax2)
      call check(error, size(partlist%z) == npartmax2)
      call check(error, size(partlist%trelease) == npartmax2)
      call check(error, size(partlist%tstop) == npartmax2)
      call check(error, size(partlist%ttrack) == npartmax2)
      call check(error, size(partlist%istopweaksink) == npartmax2)
      call check(error, size(partlist%istopzone) == npartmax2)
      call check(error, size(partlist%istatus) == npartmax2)
      call check(error, size(partlist%irpt) == npartmax2)
      call check(error, size(partlist%iTrackingDomain, 1) == npartmax2)
      call check(error, size(partlist%iTrackingDomain, 2) == 5)
      call check(error, size(partlist%iTrackingDomainBoundary, 1) == npartmax2)
      call check(error, size(partlist%iTrackingDomainBoundary, 2) == 5)
      if (allocated(error)) return

    end subroutine test_resize_particle_list
end module test_prt_particle
