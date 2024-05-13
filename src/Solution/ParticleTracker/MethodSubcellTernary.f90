module MethodSubcellTernaryModule
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DSAME, DONE, DHALF
  use ErrorUtilModule, only: pstop
  use GeomUtilModule, only: skew
  use MethodModule, only: MethodType
  use CellModule, only: CellType
  use SubcellModule, only: SubcellType
  use SubcellTriModule, only: SubcellTriType, create_subcell_tri
  use ParticleModule, only: ParticleType, get_particle_id
  use TrackModule, only: TrackFileControlType
  use TernarySolveTrack, only: traverse_triangle, step_analytical, canonical
  use PrtFmiModule, only: PrtFmiType
  use BaseDisModule, only: DisBaseType
  implicit none

  private
  public :: MethodSubcellTernaryType
  public :: create_method_subcell_ternary

  !> @brief Ternary triangular subcell tracking method
  type, extends(MethodType) :: MethodSubcellTernaryType
    integer(I4B), public, pointer :: zeromethod
  contains
    procedure, public :: apply => apply_mst
    procedure, public :: destroy
    procedure, private :: track_subcell
  end type MethodSubcellTernaryType

contains

  !> @brief Create a new ternary subcell-method object
  subroutine create_method_subcell_ternary(method)
    ! -- dummy
    type(MethodSubcellTernaryType), pointer :: method
    ! -- local
    type(SubcellTriType), pointer :: subcell

    allocate (method)
    allocate (method%zeromethod)
    call create_subcell_tri(subcell)
    method%subcell => subcell
    method%type => method%subcell%type
    method%delegates = .false.
    method%zeromethod = 0
  end subroutine create_method_subcell_ternary

  !> @brief Destructor for a ternary subcell-method object
  subroutine destroy(this)
    class(MethodSubcellTernaryType), intent(inout) :: this
    deallocate (this%type)
  end subroutine destroy

  !> @brief Apply the ternary subcell method
  subroutine apply_mst(this, particle, tmax)
    class(MethodSubcellTernaryType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    real(DP), intent(in) :: tmax

    select type (subcell => this%subcell)
    type is (SubcellTriType)
      call this%track_subcell(subcell, particle, tmax)
    end select
  end subroutine apply_mst

  !> @brief Track a particle across a triangular subcell using the ternary method
  subroutine track_subcell(this, subcell, particle, tmax)
    ! dummy
    class(MethodSubcellTernaryType), intent(inout) :: this
    class(SubcellTriType), intent(in) :: subcell
    type(ParticleType), pointer, intent(inout) :: particle
    real(DP), intent(in) :: tmax
    ! local
    integer(I4B) :: exitFace
    logical(LGP) :: lbary
    real(DP) :: x0
    real(DP) :: y0
    real(DP) :: x1
    real(DP) :: y1
    real(DP) :: x2
    real(DP) :: y2
    real(DP) :: v0x
    real(DP) :: v0y
    real(DP) :: v1x
    real(DP) :: v1y
    real(DP) :: v2x
    real(DP) :: v2y
    real(DP) :: xi
    real(DP) :: yi
    real(DP) :: zi
    real(DP) :: zirel
    real(DP) :: ztop
    real(DP) :: zbot
    real(DP) :: dz
    real(DP) :: rxx
    real(DP) :: rxy
    real(DP) :: ryx
    real(DP) :: ryy
    real(DP) :: sxx
    real(DP) :: sxy
    real(DP) :: syy
    real(DP) :: rot(2, 2), res(2), loc(2)
    real(DP) :: alp
    real(DP) :: bet
    real(DP) :: alp0
    real(DP) :: bet0
    real(DP) :: alp1
    real(DP) :: bet1
    real(DP) :: alp2
    real(DP) :: bet2
    real(DP) :: alpi
    real(DP) :: beti
    real(DP) :: gami
    real(DP) :: vzbot
    real(DP) :: vztop
    real(DP) :: vzi
    real(DP) :: vziodz
    real(DP) :: az
    real(DP) :: dtexitz
    real(DP) :: dt
    real(DP) :: t
    real(DP) :: t0
    real(DP) :: dtexitxy
    real(DP) :: texit
    real(DP) :: x
    real(DP) :: y
    real(DP) :: z
    integer(I4B) :: izstatus
    integer(I4B) :: itopbotexit
    integer(I4B) :: ntmax
    integer(I4B) :: nsave
    integer(I4B) :: isolv
    integer(I4B) :: itrifaceenter
    integer(I4B) :: itrifaceexit
    real(DP) :: tol
    real(DP) :: step
    real(DP) :: dtexit
    real(DP) :: alpexit
    real(DP) :: betexit
    integer(I4B) :: reason
    integer(I4B) :: i
    integer(I4B) :: tslice(2)
    real(DP) :: tolnudge, lolimit, hilimit, delta

    lbary = .true. ! todo: remove before initial release
    ntmax = 10000
    nsave = 1 ! needed???
    isolv = this%zeromethod
    tol = 1d-7
    step = 1e-3 ! needed only for euler
    reason = -1

    ! -- Set some local variables for convenience
    xi = particle%x
    yi = particle%y
    zi = particle%z
    x0 = subcell%x0
    y0 = subcell%y0
    x1 = subcell%x1
    y1 = subcell%y1
    x2 = subcell%x2
    y2 = subcell%y2
    v0x = subcell%v0x
    v0y = subcell%v0y
    v1x = subcell%v1x
    v1y = subcell%v1y
    v2x = subcell%v2x
    v2y = subcell%v2y
    zbot = subcell%zbot
    ztop = subcell%ztop
    dz = subcell%dz
    vzbot = subcell%vzbot
    vztop = subcell%vztop

    ! -- Transform coordinates to "canonical" configuration
    call canonical(x0, y0, x1, y1, x2, y2, &
                   v0x, v0y, v1x, v1y, v2x, v2y, &
                   xi, yi, &
                   rxx, rxy, ryx, ryy, &
                   sxx, sxy, syy, &
                   alp0, bet0, alp1, bet1, alp2, bet2, alpi, beti, &
                   lbary)

    ! Nudge particle if necessary so that it begins within the subcell;
    ! no closer than distance tolnudge from any edge of the subcell
    ! todo before initial release: pull into separate routine
    gami = 1d0 - (alpi + beti)
    tolnudge = DSAME
    lolimit = tolnudge
    hilimit = DONE - 2d0 * tolnudge
    ! Check alpha coordinate against lower limit
    if (alpi < lolimit) then
      ! Alpha is too low, so nudge alpha to lower limit; this is a move
      ! parallel to the "alpha axis," which also changes gamma
      alpi = lolimit
      gami = DONE - alpi - beti
      ! Check beta coordinate against lower limit (which in this
      ! case is equivalent to checking gamma coordinate against
      ! upper limit)
      if (beti < lolimit) then
        ! Beta is too low (gamma is too high), so nudge beta to lower limit;
        ! this is a move parallel to the "beta axis," which also changes gamma
        beti = lolimit
        gami = hilimit
      ! Check beta coordinate against upper limit (which in this
      ! case is equivalent to checking gamma coordinate against
      ! lower limit)
      else if (beti > hilimit) then
        ! Beta is too high (gamma is too low), so nudge beta to lower limit;
        ! this is a move parallel to the "beta axis," which also changes gamma
        beti = hilimit
        gami = lolimit
      end if
    end if
    ! Check beta coordinate against lower limit. (If alpha coordinate
    ! was nudged to lower limit, beta and gamma coordinates have also
    ! been adjusted as necessary to place particle within subcell, and
    ! subsequent checks on beta and gamma will evaluate to false, and
    ! no further adjustments will be made.)
    if (beti < lolimit) then
      ! Beta is too low, so nudge beta to lower limit; this is a move
      ! parallel to the "beta axis," which also changes gamma
      beti = lolimit
      gami = DONE - alpi - beti
      ! Check alpha coordinate against lower limit (which in this
      ! case is equivalent to checking gamma coordinate against
      ! upper limit)
      if (alpi < lolimit) then
        ! Alpha is too low (gamma is too high), so nudge alpha to lower limit;
        ! this is a move parallel to the "alpha axis," which also changes gamma
        alpi = lolimit
        gami = hilimit
      ! Check alpha coordinate against upper limit (which in this
      ! case is equivalent to checking gamma coordinate against
      ! lower limit)
      else if (alpi > hilimit) then
        ! Alpha is too high (gamma is too low), so nudge alpha to lower limit;
        ! this is a move parallel to the "alpha axis," which also changes gamma
        alpi = hilimit
        gami = lolimit
      end if
    end if
    ! Check gamma coordinate against lower limit.(If alpha and/or beta
    ! coordinate was nudged to lower limit, gamma coordinate has also
    ! been adjusted as necessary to place particle within subcell, and
    ! subsequent check on gamma will evaluate to false, and no further
    ! adjustment will be made.)
    if (gami < lolimit) then
      ! Gamma is too low, so nudge gamma to lower limit; this is a move
      ! parallel to the "gamma axis," which also changes alpha and beta
      delta = DHALF * (lolimit - gami)
      gami = DSAME
      alpi = alpi - delta
      beti = beti - delta
    end if

    ! -- Do calculations related to analytical z solution, could possibly
    !    be done just once for each cell, todo: profile?
    zirel = (zi - zbot) / dz
    call calculate_dt(vzbot, vztop, dz, zirel, vzi, &
                      az, dtexitz, izstatus, &
                      itopbotexit)
    vziodz = vzi / dz

    ! -- Traverse triangular subcell
    itrifaceenter = particle%iboundary(3) - 1
    if (itrifaceenter .eq. -1) itrifaceenter = 999
    call traverse_triangle(isolv, tol, step, &
                           dtexitxy, alpexit, betexit, &
                           itrifaceenter, itrifaceexit, &
                           rxx, rxy, ryx, ryy, &
                           alp0, bet0, alp1, bet1, alp2, bet2, alpi, beti, &
                           vziodz, az, lbary)

    ! -- Check for no exit face, todo before initial release: terminate with istatus=9
    if ((itopbotexit .eq. 0) .and. (itrifaceexit .eq. 0)) then
      print *, "Subcell with no exit face: particle", get_particle_id(particle), &
        "cell", particle%idomain(2)
      call pstop(1)
    end if

    ! -- Determine (earliest) exit face and corresponding travel time to exit
    if (itopbotexit .eq. 0) then
      ! -- Exits through triangle face first
      exitFace = itrifaceexit
      dtexit = dtexitxy
    else if (itrifaceexit .eq. 0) then
      ! -- Exits through top/bottom first
      exitFace = 45
      dtexit = dtexitz
    else if (dtexitz .lt. dtexitxy) then
      ! -- Exits through top/bottom first
      exitFace = 45
      dtexit = dtexitz
    else
      ! -- Exits through triangle face first
      exitFace = itrifaceexit
      dtexit = dtexitxy
    end if
    if (exitFace .eq. 45) then
      if (itopbotexit .eq. -1) then
        exitFace = 4
      else
        exitFace = 5
      end if
    end if

    ! -- Compute exit time, irrespective of tmax for now
    texit = particle%ttrack + dtexit
    t0 = particle%ttrack

    ! -- Select user tracking times to solve. If this is the first time step
    !    of the simulation, include all times before it begins; if it is the
    !    last time step, include all times after it ends. Otherwise take the
    !    times within the current period and time step only
    !    todo AMP: consider correctness/duplication, maybe pull into separate
    !    routine if duplicative
    ! AMP kluge note: Why are tracking times before first time step begins included?

    call this%tracktimes%try_advance()
    tslice = this%tracktimes%selection
    if (all(tslice > 0)) then
      do i = tslice(1), tslice(2)
        t = this%tracktimes%times(i)
        if (t < particle%ttrack .or. t >= texit .or. t >= tmax) cycle
        dt = t - t0
        !!call step_analytical(dt, alp, bet)
        !!loc = (/alp, bet/)
        !!if (lbary) loc = skew(loc, (/sxx, sxy, syy/), invert=.true.)
        !!rot = reshape((/rxx, rxy, ryx, ryy/), shape(rot))
        !!res = matmul(rot, loc) ! rotate vector
        !!x = res(1) + x0
        !!y = res(2) + y0
        !!if (izstatus .eq. 2) then
        !!  ! -- vz uniformly zero
        !!  z = zi
        !!else if (izstatus .eq. 1) then
        !!  ! -- vz uniform, nonzero
        !!  z = zi + vzi * dt
        !!else
        !!  ! -- vz nonuniform
        !!  z = zbot + (vzi * dexp(az * dt) - vzbot) / az
        !!end if
        ! AMP kluge note: refactored into new subroutine
        call calculate_xyz_position(dt, rxx, rxy, ryx, ryy, sxx, sxy, syy, &
                                    izstatus, lbary, x, y, z)
        particle%x = x
        particle%y = y
        particle%z = z
        particle%ttrack = t
        particle%istatus = 1
        call this%save(particle, reason=5)
      end do
    end if

    ! Compute final time, taking into account tmax, and set final
    ! particle status
    if (texit .gt. tmax) then
      ! -- The computed exit time is greater than the maximum time, so set
      ! -- final time for particle trajectory equal to maximum time.
      t = tmax
      dt = t - t0
      exitFace = 0
      particle%istatus = 1
      particle%advancing = .false.
      reason = 2 ! timestep end
    else
      ! -- The computed exit time is less than or equal to the maximum time,
      ! -- so set final time for particle trajectory equal to exit time.
      t = texit
      dt = dtexit
      reason = 1 ! cell transition
    end if

    ! -- Calculate final local (unscaled) subcell coordinates
    !!call step_analytical(dt, alp, bet)
    !!if (exitFace .eq. 1) then
    !!  bet = 0d0
    !!else if (exitFace .eq. 2) then
    !!  alp = 1d0 - bet
    !!else if (exitFace .eq. 3) then
    !!  alp = 0d0
    !!end if
    !!loc = (/alp, bet/)
    !!if (lbary) loc = skew(loc, (/sxx, sxy, syy/), invert=.true.)
    !!rot = reshape((/rxx, rxy, ryx, ryy/), shape(rot))
    !!res = matmul(rot, loc) ! rotate vector
    !!x = res(1) + x0
    !!y = res(2) + y0
    !!if (exitFace .eq. 4) then
    !!  z = zbot
    !!else if (exitFace .eq. 5) then
    !!  z = ztop
    !!else
    !!  if (izstatus .eq. 2) then
    !!    ! -- vz uniformly zero
    !!    z = zi
    !!  else if (izstatus .eq. 1) then
    !!    ! -- vz uniform, nonzero
    !!    z = zi + vzi * dt
    !!  else
    !!    ! -- vz nonuniform
    !!    z = zbot + (vzi * dexp(az * dt) - vzbot) / az
    !!  end if
    !!end if
    ! AMP kluge note: refactored into new subroutine
    call calculate_xyz_position(dt, rxx, rxy, ryx, ryy, sxx, sxy, syy, &
                                izstatus, lbary, x, y, z, exitface)

    ! -- Set final particle location in local (unscaled) subcell coordinates,
    ! -- final time for particle trajectory, and exit face
    particle%x = x
    particle%y = y
    particle%z = z
    particle%ttrack = t
    particle%iboundary(3) = exitFace

    ! -- Save particle track record
    if (reason > -1) &
      call this%save(particle, reason=reason) ! reason=2: timestep
  end subroutine track_subcell

  !> @brief Do calculations related to analytical z solution
  !!
  !! This subroutine consists partly or entirely of code written by
  !! David W. Pollock of the USGS for MODPATH 7. The authors of the present
  !! code are responsible for its appropriate application in this context
  !! and for any modifications or errors.
  !<
  subroutine calculate_dt(v1, v2, dx, xL, v, dvdx, &
                          dt, status, itopbotexit)
    real(DP) :: v1
    real(DP) :: v2
    real(DP) :: dx
    real(DP) :: xL
    real(DP) :: v
    real(DP) :: dvdx
    real(DP) :: dt
    real(DP) :: v2a
    real(DP) :: v1a
    real(DP) :: dv
    real(DP) :: dva
    real(DP) :: vv
    real(DP) :: vvv
    real(DP) :: zro
    real(DP) :: zrom
    real(DP) :: x
    real(DP) :: tol
    real(DP) :: vr1
    real(DP) :: vr2
    real(DP) :: vr
    real(DP) :: v1v2
    integer(I4B) :: status
    integer(I4B) :: itopbotexit
    logical(LGP) :: noOutflow

    ! Initialize variables
    status = -1
    dt = 1.0d+20
    v2a = v2
    if (v2a .lt. 0d0) v2a = -v2a
    v1a = v1
    if (v1a .lt. 0d0) v1a = -v1a
    dv = v2 - v1
    dva = dv
    if (dva .lt. 0d0) dva = -dva

    ! Check for a uniform zero velocity in this direction.
    ! If so, set status = 2 and return (dt = 1.0d+20).
    tol = 1.0d-15
    if ((v2a .lt. tol) .and. (v1a .lt. tol)) then
      v = 0d0
      dvdx = 0d0
      status = 2
      itopbotexit = 0
      return
    end if

    ! Check for uniform non-zero velocity in this direction.
    ! If so, set compute dt using the constant velocity,
    ! set status = 1 and return.
    vv = v1a
    if (v2a .gt. vv) vv = v2a
    vvv = dva / vv
    if (vvv .lt. 1.0d-4) then
      zro = tol
      zrom = -zro
      v = v1
      x = xL * dx
      if (v1 .gt. zro) then
        dt = (dx - x) / v1
        itopbotexit = -2
      end if
      if (v1 .lt. zrom) then
        dt = -x / v1
        itopbotexit = -1
      end if
      dvdx = 0d0
      status = 1
      return
    end if

    ! Velocity has a linear variation.
    ! Compute velocity corresponding to particle position
    dvdx = dv / dx
    v = (1.0d0 - xL) * v1 + xL * v2

    ! If flow is into the cell from both sides there is no outflow.
    ! In that case, set status = 3 and return
    noOutflow = .true.
    if (v1 .lt. 0d0) noOutflow = .false.
    if (v2 .gt. 0d0) noOutflow = .false.
    if (noOutflow) then
      status = 3
      itopbotexit = 0
      return
    end if

    ! If there is a divide in the cell for this flow direction, check to see if the
    ! particle is located exactly on the divide. If it is, move it very slightly to
    ! get it off the divide. This avoids possible numerical problems related to
    ! stagnation points.
    if ((v1 .le. 0d0) .and. (v2 .ge. 0d0)) then
      if (abs(v) .le. 0d0) then
        v = 1.0d-20
        if (v2 .le. 0d0) v = -v
      end if
    end if

    ! If there is a flow divide, find out what side of the divide the particle
    ! is on and set the value of vr appropriately to reflect that location.
    vr1 = v1 / v
    vr2 = v2 / v
    vr = vr1
    itopbotexit = -1
    if (vr .le. 0d0) then
      vr = vr2
      itopbotexit = -2
    end if

    ! Check if velocity is in the same direction throughout cell (i.e. no flow divide).
    ! Check if product v1*v2 > 0 then the velocity is in the same direction throughout
    ! the cell (i.e. no flow divide). If so, set vr to reflect appropriate direction.
    v1v2 = v1 * v2
    if (v1v2 .gt. 0d0) then
      if (v .gt. 0d0) then
        vr = vr2
        itopbotexit = -2
      end if
      if (v .lt. 0d0) then
        vr = vr1
        itopbotexit = -1
      end if
    end if

    ! Compute travel time to exit face. Return with status = 0
    dt = log(vr) / dvdx
    status = 0
  end subroutine calculate_dt

  subroutine calculate_xyz_position(dt, rxx, rxy, ryx, ryy, sxx, sxy, syy, &
                                    izstatus, lbary, x, y, z, exitFace)
    ! dummy
    real(DP) :: dt
    real(DP) :: rxx
    real(DP) :: rxy
    real(DP) :: ryx
    real(DP) :: ryy
    real(DP) :: sxx
    real(DP) :: sxy
    real(DP) :: syy
    integer(I4B) :: izstatus
    integer(I4B) :: lbary
    real(DP) :: x
    real(DP) :: y
    real(DP) :: z
    integer(I4B), optional :: exitFace
    ! local
    real(DP) :: rot(2, 2), res(2), loc(2)
    real(DP) :: alp
    real(DP) :: bet

    call step_analytical(dt, alp, bet)
    
    if (present(exitface)) then
      ! exitFace is specified, so set corresponding coordinate
      ! exactly
      if (exitFace .eq. 1) then
        bet = 0d0
      else if (exitFace .eq. 2) then
        alp = DONE - bet
      else if (exitFace .eq. 3) then
        alp = 0d0
      else if (exitFace .eq. 4) then
        z = zbot
      else ! exitFace == 5
        z = ztop
      end if
    else
      ! exitFace is not specified, so accept calculated values
      ! of alpha and beta, and calculate z
      if (izstatus .eq. 2) then
        ! -- vz uniformly zero
        z = zi
      else if (izstatus .eq. 1) then
        ! -- vz uniform, nonzero
        z = zi + vzi * dt
      else
        ! -- vz nonuniform
        z = zbot + (vzi * dexp(az * dt) - vzbot) / az
      end if
    end if
    
    ! transform (alp, beta) to (x, y)
    loc = (/alp, bet/)
    if (lbary) loc = skew(loc, (/sxx, sxy, syy/), invert=.true.)
    rot = reshape((/rxx, rxy, ryx, ryy/), shape(rot))
    res = matmul(rot, loc) ! rotate vector
    x = res(1) + x0
    y = res(2) + y0
    
  end subroutine calculate_xyz_position
      

end module MethodSubcellTernaryModule
