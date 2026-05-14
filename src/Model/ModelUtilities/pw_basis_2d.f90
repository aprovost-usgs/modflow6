!====================== pw_basis_2d.f90 ======================
module pw_basis_2d
  use iso_fortran_env, only: dp => real64
  implicit none

  ! Explicit interface to LAPACK DSYEV (symmetric eigendecomposition)
  interface
     subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
       import dp
       character(len=1), intent(in) :: jobz, uplo
       integer,          intent(in) :: n, lda, lwork
       real(dp),         intent(inout) :: a(lda, n)
       real(dp),         intent(out)   :: w(n)
       real(dp),         intent(inout) :: work(*)
       integer,          intent(out)   :: info
     end subroutine dsyev
  end interface

contains

  pure function norm2(v) result(n)
    ! Euclidean norm (2D)
    real(dp), intent(in) :: v(:)
    real(dp) :: n
    n = sqrt(dot_product(v, v))
  end function norm2

!!  subroutine projection_weighted_basis_2d(V, M, S, U, lambda, B, info)
  subroutine projection_weighted_basis_2d(VUNIT, VSCL, M, U, USCL, info)
    ! Compute a projection-weighted basis in 2D:
    ! Inputs:
    !   VUNIT(2,M)    : original 2D UNIT vectors (each column is a unit vector)
    !   VSCL(M)   : magnitude for scaling each unit vector
    !   V(2,M)    : original vectors scaled (VUNIT * VSCL)
    !   M         : number of original vectors
    ! Outputs:
    !   B(2,2)    : scaled basis using projection-weighted magnitudes
    !   info      : LAPACK status (0 = success)
    ! Other:
    !   S(2,2)    : frame operator sum_i v_i v_i^T
    !   U(2,2)    : eigenvectors (columns), ordered by descending eigenvalue
    !   USCL(M)   : magnitude for scaling each eigenvector
    !   lambda(2) : eigenvalues (descending)
    implicit none
    integer,  intent(in)  :: M
!!    real(dp), intent(in)   :: V(2, M)
!!    real(dp), intent(in)  :: V(2, M), VSCL(M)
    real(dp), intent(in)  :: VUNIT(2, M), VSCL(M)
!!    real(dp), intent(out) :: S(2,2), U(2,2), lambda(2), B(2,2)
    real(dp), intent(out) :: U(2,2), USCL(2)
    integer,  intent(out) :: info

    real(dp)               :: V(2, M), S(2,2), lambda(2)
    real(dp)               :: A(2,2), w(2)
    real(dp), allocatable  :: work(:)
    integer                :: lwork, i, j
    real(dp)               :: numer, denom, tiny
    
    ! Build S = sum v_i v_i^T
    S = 0.0_dp
    do i = 1, M
       V(:, i) = VUNIT(:, i) * VSCL(i)
       S(1,1) = S(1,1) + V(1,i)*V(1,i)
       S(1,2) = S(1,2) + V(1,i)*V(2,i)
       S(2,1) = S(2,1) + V(2,i)*V(1,i)
       S(2,2) = S(2,2) + V(2,i)*V(2,i)
    end do

    ! Copy S to A for LAPACK (DSYEV overwrites input)
    A = S

    ! Workspace query
    lwork = -1
    allocate(work(1))
    call dsyev('V','U', 2, A, 2, w, work, lwork, info)
    if (info /= 0) then
       deallocate(work)
!!       U = 0.0_dp; lambda = 0.0_dp; B = 0.0_dp
       U = 0.0_dp; lambda = 0.0_dp
       return
    end if
    lwork = int(work(1))
    deallocate(work)
    allocate(work(lwork))

    ! Actual eigendecomposition (DSYEV returns ascending eigenvalues)
    call dsyev('V','U', 2, A, 2, w, work, lwork, info)
    deallocate(work)
    if (info /= 0) then
!!       U = 0.0_dp; lambda = 0.0_dp; B = 0.0_dp
       U = 0.0_dp; lambda = 0.0_dp
       return
    end if

    ! Reorder eigenpairs to descending order
    call reorder_descending_2d(w, A, lambda, U)

    ! Projection-weighted magnitudes
    tiny = 100.0_dp * epsilon(1.0_dp)
    do j = 1, 2
       numer = 0.0_dp
       denom = 0.0_dp
       do i = 1, M
!!          denom = denom + abs(dot_product(U(:,j), V(:,i)))
          denom = denom + abs(dot_product(U(:,j), VUNIT(:,i)))
!!          numer = numer + abs(dot_product(U(:,j), V(:,i))) * norm2(V(:,i))
!!          numer = numer + dot_product(U(:,j), V(:,i)) * norm2(V(:,i))
          numer = numer + dot_product(U(:,j), V(:,i))
!!          numer = numer + abs(dot_product(U(:,j), V(:,i))) * VSCL(i)
       end do
       if (denom > tiny) then
!!          B(:,j) = (numer / denom) * U(:,j)
          USCL(j) = numer / denom
       else
!!          B(:,j) = 0.0_dp
          USCL(j) = 0.0_dp
       end if
    end do

  contains

    subroutine reorder_descending_2d(w_in, A_in, w_out, U_out)
      ! Sort eigenpairs by descending eigenvalue for 2x2 case.
      real(dp), intent(in)  :: w_in(2), A_in(2,2)
      real(dp), intent(out) :: w_out(2), U_out(2,2)
      integer :: p1, p2

      if (w_in(1) >= w_in(2)) then
         p1 = 1; p2 = 2
      else
         p1 = 2; p2 = 1
      end if

      w_out(1)   = w_in(p1)
      w_out(2)   = w_in(p2)
      U_out(:,1) = A_in(:,p1)
      U_out(:,2) = A_in(:,p2)
    end subroutine reorder_descending_2d

  end subroutine projection_weighted_basis_2d

end module pw_basis_2d
!====================== end module ======================