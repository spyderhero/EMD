#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"


subroutine qlm_compute_charge (CCTK_ARGUMENTS, hn)
  implicit none
  integer :: hn
  integer :: i, j
  CCTK_REAL :: Ex, Ey, Ez
  CCTK_REAL :: dX_dtheta(3), dX_dphi(3), dS(3)
  CCTK_REAL :: sqrtgamma, gamma_det

  CCTK_REAL charge_local

  charge_local = 0.0

  do j = 1, qlm_nphi(hn)-1
    do i = 1, qlm_ntheta(hn)-1

      ! Tangent vectors in 3-space
      dX_dtheta(1) = qlm_x(i+1,j,hn) - qlm_x(i,j,hn)
      dX_dtheta(2) = qlm_y(i+1,j,hn) - qlm_y(i,j,hn)
      dX_dtheta(3) = qlm_z(i+1,j,hn) - qlm_z(i,j,hn)

      dX_dphi(1) = qlm_x(i,j+1,hn) - qlm_x(i,j,hn)
      dX_dphi(2) = qlm_y(i,j+1,hn) - qlm_y(i,j,hn)
      dX_dphi(3) = qlm_z(i,j+1,hn) - qlm_z(i,j,hn)

      ! Cross product → surface element vector
      dS(1) = dX_dtheta(2)*dX_dphi(3) - dX_dtheta(3)*dX_dphi(2)
      dS(2) = dX_dtheta(3)*dX_dphi(1) - dX_dtheta(1)*dX_dphi(3)
      dS(3) = dX_dtheta(1)*dX_dphi(2) - dX_dtheta(2)*dX_dphi(1)

      ! Compute sqrt(gamma) from qlm_gxx etc.
      gamma_det = qlm_gxx(i,j,hn)*qlm_gyy(i,j,hn)*qlm_gzz(i,j,hn) &
                  + 2.0*qlm_gxy(i,j,hn)*qlm_gyz(i,j,hn)*qlm_gxz(i,j,hn) &
                  - qlm_gxx(i,j,hn)*qlm_gyz(i,j,hn)**2 &
                  - qlm_gyy(i,j,hn)*qlm_gxz(i,j,hn)**2 &
                  - qlm_gzz(i,j,hn)*qlm_gxy(i,j,hn)**2

      sqrtgamma = sqrt(gamma_det)

      ! Multiply surface element by sqrt(gamma)
      dS = dS * sqrtgamma

      ! Electric field at this point
      Ex = qlm_ex(i,j,hn)
      Ey = qlm_ey(i,j,hn)
      Ez = qlm_ez(i,j,hn)

      ! Flux contribution
      charge_local = charge_local + (Ex*dS(1) + Ey*dS(2) + Ez*dS(3))

    end do
  end do

  ! Divide by 4π
  qlm_charge(hn) = charge_local / (4.0*acos(-1.0))

end subroutine qlm_compute_charge
