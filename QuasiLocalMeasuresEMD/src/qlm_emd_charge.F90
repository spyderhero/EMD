#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"


subroutine qlm_emd_compute_charge (CCTK_ARGUMENTS, hn)
  use cctk
  use constants
  use qlm_emd_derivs
  use qlm_emd_variables
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: hn
  integer :: i, j
  CCTK_REAL :: gg(3,3)
  CCTK_REAL :: calc_Ex, calc_Ey, calc_Ez
  CCTK_REAL :: dX_dtheta(3), dX_dphi(3), dS(3)
  CCTK_REAL :: sqrtgamma, gamma_det

  CCTK_REAL charge_local

  character :: msg*1000

  charge_local = 0.0

  do j = 1, qlm_emd_nphi(hn)-1
    do i = 1, qlm_emd_ntheta(hn)-1

      ! Tangent vectors in 3-space
      dX_dtheta(1) = qlm_emd_x(i+1,j,hn) - qlm_emd_x(i,j,hn)
      dX_dtheta(2) = qlm_emd_y(i+1,j,hn) - qlm_emd_y(i,j,hn)
      dX_dtheta(3) = qlm_emd_z(i+1,j,hn) - qlm_emd_z(i,j,hn)

      dX_dphi(1) = qlm_emd_x(i,j+1,hn) - qlm_emd_x(i,j,hn)
      dX_dphi(2) = qlm_emd_y(i,j+1,hn) - qlm_emd_y(i,j,hn)
      dX_dphi(3) = qlm_emd_z(i,j+1,hn) - qlm_emd_z(i,j,hn)

      ! Cross product → surface element vector
      dS(1) = dX_dtheta(2)*dX_dphi(3) - dX_dtheta(3)*dX_dphi(2)
      dS(2) = dX_dtheta(3)*dX_dphi(1) - dX_dtheta(1)*dX_dphi(3)
      dS(3) = dX_dtheta(1)*dX_dphi(2) - dX_dtheta(2)*dX_dphi(1)

      ! Compute sqrt(gamma) from qlm_emd_gxx etc.
      gg(1,1) = qlm_emd_gxx(i,j)
      gg(1,2) = qlm_emd_gxy(i,j)
      gg(1,3) = qlm_emd_gxz(i,j)
      gg(2,2) = qlm_emd_gyy(i,j)
      gg(2,3) = qlm_emd_gyz(i,j)
      gg(3,3) = qlm_emd_gzz(i,j)
      gg(2,1) = gg(1,2)
      gg(3,1) = gg(1,3)
      gg(3,2) = gg(2,3)
      gamma_det = gg(1,1)*gg(2,2)*gg(3,3) &
                  + 2.0*gg(1,2)*gg(2,3)*gg(1,3) &
                  - gg(1,1)*gg(2,3)**2 &
                  - gg(2,2)*gg(1,3)**2 &
                  - gg(3,3)*gg(1,2)**2

      sqrtgamma = sqrt(gamma_det)

      ! Multiply surface element by sqrt(gamma)
      dS = dS * sqrtgamma

      ! Electric field at this point
      calc_Ex = qlm_emd_ex(i,j)
      calc_Ey = qlm_emd_ey(i,j)
      calc_Ez = qlm_emd_ez(i,j)

      ! Flux contribution
      charge_local = charge_local + (calc_Ex*dS(1) + calc_Ey*dS(2) + calc_Ez*dS(3))

    end do
  end do

  
  ! Divide by 4π
  qlm_emd_electric_charge(hn) = charge_local / (4.0*pi)

  write (msg, '("   Electric charge Qe:            ",g14.6)') qlm_emd_electric_charge(hn)
  call CCTK_INFO (msg)

end subroutine qlm_emd_compute_charge
