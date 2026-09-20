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
  CCTK_REAL :: alpha, gg(3,3)
  CCTK_REAL :: E_f(3), dA(3,3), B_f(3)
  CCTK_REAL :: dX_dtheta(3), dX_dphi(3), dS(3)
  CCTK_REAL :: sqrtgamma, gamma_det
  CCTK_REAL :: phi_f, dilaton_factor
  ! Dilatonic electric charge:
  ! Q_EMD = (1/4*pi) int exp(-2*coupling_constant*phi1) E^i dS_i.
  ! Here phi1 is qlm_emd_phi1(i,j), i.e. ScalarBase::phi1 interpolated to the QLM surface.

  CCTK_REAL charge_electric_local, charge_magnetic_local, charge_electric_dilatonic_local


  character :: msg*1000

  charge_electric_local = 0.0d0
  charge_electric_dilatonic_local = 0.0d0  
  charge_magnetic_local = 0.0d0

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
      alpha = qlm_emd_alpha(i,j)
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

      ! Electric field and dilaton at this point
      E_f(1) = qlm_emd_ex(i,j)
      E_f(2) = qlm_emd_ey(i,j)
      E_f(3) = qlm_emd_ez(i,j)
      phi_f = qlm_emd_phi1(i,j)
      dilaton_factor = exp(-2.0d0 * coupling_constant * phi_f)

      dA(1,1) = qlm_emd_daxx(i,j)
      dA(1,2) = qlm_emd_daxy(i,j)
      dA(1,3) = qlm_emd_daxz(i,j)
      dA(2,1) = qlm_emd_dayx(i,j)
      dA(2,2) = qlm_emd_dayy(i,j)
      dA(2,3) = qlm_emd_dayz(i,j)
      dA(3,1) = qlm_emd_dazx(i,j)
      dA(3,2) = qlm_emd_dazy(i,j)
      dA(3,3) = qlm_emd_dazz(i,j)

      ! Compute magnetic field
      B_f(1) = - alpha * (dA(3,2) - dA(2,3))
      B_f(2) = - alpha * (dA(1,3) - dA(3,1))
      B_f(3) = - alpha * (dA(2,1) - dA(1,2))

      ! Flux contribution
      charge_electric_local = charge_electric_local + (E_f(1)*dS(1) + E_f(2)*dS(2) + E_f(3)*dS(3))
      charge_electric_dilatonic_local = charge_electric_dilatonic_local + &
           dilaton_factor * (E_f(1)*dS(1) + E_f(2)*dS(2) + E_f(3)*dS(3))
      charge_magnetic_local = charge_magnetic_local + (B_f(1)*dS(1) + B_f(2)*dS(2) + B_f(3)*dS(3))

    end do
  end do

  ! Divide by 4π
  qlm_emd_electric_charge(hn) = charge_electric_local / (4.0*pi)
  qlm_emd_electric_dilatonic_charge(hn) = charge_electric_dilatonic_local / (4.0*pi)
  qlm_emd_magnetic_charge(hn) = - charge_magnetic_local / (4.0*pi)

  qlm_emd_hn(hn) = dble(hn)
  qlm_emd_surface_index_value(hn) = dble(surface_index(hn))

  write (msg, '("   QLM surface: hn=",i4," surface_index=",i4)') hn, surface_index(hn)
  call CCTK_INFO (msg)
  write (msg, '("   Electric charge Qe:            ",g14.6)') qlm_emd_electric_charge(hn)
  call CCTK_INFO (msg)
  write (msg, '("   Electric and dilatonic charge Q_EMD:    ",g14.6)') qlm_emd_electric_dilatonic_charge(hn)
  call CCTK_INFO (msg)
  write (msg, '("   Magnetic charge Qm:            ",g14.6)') qlm_emd_magnetic_charge(hn)
  call CCTK_INFO (msg)
  

end subroutine qlm_emd_compute_charge
