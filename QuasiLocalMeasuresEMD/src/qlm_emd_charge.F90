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
  integer :: a, b, c, i, j, m
  CCTK_REAL :: alpha, gg(3,3), gu(3,3), d1_gg(3,3,3)
  CCTK_REAL :: cf1(3,3,3), cf2(3,3,3)
  CCTK_REAL :: E_f(3), A_f(3), dA(3,3), cdA(3,3), B_f(3)
  CCTK_REAL :: dX_dtheta(3), dX_dphi(3), dS(3)
  CCTK_REAL :: sqrtgamma, gamma_det

  CCTK_REAL charge_electric_local, charge_magnetic_local

  character :: msg*1000

  charge_electric_local = 0.0
  charge_magnetic_local = 0.0

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

      ! Compute upper metric.
      gu(1,1) = (gg(2,2) * gg(3,3) - gg(2,3) ** 2     ) / gamma_det
      gu(2,2) = (gg(1,1) * gg(3,3) - gg(1,3) ** 2     ) / gamma_det
      gu(3,3) = (gg(1,1) * gg(2,2) - gg(1,2) ** 2     ) / gamma_det
      gu(1,2) = (gg(1,3) * gg(2,3) - gg(1,2) * gg(3,3)) / gamma_det
      gu(1,3) = (gg(1,2) * gg(2,3) - gg(1,3) * gg(2,2)) / gamma_det
      gu(2,3) = (gg(1,3) * gg(1,2) - gg(2,3) * gg(1,1)) / gamma_det
      gu(2,1) = gu(1,2)
      gu(3,1) = gu(1,3)
      gu(3,2) = gu(2,3)

      d1_gg(1,1,1) = qlm_emd_dgxxx(i,j)
      d1_gg(1,2,1) = qlm_emd_dgxyx(i,j)
      d1_gg(1,3,1) = qlm_emd_dgxzx(i,j)
      d1_gg(1,1,2) = qlm_emd_dgxxy(i,j)
      d1_gg(1,2,2) = qlm_emd_dgxyy(i,j)
      d1_gg(1,3,2) = qlm_emd_dgxzy(i,j)
      d1_gg(1,1,3) = qlm_emd_dgxxz(i,j)
      d1_gg(1,2,3) = qlm_emd_dgxyz(i,j)
      d1_gg(1,3,3) = qlm_emd_dgxzz(i,j)

      d1_gg(2,1,1) = d1_gg(1,2,1)
      d1_gg(2,2,1) = qlm_emd_dgyyx(i,j)
      d1_gg(2,3,1) = qlm_emd_dgyzx(i,j)
      d1_gg(2,1,2) = d1_gg(1,2,2)
      d1_gg(2,2,2) = qlm_emd_dgyyy(i,j)
      d1_gg(2,3,2) = qlm_emd_dgyzy(i,j)
      d1_gg(2,1,3) = d1_gg(1,2,3)
      d1_gg(2,2,3) = qlm_emd_dgyyz(i,j)
      d1_gg(2,3,3) = qlm_emd_dgyzz(i,j)

      d1_gg(3,1,1) = d1_gg(1,3,1)
      d1_gg(3,2,1) = d1_gg(2,3,1)
      d1_gg(3,3,1) = qlm_emd_dgzzx(i,j)
      d1_gg(3,1,2) = d1_gg(1,3,2)
      d1_gg(3,2,2) = d1_gg(2,3,2)
      d1_gg(3,3,2) = qlm_emd_dgzzy(i,j)
      d1_gg(3,1,3) = d1_gg(1,3,3)
      d1_gg(3,2,3) = d1_gg(2,3,3)
      d1_gg(3,3,3) = qlm_emd_dgzzz(i,j)

      ! Compute Christoffel.
      cf1 = 0
      do a = 1, 3
        do b = 1, 3
          do c = b, 3
            cf1(a,b,c) = 0.5d0 * (d1_gg(a,b,c) + d1_gg(a,c,b) - d1_gg(b,c,a))
          end do
        end do
      end do
      cf1(:,2,1) = cf1(:,1,2)
      cf1(:,3,1) = cf1(:,1,3)
      cf1(:,3,2) = cf1(:,2,3)

      cf2 = 0
      do a = 1, 3
        do b = 1, 3
          do c = b, 3
            do m = 1, 3
              cf2(a,b,c) = cf2(a,b,c) + gu(a,m) * cf1(m,b,c)
            end do
          end do
        end do
      end do
      cf2(:,2,1) = cf2(:,1,2)
      cf2(:,3,1) = cf2(:,1,3)
      cf2(:,3,2) = cf2(:,2,3)

      ! Electric field at this point
      E_f(1) = qlm_emd_ex(i,j)
      E_f(2) = qlm_emd_ey(i,j)
      E_f(3) = qlm_emd_ez(i,j)

      ! Vector potential at this point
      A_f(1) = qlm_emd_ax(i,j)
      A_f(2) = qlm_emd_ay(i,j)
      A_f(3) = qlm_emd_az(i,j)

      dA(1,1) = qlm_emd_daxx(i,j)
      dA(1,2) = qlm_emd_daxy(i,j)
      dA(1,3) = qlm_emd_daxz(i,j)
      dA(2,1) = qlm_emd_dayx(i,j)
      dA(2,2) = qlm_emd_dayy(i,j)
      dA(2,3) = qlm_emd_dayz(i,j)
      dA(3,1) = qlm_emd_dazx(i,j)
      dA(3,2) = qlm_emd_dazy(i,j)
      dA(3,3) = qlm_emd_dazz(i,j)

      ! Compite covariant derivatives
      cdA = dA
      do a = 1, 3
        do b = 1, 3
          do m = 1, 3
            cdA(a,b) = cdA(a,b) - cf2(m,a,b) * A_f(m)
        end do
        end do
      end do

      ! Compute magnetic field
      B_f(1) = - alpha * (cdA(3,2) - cdA(2,3))
      B_f(2) = - alpha * (cdA(1,3) - cdA(3,1))
      B_f(3) = - alpha * (cdA(2,1) - cdA(1,2))

      ! Flux contribution
      charge_electric_local = charge_electric_local + (E_f(1)*dS(1) + E_f(2)*dS(2) + E_f(3)*dS(3))
      charge_magnetic_local = charge_magnetic_local + (B_f(1)*dS(1) + B_f(2)*dS(2) + B_f(3)*dS(3))

    end do
  end do

  
  ! Divide by 4π
  qlm_emd_electric_charge(hn) = charge_electric_local / (4.0*pi)
  qlm_emd_magnetic_charge(hn) = - charge_magnetic_local / (4.0*pi)

  write (msg, '("   Electric charge Qe:            ",g14.6)') qlm_emd_electric_charge(hn)
  call CCTK_INFO (msg)
  write (msg, '("   Magnetic charge Qm:            ",g14.6)') qlm_emd_magnetic_charge(hn)
  call CCTK_INFO (msg)

end subroutine qlm_emd_compute_charge
