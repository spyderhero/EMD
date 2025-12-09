#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine qlm_emd_calc_twometric (CCTK_ARGUMENTS, hn)
  use adm_metric
  use cctk
  use qlm_emd_boundary
  use qlm_emd_derivs
  use qlm_emd_variables
  use ricci
  use ricci2
  use tensor
  use tensor2
  use tensor4
  
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: hn
  
  CCTK_REAL :: gg(3,3), dgg(3,3,3)
  CCTK_REAL :: ee(3,2), dee(3,2,2)
  CCTK_REAL :: qq(2,2), dqq(2,2,2)
 
  CCTK_REAL :: delta_space(2)
  
  integer   :: i, j
  integer   :: a, b, c, d, e, f
  
  if (veryverbose/=0) then
     call CCTK_INFO ("Calculating two-metric")
  end if
  
  delta_space(:) = (/ qlm_emd_delta_theta(hn), qlm_emd_delta_phi(hn) /)
  
  ! Calculate the two-metric
  do j = 1+qlm_emd_nghostsphi(hn), qlm_emd_nphi(hn)-qlm_emd_nghostsphi(hn)
     do i = 1+qlm_emd_nghoststheta(hn), qlm_emd_ntheta(hn)-qlm_emd_nghoststheta(hn)
        
        gg(1,1) = qlm_emd_gxx(i,j)
        gg(1,2) = qlm_emd_gxy(i,j)
        gg(1,3) = qlm_emd_gxz(i,j)
        gg(2,2) = qlm_emd_gyy(i,j)
        gg(2,3) = qlm_emd_gyz(i,j)
        gg(3,3) = qlm_emd_gzz(i,j)
        gg(2,1) = gg(1,2)
        gg(3,1) = gg(1,3)
        gg(3,2) = gg(2,3)
        
        dgg(1,1,1) = qlm_emd_dgxxx(i,j)
        dgg(1,2,1) = qlm_emd_dgxyx(i,j)
        dgg(1,3,1) = qlm_emd_dgxzx(i,j)
        dgg(2,2,1) = qlm_emd_dgyyx(i,j)
        dgg(2,3,1) = qlm_emd_dgyzx(i,j)
        dgg(3,3,1) = qlm_emd_dgzzx(i,j)
        dgg(1,1,2) = qlm_emd_dgxxy(i,j)
        dgg(1,2,2) = qlm_emd_dgxyy(i,j)
        dgg(1,3,2) = qlm_emd_dgxzy(i,j)
        dgg(2,2,2) = qlm_emd_dgyyy(i,j)
        dgg(2,3,2) = qlm_emd_dgyzy(i,j)
        dgg(3,3,2) = qlm_emd_dgzzy(i,j)
        dgg(1,1,3) = qlm_emd_dgxxz(i,j)
        dgg(1,2,3) = qlm_emd_dgxyz(i,j)
        dgg(1,3,3) = qlm_emd_dgxzz(i,j)
        dgg(2,2,3) = qlm_emd_dgyyz(i,j)
        dgg(2,3,3) = qlm_emd_dgyzz(i,j)
        dgg(3,3,3) = qlm_emd_dgzzz(i,j)
        dgg(2,1,:) = dgg(1,2,:)
        dgg(3,1,:) = dgg(1,3,:)
        dgg(3,2,:) = dgg(2,3,:)
        
        ee(1,1:2) = deriv (qlm_emd_x(:,:,hn), i, j, delta_space)
        ee(2,1:2) = deriv (qlm_emd_y(:,:,hn), i, j, delta_space)
        ee(3,1:2) = deriv (qlm_emd_z(:,:,hn), i, j, delta_space)
        
        dee(1,1:2,1:2) = deriv2 (qlm_emd_x(:,:,hn), i, j, delta_space)
        dee(2,1:2,1:2) = deriv2 (qlm_emd_y(:,:,hn), i, j, delta_space)
        dee(3,1:2,1:2) = deriv2 (qlm_emd_z(:,:,hn), i, j, delta_space)
        
        do a=1,2
           do b=1,2
              qq(a,b) = 0
              do c=1,3
                 do d=1,3
                    qq(a,b) = qq(a,b) + gg(c,d) * ee(c,a) * ee(d,b)
                 end do
              end do
           end do
        end do
        
        do a=1,2
           do b=1,2
              do c=1,2
                 dqq(a,b,c) = 0
                 do d=1,3
                    do e=1,3
                       do f=1,3
                          dqq(a,b,c) = dqq(a,b,c) + dgg(d,e,f) * ee(d,a) * ee(e,b) * ee(f,c)
                       end do
                       dqq(a,b,c) = dqq(a,b,c) + gg(d,e) * dee(d,a,c) * ee(e,b)
                       dqq(a,b,c) = dqq(a,b,c) + gg(d,e) * ee(d,a) * dee(e,b,c)
                    end do
                 end do
              end do
           end do
        end do
        
        ! Could also calculate this as:
        !    q^ab = m^a mbar^b + mbar^a m^b
        qlm_emd_qtt(i,j,hn) = qq(1,1)
        qlm_emd_qtp(i,j,hn) = qq(1,2)
        qlm_emd_qpp(i,j,hn) = qq(2,2)
        
        qlm_emd_dqttt(i,j) = dqq(1,1,1)
        qlm_emd_dqtpt(i,j) = dqq(1,2,1)
        qlm_emd_dqppt(i,j) = dqq(2,2,1)
        qlm_emd_dqttp(i,j) = dqq(1,1,2)
        qlm_emd_dqtpp(i,j) = dqq(1,2,2)
        qlm_emd_dqppp(i,j) = dqq(2,2,2)
        
     end do
  end do
  
  call emd_set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_qtt(:,:,hn), +1)
  call emd_set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_qtp(:,:,hn), +1)
  call emd_set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_qpp(:,:,hn), +1)
  
  call emd_set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_dqttt(:,:), -1)
  call emd_set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_dqtpt(:,:), -1)
  call emd_set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_dqppt(:,:), -1)
  call emd_set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_dqttp(:,:), -1)
  call emd_set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_dqtpp(:,:), -1)
  call emd_set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_dqppp(:,:), -1)
  
end subroutine qlm_emd_calc_twometric
