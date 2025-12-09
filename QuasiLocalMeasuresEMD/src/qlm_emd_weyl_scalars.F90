#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine qlm_emd_calc_weyl_scalars (CCTK_ARGUMENTS, hn)
  use qlm_emd_boundary
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: hn

  call qlm_emd_calc_weyl_scalars1 (CCTK_PASS_FTOF, hn)

  call set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_psi0(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_psi1(:,:,hn), -1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_psi2(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_psi3(:,:,hn), -1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_psi4(:,:,hn), +1)
  
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_i(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_j(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_s(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_sdiff(:,:,hn), +1)
  
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_phi00(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_phi11(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_phi01(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_phi12(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_phi10(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_phi21(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_phi02(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_phi22(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_phi20(:,:,hn), +1)
  
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_lambda(:,:,hn), +1)
  
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_lie_n_theta_l(:,:,hn), +1)
 
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_rsc(:,:,hn), +1)
  
end subroutine qlm_emd_calc_weyl_scalars
