#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine qlm_emd_calc_tetrad (CCTK_ARGUMENTS, hn)
  use qlm_emd_boundary
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: hn

  call qlm_emd_calc_tetrad1 (CCTK_PASS_FTOF, hn)
  
  call emd_set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_l0(:,:,hn), +1)
  call emd_set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_l1(:,:,hn), +1)
  call emd_set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_l2(:,:,hn), +1)
  call emd_set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_l3(:,:,hn), +1)
  
  call emd_set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_n0(:,:,hn), +1)
  call emd_set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_n1(:,:,hn), +1)
  call emd_set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_n2(:,:,hn), +1)
  call emd_set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_n3(:,:,hn), +1)
  
  call emd_set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_m0(:,:,hn), +1)
  call emd_set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_m1(:,:,hn), +1)
  call emd_set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_m2(:,:,hn), +1)
  call emd_set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_m3(:,:,hn), +1)
  
end subroutine qlm_emd_calc_tetrad
