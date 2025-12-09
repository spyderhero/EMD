#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine qlm_emd_calc_newman_penrose (CCTK_ARGUMENTS, hn)
  use qlm_emd_boundary
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: hn

  call qlm_emd_calc_newman_penrose1 (CCTK_PASS_FTOF, hn)

  call set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_npkappa  (:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_nptau    (:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_npsigma  (:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_nprho    (:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_npepsilon(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_npgamma  (:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_npbeta   (:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_npalpha  (:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_nppi     (:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_npnu     (:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_npmu     (:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_emd_nplambda (:,:,hn), +1)
  
end subroutine qlm_emd_calc_newman_penrose
