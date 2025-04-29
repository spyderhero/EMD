#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

subroutine GibbonsMaeda_InitSymBound( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT ierr

  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "GibbonsMaedaEvolve::rhs_Ex" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "GibbonsMaedaEvolve::rhs_Ey" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "GibbonsMaedaEvolve::rhs_Ez" )

  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "GibbonsMaedaEvolve::rhs_Ax" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "GibbonsMaedaEvolve::rhs_Ay" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "GibbonsMaedaEvolve::rhs_Az" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "GibbonsMaedaEvolve::rhs_Zeta" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "GibbonsMaedaEvolve::rhs_Aphi" )
  
  call SetCartSymVN( ierr, cctkGH, (/ -1, -1, -1/), "GibbonsMaedaEvolve::rhs_phi1" )
  call SetCartSymVN( ierr, cctkGH, (/ -1, -1, -1/), "GibbonsMaedaEvolve::rhs_phi2" )
  call SetCartSymVN( ierr, cctkGH, (/ -1, -1, -1/), "GibbonsMaedaEvolve::rhs_Kphi1" )
  call SetCartSymVN( ierr, cctkGH, (/ -1, -1, -1/), "GibbonsMaedaEvolve::rhs_Kphi2" )
  
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "GibbonsMaedaEvolve::jrPhi_gf" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "GibbonsMaedaEvolve::jrF_gf" )
  
  call SetCartSymVN( ierr, cctkGH, (/ -1, 1, 1/), "GibbonsMaedaEvolve::SxrPhi_gf" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, -1, 1/), "GibbonsMaedaEvolve::SyrPhi_gf" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, -1/), "GibbonsMaedaEvolve::SzrPhi_gf" )
  
  call SetCartSymVN( ierr, cctkGH, (/ -1, 1, 1/), "GibbonsMaedaEvolve::SxrF_gf" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, -1, 1/), "GibbonsMaedaEvolve::SyrF_gf" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, -1/), "GibbonsMaedaEvolve::SzrF_gf" )

end subroutine GibbonsMaeda_InitSymBound
