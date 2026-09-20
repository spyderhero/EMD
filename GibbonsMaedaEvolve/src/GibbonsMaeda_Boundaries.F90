#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine GibbonsMaeda_Boundaries( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT ierr

  CCTK_INT, parameter :: one = 1

  ! The outgoing (radiative) boundary conditions are being handled from calls to
  ! the NewRad infraestructure. Here we register all BCs as 'none', which
  ! enforces all the symmetry BCs.

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,      &
       "ProcaBase::Ei", "none")
  if (ierr < 0)                                                            &
       call CCTK_WARN(0, "Failed to register BC for ProcaBase::Ei!")

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,      &
       "ProcaBase::Ai", "none")
  if (ierr < 0)                                                            &
       call CCTK_WARN(0, "Failed to register BC for ProcaBase::Ai!")

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, one, -one,        &
       "ProcaBase::Aphi", "none")
  if (ierr < 0)                                                            &
       call CCTK_WARN(0, "Failed to register BC for ProcaBase::Aphi!")

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, one, -one,        &
       "ProcaBase::Zeta", "none")
  if (ierr < 0)                                                            &
       call CCTK_WARN(0, "Failed to register BC for ProcaBase::Zeta!")
  
  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, one, -one,        &
        "ScalarBase::phi1", "none")
  if (ierr < 0)                                                            &
        call CCTK_WARN(0, "Failed to register BC for ScalarBase::phi1!")

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, one, -one,        &
        "ScalarBase::phi2", "none")
  if (ierr < 0)                                                            &
        call CCTK_WARN(0, "Failed to register BC for ScalarBase::phi2!")

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, one, -one,        &
        "ScalarBase::Kphi1", "none")
  if (ierr < 0)                                                            &
        call CCTK_WARN(0, "Failed to register BC for ScalarBase::Kphi1!")

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, one, -one,        &
        "ScalarBase::Kphi2", "none")
  if (ierr < 0)                                                            &
        call CCTK_WARN(0, "Failed to register BC for ScalarBase::Kphi2!")
  
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,        &
        "GibbonsMaedaEvolve::EMD_jr_phi_gfs", "none")
  if (ierr < 0)                                                            &
        call CCTK_WARN(0, "Failed to register BC for GibbonsMaedaEvolve::EMD_jr_phi_gfs!")
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,        &
        "GibbonsMaedaEvolve::EMD_Sir_phi_gfs", "none")
  if (ierr < 0)                                                            &
        call CCTK_WARN(0, "Failed to register BC for GibbonsMaedaEvolve::EMD_Sir_phi_gfs!")
        
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,        &
        "GibbonsMaedaEvolve::EMD_jr_F_gfs", "none")
  if (ierr < 0)                                                            &
        call CCTK_WARN(0, "Failed to register BC for GibbonsMaedaEvolve::EMD_jr_F_gfs!")
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,        &
        "GibbonsMaedaEvolve::EMD_Sir_F_gfs", "none")
  if (ierr < 0)                                                            &
        call CCTK_WARN(0, "Failed to register BC for GibbonsMaedaEvolve::EMD_Sir_F_gfs!")

end subroutine GibbonsMaeda_Boundaries
