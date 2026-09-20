
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void GibbonsMaeda_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr = 0, group, rhs, var;

  // register evolution and rhs gridfunction groups with MoL

  /* metric and extrinsic curvature */
  group = CCTK_GroupIndex("ADMBase::lapse");
  ierr += MoLRegisterSaveAndRestoreGroup(group);
  group = CCTK_GroupIndex("ADMBase::shift");
  ierr += MoLRegisterSaveAndRestoreGroup(group);
  group = CCTK_GroupIndex("ADMBase::metric");
  ierr += MoLRegisterSaveAndRestoreGroup(group);
  group = CCTK_GroupIndex("ADMBase::curv");
  ierr += MoLRegisterSaveAndRestoreGroup(group);

  /* Ei and rhs_Ei */
  group = CCTK_GroupIndex("ProcaBase::Ei");
  rhs   = CCTK_GroupIndex("GibbonsMaedaEvolve::rhs_Ei");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  /* Ai and rhs_Ai */
  group = CCTK_GroupIndex("ProcaBase::Ai");
  rhs   = CCTK_GroupIndex("GibbonsMaedaEvolve::rhs_Ai");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  /* Aphi and rhs_Aphi */
  var   = CCTK_VarIndex("ProcaBase::Aphi");
  rhs   = CCTK_VarIndex("GibbonsMaedaEvolve::rhs_Aphi");
  ierr += MoLRegisterEvolved(var, rhs);

  /* Zeta and rhs_Zeta */
  var   = CCTK_VarIndex("ProcaBase::Zeta");
  rhs   = CCTK_VarIndex("GibbonsMaedaEvolve::rhs_Zeta");
  ierr += MoLRegisterEvolved(var, rhs);
  
  /* phi1 and rhs_phi1 */
  var   = CCTK_VarIndex("ScalarBase::phi1");
  rhs   = CCTK_VarIndex("GibbonsMaedaEvolve::rhs_phi1");
  ierr += MoLRegisterEvolved(var, rhs);

  /* phi2 and rhs_phi2 */
  var   = CCTK_VarIndex("ScalarBase::phi2");
  rhs   = CCTK_VarIndex("GibbonsMaedaEvolve::rhs_phi2");
  ierr += MoLRegisterEvolved(var, rhs);

  /* Kphi1 and rhs_Kphi1 */
  var   = CCTK_VarIndex("ScalarBase::Kphi1");
  rhs   = CCTK_VarIndex("GibbonsMaedaEvolve::rhs_Kphi1");
  ierr += MoLRegisterEvolved(var, rhs);

  /* Kphi2 and rhs_Kphi2 */
  var   = CCTK_VarIndex("ScalarBase::Kphi2");
  rhs   = CCTK_VarIndex("GibbonsMaedaEvolve::rhs_Kphi2");
  ierr += MoLRegisterEvolved(var, rhs);

  if (ierr) CCTK_ERROR("Problems registering with MoL");

}
