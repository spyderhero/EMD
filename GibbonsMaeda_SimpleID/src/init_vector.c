#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include <math.h>

// ============================================================
// Minimal Einstein Toolkit initial-condition thorn
// Maxwell vector potential from Poisson (Laplace) solver
// Equation solved (interior):  ∇^2 A_i = 0
// Boundary condition: magnetic monopole–like A_i
// Gauge: Lorenz, Aphi = 0
// ============================================================

// ------------------------------------------------------------
// Boundary condition: monopole-like vector potential
// (only used on outer boundary)
// ------------------------------------------------------------
static inline void monopole_bc(
  const CCTK_REAL x,
  const CCTK_REAL y,
  const CCTK_REAL z,
  const CCTK_REAL Qm,
  CCTK_REAL *Ax,
  CCTK_REAL *Ay,
  CCTK_REAL *Az)
{
  const CCTK_REAL r2 = x*x + y*y + z*z + 1e-14;
  *Ax = -Qm * y / r2;
  *Ay =  Qm * x / r2;
  *Az =  0.0;
}

// ------------------------------------------------------------
// Main initial-data routine
// ------------------------------------------------------------
void Vector_Potential_InitData(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const int imin = 0;
  const int imax = cctk_lsh[0];
  const int jmin = 0;
  const int jmax = cctk_lsh[1];
  const int kmin = 0;
  const int kmax = cctk_lsh[2];

  const CCTK_REAL dx2 = CCTK_DELTA_SPACE(0)*CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy2 = CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz2 = CCTK_DELTA_SPACE(2)*CCTK_DELTA_SPACE(2);

  const CCTK_REAL denom = 2.0*(1.0/dx2 + 1.0/dy2 + 1.0/dz2);

  // ==========================================================
  // (1) Initial guess: zero everywhere
  // ==========================================================
  for (int k = kmin; k <= kmax; ++k)
  for (int j = jmin; j <= jmax; ++j)
  for (int i = imin; i <= imax; ++i)
  {
    const int idx = CCTK_GFINDEX3D(cctkGH,i,j,k);
    Ax[idx] = 0.0;
    Ay[idx] = 0.0;
    Az[idx] = 0.0;
  }

  // ==========================================================
  // (2) Impose boundary condition once
  // ==========================================================
  for (int k = kmin; k <= kmax; ++k)
  for (int j = jmin; j <= jmax; ++j)
  for (int i = imin; i <= imax; ++i)
  {
    if (i==imin || i==imax || j==jmin || j==jmax || k==kmin || k==kmax)
    {
      const int idx = CCTK_GFINDEX3D(cctkGH,i,j,k);
      monopole_bc(x[idx], y[idx], z[idx], par_qm_plus + par_qm_minus,
                  &Ax[idx], &Ay[idx], &Az[idx]);
    }
  }

  // ==========================================================
  // (3) Jacobi iteration: solve ∇^2 A_i = 0
  // ==========================================================
  for (int it = 0; it < poisson_iterations; ++it)
  {
    for (int k = kmin+1; k <= kmax-1; ++k)
    for (int j = jmin+1; j <= jmax-1; ++j)
    for (int i = imin+1; i <= imax-1; ++i)
    {
      const int idx = CCTK_GFINDEX3D(cctkGH,i,j,k);
      const int ip  = CCTK_GFINDEX3D(cctkGH,i+1,j,k);
      const int im  = CCTK_GFINDEX3D(cctkGH,i-1,j,k);
      const int jp  = CCTK_GFINDEX3D(cctkGH,i,j+1,k);
      const int jm  = CCTK_GFINDEX3D(cctkGH,i,j-1,k);
      const int kp  = CCTK_GFINDEX3D(cctkGH,i,j,k+1);
      const int km  = CCTK_GFINDEX3D(cctkGH,i,j,k-1);

      Ax[idx] = ((Ax[ip]+Ax[im])/dx2
                + (Ax[jp]+Ax[jm])/dy2
                + (Ax[kp]+Ax[km])/dz2) / denom;

      Ay[idx] = ((Ay[ip]+Ay[im])/dx2
                + (Ay[jp]+Ay[jm])/dy2
                + (Ay[kp]+Ay[km])/dz2) / denom;

      Az[idx] = ((Az[ip]+Az[im])/dx2
                + (Az[jp]+Az[jm])/dy2
                + (Az[kp]+Az[km])/dz2) / denom;
    }
  }

  if (swap_xz) {
    /* Swap the x and z components of all tensors */
    SWAP (Ax[ind], Az[ind]);
  } /* if swap_xz */
}
