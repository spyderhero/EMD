
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"



/* Swap two variables */
static inline
void swap (CCTK_REAL * restrict const a, CCTK_REAL * restrict const b)
{
  CCTK_REAL const t = *a; *a=*b; *b=t;
}
#undef SWAP
#define SWAP(a,b) (swap(&(a),&(b)))

/* -------------------------------------------------------------------*/
void GM_GHS (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int imin[3], imax[3];

  for (int d = 0; d < 3; ++ d)
  {
    /*
    imin[d] = 0           + (cctk_bbox[2*d  ] ? 0 : cctk_nghostzones[d]);
    imax[d] = cctk_lsh[d] - (cctk_bbox[2*d+1] ? 0 : cctk_nghostzones[d]);
    */
    imin[d] = 0;
    imax[d] = cctk_lsh[d];
  }


  for (int i = imin[0]; i < imax[0]; ++i) {
    for (int j = imin[1]; j < imax[1]; ++j) {
      for (int k = imin[2]; k < imax[2]; ++k) {

        const int ind = CCTK_GFINDEX3D (cctkGH, i, j, k);

        CCTK_REAL x1, y1, z1;
        x1 = x[ind];
        y1 = y[ind];
        z1 = z[ind];

        /* We implement swapping the x and z coordinates as follows.
           The bulk of the code that performs the actual calculations
           is unchanged.  This code looks only at local variables.
           Before the bulk --i.e., here-- we swap all x and z tensor
           components, and after the code --i.e., at the end of this
           main loop-- we swap everything back.  */
        if (swap_xz) {
          /* Swap the x and z coordinates */
          SWAP (x1, z1);
        }

        CCTK_REAL r_plus
          = sqrt(pow(x1 - par_b, 2) + pow(y1, 2) + pow(z1, 2));
        
        CCTK_REAL rq
          = pow(par_q_plus, 2) * exp(2 * phi1_0) / 2 / par_m_plus; 

        CCTK_REAL F = ( pow(par_m_plus, 2) + 2 * par_m_plus 
                         * (2 * r_plus - rq) + pow(2 * r_plus + rq, 2) )
                         * pow( par_m_plus + 2 * r_plus - rq, 2) / (16 * pow(r_plus, 4)) ;

        gxx[ind] = F;
        gxy[ind] = 0;
        gxz[ind] = 0;
        gyy[ind] = F;
        gyz[ind] = 0;
        gzz[ind] = F;

        kxx[ind] = 0;
        kxy[ind] = 0;
        kxz[ind] = 0;
        kyy[ind] = 0;
        kyz[ind] = 0;
        kzz[ind] = 0;


        /* Scalar terms */
        
        phi1[ind] = phi1_0 + 0.5 * log( pow(par_m_plus + 2 * r_plus - rq, 2) / ( pow(par_m_plus, 2) 
                    + 2 * par_m_plus * (2 * r_plus - rq) + pow(2 * r_plus + rq,2) ));
        Kphi1[ind]  = 0;
        
        phi2[ind] = 0;
        Kphi2[ind]  = 0;
        
        /* EMG terms */

        Zeta[ind]  = 0;

        Ax[ind]    = 0;
        Ay[ind]    = 0;
        Az[ind]    = 0;

        
        // lapse
        if ( CCTK_EQUALS(initial_lapse, "GM_GHS") ) {
          alp[ind] = fabs( (-par_m_plus + 2 * r_plus + rq) / sqrt( pow(par_m_plus, 2)
                     + 2*par_m_plus*(2 *r_plus - rq) + pow(2*r_plus + rq, 2) ) );
        }
        
        Aphi[ind]  = - (Aphi0 - 4 * exp(2*phi1_0) * par_q_plus * r_plus
                     / ( pow(par_m_plus, 2) + 2 * par_m_plus * (2*r_plus -rq) + pow(2*r_plus + rq, 2) ))
                     / alp[ind];
        

        Ex[ind]    = (4 * exp(2 * phi1_0) * par_q_plus * (par_m_plus + 2 * r_plus - rq)
                     * (-par_m_plus + 2 * r_plus + rq) * (x1-par_b)
                     / (r_plus * pow( pow(par_m_plus, 2) + 2*par_m_plus*(2 *r_plus - rq)
                     + pow(2*r_plus + rq, 2), 2))) / alp[ind] /F;

        Ey[ind]    = (4 * exp(2 * phi1_0) * par_q_plus * (par_m_plus + 2 * r_plus - rq)
                     * (-par_m_plus + 2 * r_plus + rq) * y1
                     / (r_plus * pow( pow(par_m_plus, 2) + 2*par_m_plus*(2 *r_plus - rq)
                     + pow(2*r_plus + rq, 2), 2))) / alp[ind] /F;

        Ez[ind]    = (4 * exp(2 * phi1_0) * par_q_plus * (par_m_plus + 2 * r_plus - rq)
                     * (-par_m_plus + 2 * r_plus + rq) * z1
                     / (r_plus * pow( pow(par_m_plus, 2) + 2*par_m_plus*(2 *r_plus - rq)
                     + pow(2*r_plus + rq, 2), 2))) / alp[ind] /F;


        if (swap_xz) {
          /* Swap the x and z components of all tensors */
          SWAP (gxx[ind], gzz[ind]);
          SWAP (gxy[ind], gyz[ind]);
          SWAP (kxx[ind], kzz[ind]);
          SWAP (kxy[ind], kyz[ind]);

          SWAP (Ax[ind], Az[ind]);
          SWAP (Ex[ind], Ez[ind]);

        } /* if swap_xz */


      } /* for k */
    }   /* for j */
  }     /* for i */

}

