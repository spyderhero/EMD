/* TwoPunctures:  File  "TwoPunctures.h"*/

#define StencilSize 19
#define N_PlaneRelax 1
#define NRELAX 200
#define Step_Relax 1

typedef struct DERIVS
{
  CCTK_REAL *d0, *d1, *d2, *d3, *d11, *d12, *d13, *d22, *d23, *d33;
} derivs;

/*
Files of "TwoPunctures":
	TwoPunctures.c
	FuncAndJacobian.c
	CoordTransf.c
	Equations.c
	Newton.c
	utilities.c (see utilities.h)
**************************
*/

/* Routines in  "TwoPunctures.c"*/
CCTK_REAL TestSolution (CCTK_REAL A, CCTK_REAL B, CCTK_REAL X, CCTK_REAL R, CCTK_REAL phi);
void TestVector_w (CCTK_REAL *par, int nvar, int n1, int n2, int n3, CCTK_REAL *w);

/* Routines in  "FuncAndJacobian.c"*/
int GibbonsMaeda_Index (int ivar, int i, int j, int k, int nvar, int n1, int n2, int n3);
void GibbonsMaeda_allocate_derivs (derivs * v, int n);
void GibbonsMaeda_free_derivs (derivs * v, int n);
void GibbonsMaeda_Derivatives_AB3 (int nvar, int n1, int n2, int n3, derivs v);
void GibbonsMaeda_F_of_v (CCTK_POINTER_TO_CONST cctkGH,
       int nvar, int n1, int n2, int n3, derivs v,
	     CCTK_REAL *F, derivs u);
void GibbonsMaeda_J_times_dv (int nvar, int n1, int n2, int n3, derivs dv,
		 CCTK_REAL *Jdv, derivs u);
void GibbonsMaeda_JFD_times_dv (int i, int j, int k, int nvar, int n1,
		   int n2, int n3, derivs dv, derivs u, CCTK_REAL *values);
void GibbonsMaeda_SetMatrix_JFD (int nvar, int n1, int n2, int n3,
		    derivs u, int *ncols, int **cols, CCTK_REAL **Matrix);
CCTK_REAL GibbonsMaeda_PunctEvalAtArbitPosition (CCTK_REAL *v, int ivar, CCTK_REAL A, CCTK_REAL B, CCTK_REAL phi,
				 int nvar, int n1, int n2, int n3);
void GibbonsMaeda_calculate_derivs (int i, int j, int k, int ivar, int nvar, int n1,
		       int n2, int n3, derivs v, derivs vv);
CCTK_REAL GibbonsMaeda_interpol (CCTK_REAL a, CCTK_REAL b, CCTK_REAL c, derivs v);
CCTK_REAL GibbonsMaeda_PunctTaylorExpandAtArbitPosition (int ivar, int nvar, int n1,
                                         int n2, int n3, derivs v, CCTK_REAL x,
                                         CCTK_REAL y, CCTK_REAL z);
CCTK_REAL GibbonsMaeda_PunctIntPolAtArbitPosition (int ivar, int nvar, int n1,
				   int n2, int n3, derivs v, CCTK_REAL x,
				   CCTK_REAL y, CCTK_REAL z);
void GibbonsMaeda_SpecCoef(int n1, int n2, int n3, int ivar, CCTK_REAL *v, CCTK_REAL *cf);
CCTK_REAL GibbonsMaeda_PunctEvalAtArbitPositionFast (CCTK_REAL *v, int ivar, CCTK_REAL A, CCTK_REAL B, CCTK_REAL phi,
                                 int nvar, int n1, int n2, int n3);
CCTK_REAL GibbonsMaeda_PunctIntPolAtArbitPositionFast (int ivar, int nvar, int n1,
                                   int n2, int n3, derivs v, CCTK_REAL x,
                                   CCTK_REAL y, CCTK_REAL z);


/* Routines in  "CoordTransf.c"*/
void GibbonsMaeda_AB_To_XR (int nvar, CCTK_REAL A, CCTK_REAL B, CCTK_REAL *X,
	       CCTK_REAL *R, derivs U);
void GibbonsMaeda_C_To_c (int nvar, CCTK_REAL X, CCTK_REAL R, CCTK_REAL *x,
	     CCTK_REAL *r, derivs U);
void GibbonsMaeda_rx3_To_xyz (int nvar, CCTK_REAL x, CCTK_REAL r, CCTK_REAL phi, CCTK_REAL *y,
		 CCTK_REAL *z, derivs U);

/* Routines in  "Equations.c"*/
CCTK_REAL GibbonsMaeda_BY_KKofxyz (CCTK_REAL x, CCTK_REAL y, CCTK_REAL z);
void GibbonsMaeda_BY_Aijofxyz (CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL Aij[3][3]);
void GibbonsMaeda_Calc_E (CCTK_REAL x1, CCTK_REAL y1, CCTK_REAL z1, CCTK_REAL Ei[3], CCTK_REAL *E2);
void GibbonsMaeda_NonLinEquations (CCTK_REAL rho_adm,
          CCTK_REAL A, CCTK_REAL B, CCTK_REAL X, CCTK_REAL R,
		      CCTK_REAL x, CCTK_REAL r, CCTK_REAL phi,
		      CCTK_REAL y, CCTK_REAL z, derivs U, CCTK_REAL *values);
void GibbonsMaeda_LinEquations (CCTK_REAL A, CCTK_REAL B, CCTK_REAL X, CCTK_REAL R,
		   CCTK_REAL x, CCTK_REAL r, CCTK_REAL phi,
		   CCTK_REAL y, CCTK_REAL z, derivs dU, derivs U, CCTK_REAL *values);

/* Routines in  "Newton.c"*/
void GibbonsMaeda_TestRelax (CCTK_POINTER_TO_CONST cctkGH,
                int nvar, int n1, int n2, int n3, derivs v, CCTK_REAL *dv);
void GibbonsMaeda_Newton (CCTK_POINTER_TO_CONST cctkGH,
             int nvar, int n1, int n2, int n3, derivs v,
	           CCTK_REAL tol, int itmax);


/* 
 27: -1.325691774825335e-03
 37: -1.325691778944117e-03
 47: -1.325691778942711e-03
 
 17: -1.510625972641537e-03
 21: -1.511443006977708e-03
 27: -1.511440785153687e-03
 37: -1.511440809549005e-03
 39: -1.511440809597588e-03
 */
