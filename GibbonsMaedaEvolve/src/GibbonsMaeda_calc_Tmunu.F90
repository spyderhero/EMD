#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

subroutine GibbonsMaeda_calc_Tmunu( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  ! Fundamental variables
  CCTK_REAL                alph, beta(3)
  CCTK_REAL                gg(3,3), gu(3,3), detgg
  CCTK_REAL                lE(3), lB(3), lA(3), lAphi
  CCTK_REAL                Ed(3), Bd(3)
  CCTK_REAL                lphi1, lKphi1
  CCTK_REAL                lphi2, lKphi2
  CCTK_REAL                Tab(4,4)

  ! First derivatives
  CCTK_REAL                d1_lA(3,3)
  CCTK_REAL                d1_lphi1(3)
  CCTK_REAL                d1_lphi2(3)

  ! Auxiliary variables
  CCTK_REAL                eps_lc_d(3,3,3), eps_lc_u(3,3,3)

  ! Matter variables
  CCTK_REAL                srcE, srcjdi(3), srcSij(3,3)
  CCTK_REAL                srcE_p, srcjdi_p(3), srcSij_p(3,3)
  CCTK_REAL                srcE_F, srcjdi_F(3), srcSij_F(3,3)
  CCTK_REAL                jr_p, jr_F, Sir_p(3), Sir_F(3)

  ! Misc variables
  CCTK_REAL                dx12, dy12, dz12
  CCTK_REAL                odx60, ody60, odz60
  CCTK_REAL                aux_p, aux_F
  CCTK_REAL                xx(3), rr

  CCTK_REAL, parameter ::  one  = 1
  CCTK_REAL, parameter ::  pi   = acos(-one)
  CCTK_REAL, parameter ::  pi4  = 4*pi
  CCTK_REAL, parameter ::  pi8  = 8*pi
  CCTK_INT                 i, j, k
  CCTK_INT                 a, b, c, m
  
  ! jacobian
  integer                  istat
  logical                  use_jacobian
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: lJ11, lJ12, lJ13
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: lJ21, lJ22, lJ23
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: lJ31, lJ32, lJ33
  CCTK_POINTER             lJ11_ptr, lJ12_ptr, lJ13_ptr
  CCTK_POINTER             lJ21_ptr, lJ22_ptr, lJ23_ptr
  CCTK_POINTER             lJ31_ptr, lJ32_ptr, lJ33_ptr
  CCTK_REAL                jac(3,3)

  pointer (lJ11_ptr, lJ11), (lJ12_ptr, lJ12), (lJ13_ptr, lJ13)
  pointer (lJ21_ptr, lJ21), (lJ22_ptr, lJ22), (lJ23_ptr, lJ23)
  pointer (lJ31_ptr, lJ31), (lJ32_ptr, lJ32), (lJ33_ptr, lJ33)

  call CCTK_IsFunctionAliased(istat, "MultiPatch_GetDomainSpecification")
  if (istat == 0) then
     use_jacobian = .false.
  else
     use_jacobian = .true.
  end if

  if (use_jacobian) then
     call CCTK_VarDataPtr(lJ11_ptr, cctkGH, 0, "Coordinates::J11")
     call CCTK_VarDataPtr(lJ12_ptr, cctkGH, 0, "Coordinates::J12")
     call CCTK_VarDataPtr(lJ13_ptr, cctkGH, 0, "Coordinates::J13")
     call CCTK_VarDataPtr(lJ21_ptr, cctkGH, 0, "Coordinates::J21")
     call CCTK_VarDataPtr(lJ22_ptr, cctkGH, 0, "Coordinates::J22")
     call CCTK_VarDataPtr(lJ23_ptr, cctkGH, 0, "Coordinates::J23")
     call CCTK_VarDataPtr(lJ31_ptr, cctkGH, 0, "Coordinates::J31")
     call CCTK_VarDataPtr(lJ32_ptr, cctkGH, 0, "Coordinates::J32")
     call CCTK_VarDataPtr(lJ33_ptr, cctkGH, 0, "Coordinates::J33")
  end if

  dx12 = 12*CCTK_DELTA_SPACE(1)
  dy12 = 12*CCTK_DELTA_SPACE(2)
  dz12 = 12*CCTK_DELTA_SPACE(3)

  odx60 = 1 / (60 * CCTK_DELTA_SPACE(1))
  ody60 = 1 / (60 * CCTK_DELTA_SPACE(2))
  odz60 = 1 / (60 * CCTK_DELTA_SPACE(3))

  !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(k,j,i,a,b,c,m,aux_p,aux_F,&
  !$OMP                                 alph,beta,&
  !$OMP                                 gg,gu,detgg,&
  !$OMP                                 lE,lB,lA,lAphi,&
  !$OMP                                 Ed,Bd,Tab,&
  !$OMP                                 lphi1, lKphi1,&
  !$OMP                                 lphi2, lKphi2,&
  !$OMP                                 d1_lphi1,d1_lphi2, &
  !$OMP                                 eps_lc_d,eps_lc_u,&
  !$OMP                                 d1_lA, &
  !$OMP                                 srcE, srcjdi, srcSij,&
  !$OMP                                 srcE_p, srcjdi_p, srcSij_p,&
  !$OMP                                 srcE_F, srcjdi_F, srcSij_F,&
  !$OMP                                 xx, rr,&
  !$OMP                                 jr_p, jr_F, Sir_p, Sir_F, jac)
  do k = 1+cctk_nghostzones(3), cctk_lsh(3)-cctk_nghostzones(3)
     do j = 1+cctk_nghostzones(2), cctk_lsh(2)-cctk_nghostzones(2)
        do i = 1+cctk_nghostzones(1), cctk_lsh(1)-cctk_nghostzones(1)

           !------------ Get local variables ----------

           alph      = alp(i,j,k)

           beta(1)   = betax(i,j,k)
           beta(2)   = betay(i,j,k)
           beta(3)   = betaz(i,j,k)

           gg(1,1)   = gxx(i,j,k)
           gg(1,2)   = gxy(i,j,k)
           gg(1,3)   = gxz(i,j,k)
           gg(2,2)   = gyy(i,j,k)
           gg(2,3)   = gyz(i,j,k)
           gg(3,3)   = gzz(i,j,k)
           gg(2,1)   = gg(1,2)
           gg(3,1)   = gg(1,3)
           gg(3,2)   = gg(2,3)
           
           if (use_jacobian) then
              jac(1,1) = lJ11(i,j,k)
              jac(1,2) = lJ12(i,j,k)
              jac(1,3) = lJ13(i,j,k)
              jac(2,1) = lJ21(i,j,k)
              jac(2,2) = lJ22(i,j,k)
              jac(2,3) = lJ23(i,j,k)
              jac(3,1) = lJ31(i,j,k)
              jac(3,2) = lJ32(i,j,k)
              jac(3,3) = lJ33(i,j,k)
           else
              jac      = 0.0d0
              jac(1,1) = 1.0d0
              jac(2,2) = 1.0d0
              jac(3,3) = 1.0d0
           end if

           lE(1)     = Ex(i,j,k)
           lE(2)     = Ey(i,j,k)
           lE(3)     = Ez(i,j,k)

           lA(1)     = Ax(i,j,k)
           lA(2)     = Ay(i,j,k)
           lA(3)     = Az(i,j,k)

           lAphi     = Aphi(i,j,k)

           Ed(1)     = gg(1,1) * lE(1) + gg(1,2) * lE(2) + gg(1,3) * lE(3)
           Ed(2)     = gg(2,1) * lE(1) + gg(2,2) * lE(2) + gg(2,3) * lE(3)
           Ed(3)     = gg(3,1) * lE(1) + gg(3,2) * lE(2) + gg(3,3) * lE(3)
           
           lphi1     = phi1(i,j,k)
           lphi2     = phi2(i,j,k)

           lKphi1    = Kphi1(i,j,k)
           lKphi2    = Kphi2(i,j,k)


           !------------ Invert 3-metric ----------------
           detgg   =     gg(1,1) * gg(2,2) * gg(3,3)                              &
                   + 2 * gg(1,2) * gg(1,3) * gg(2,3)                              &
                   -     gg(1,1) * gg(2,3) ** 2                                   &
                   -     gg(2,2) * gg(1,3) ** 2                                   &
                   -     gg(3,3) * gg(1,2) ** 2

           gu(1,1) = (gg(2,2) * gg(3,3) - gg(2,3) ** 2     ) / detgg
           gu(2,2) = (gg(1,1) * gg(3,3) - gg(1,3) ** 2     ) / detgg
           gu(3,3) = (gg(1,1) * gg(2,2) - gg(1,2) ** 2     ) / detgg
           gu(1,2) = (gg(1,3) * gg(2,3) - gg(1,2) * gg(3,3)) / detgg
           gu(1,3) = (gg(1,2) * gg(2,3) - gg(1,3) * gg(2,2)) / detgg
           gu(2,3) = (gg(1,3) * gg(1,2) - gg(2,3) * gg(1,1)) / detgg
           gu(2,1) = gu(1,2)
           gu(3,1) = gu(1,3)
           gu(3,2) = gu(2,3)
           !-------------------------------------------------


           !------------- Centered 1st derivatives ----------

           if (derivs_order == 4) then

             ! d1_lA(3,3)
             d1_lA(1,1) = (   -Ax(i+2,j,k) + 8*Ax(i+1,j,k)               &
                           - 8*Ax(i-1,j,k) +   Ax(i-2,j,k) ) / dx12
             d1_lA(1,2) = (   -Ax(i,j+2,k) + 8*Ax(i,j+1,k)               &
                           - 8*Ax(i,j-1,k) +   Ax(i,j-2,k) ) / dy12
             d1_lA(1,3) = (   -Ax(i,j,k+2) + 8*Ax(i,j,k+1)               &
                           - 8*Ax(i,j,k-1) +   Ax(i,j,k-2) ) / dz12

             d1_lA(2,1) = (   -Ay(i+2,j,k) + 8*Ay(i+1,j,k)               &
                           - 8*Ay(i-1,j,k) +   Ay(i-2,j,k) ) / dx12
             d1_lA(2,2) = (   -Ay(i,j+2,k) + 8*Ay(i,j+1,k)               &
                           - 8*Ay(i,j-1,k) +   Ay(i,j-2,k) ) / dy12
             d1_lA(2,3) = (   -Ay(i,j,k+2) + 8*Ay(i,j,k+1)               &
                           - 8*Ay(i,j,k-1) +   Ay(i,j,k-2) ) / dz12

             d1_lA(3,1) = (   -Az(i+2,j,k) + 8*Az(i+1,j,k)               &
                           - 8*Az(i-1,j,k) +   Az(i-2,j,k) ) / dx12
             d1_lA(3,2) = (   -Az(i,j+2,k) + 8*Az(i,j+1,k)               &
                           - 8*Az(i,j-1,k) +   Az(i,j-2,k) ) / dy12
             d1_lA(3,3) = (   -Az(i,j,k+2) + 8*Az(i,j,k+1)               &
                           - 8*Az(i,j,k-1) +   Az(i,j,k-2) ) / dz12
                           
             ! d1_lphi1(3)
             d1_lphi1(1) = (   -phi1(i+2,j,k) + 8*phi1(i+1,j,k)           &
                            - 8*phi1(i-1,j,k) +   phi1(i-2,j,k) ) / dx12

             d1_lphi1(2) = (   -phi1(i,j+2,k) + 8*phi1(i,j+1,k)           &
                            - 8*phi1(i,j-1,k) +   phi1(i,j-2,k) ) / dy12

             d1_lphi1(3) = (   -phi1(i,j,k+2) + 8*phi1(i,j,k+1)           &
                            - 8*phi1(i,j,k-1) +   phi1(i,j,k-2) ) / dz12

             ! d1_lphi2(3)
             d1_lphi2(1) = (   -phi2(i+2,j,k) + 8*phi2(i+1,j,k)           &
                            - 8*phi2(i-1,j,k) +   phi2(i-2,j,k) ) / dx12

             d1_lphi2(2) = (   -phi2(i,j+2,k) + 8*phi2(i,j+1,k)           &
                            - 8*phi2(i,j-1,k) +   phi2(i,j-2,k) ) / dy12

             d1_lphi2(3) = (   -phi2(i,j,k+2) + 8*phi2(i,j,k+1)           &
                            - 8*phi2(i,j,k-1) +   phi2(i,j,k-2) ) / dz12

           else if (derivs_order == 6) then

             ! d1_lA(3,3)
             d1_lA(1,1) = (  Ax(i+3,j,k) - 9*Ax(i+2,j,k) + 45*Ax(i+1,j,k) &
                           - Ax(i-3,j,k) + 9*Ax(i-2,j,k) - 45*Ax(i-1,j,k) ) * odx60
             d1_lA(1,2) = (  Ax(i,j+3,k) - 9*Ax(i,j+2,k) + 45*Ax(i,j+1,k) &
                           - Ax(i,j-3,k) + 9*Ax(i,j-2,k) - 45*Ax(i,j-1,k) ) * ody60
             d1_lA(1,3) = (  Ax(i,j,k+3) - 9*Ax(i,j,k+2) + 45*Ax(i,j,k+1) &
                           - Ax(i,j,k-3) + 9*Ax(i,j,k-2) - 45*Ax(i,j,k-1) ) * odz60

             d1_lA(2,1) = (  Ay(i+3,j,k) - 9*Ay(i+2,j,k) + 45*Ay(i+1,j,k) &
                           - Ay(i-3,j,k) + 9*Ay(i-2,j,k) - 45*Ay(i-1,j,k) ) * odx60
             d1_lA(2,2) = (  Ay(i,j+3,k) - 9*Ay(i,j+2,k) + 45*Ay(i,j+1,k) &
                           - Ay(i,j-3,k) + 9*Ay(i,j-2,k) - 45*Ay(i,j-1,k) ) * ody60
             d1_lA(2,3) = (  Ay(i,j,k+3) - 9*Ay(i,j,k+2) + 45*Ay(i,j,k+1) &
                           - Ay(i,j,k-3) + 9*Ay(i,j,k-2) - 45*Ay(i,j,k-1) ) * odz60

             d1_lA(3,1) = (  Az(i+3,j,k) - 9*Az(i+2,j,k) + 45*Az(i+1,j,k) &
                           - Az(i-3,j,k) + 9*Az(i-2,j,k) - 45*Az(i-1,j,k) ) * odx60
             d1_lA(3,2) = (  Az(i,j+3,k) - 9*Az(i,j+2,k) + 45*Az(i,j+1,k) &
                           - Az(i,j-3,k) + 9*Az(i,j-2,k) - 45*Az(i,j-1,k) ) * ody60
             d1_lA(3,3) = (  Az(i,j,k+3) - 9*Az(i,j,k+2) + 45*Az(i,j,k+1) &
                           - Az(i,j,k-3) + 9*Az(i,j,k-2) - 45*Az(i,j,k-1) ) * odz60
                           
             d1_lphi1(1) = (  phi1(i+3,j,k) - 9*phi1(i+2,j,k) + 45*phi1(i+1,j,k) &
                            - phi1(i-3,j,k) + 9*phi1(i-2,j,k) - 45*phi1(i-1,j,k) ) * odx60

             d1_lphi1(2) = (  phi1(i,j+3,k) - 9*phi1(i,j+2,k) + 45*phi1(i,j+1,k) &
                            - phi1(i,j-3,k) + 9*phi1(i,j-2,k) - 45*phi1(i,j-1,k) ) * ody60

             d1_lphi1(3) = (  phi1(i,j,k+3) - 9*phi1(i,j,k+2) + 45*phi1(i,j,k+1) &
                            - phi1(i,j,k-3) + 9*phi1(i,j,k-2) - 45*phi1(i,j,k-1) ) * odz60


             d1_lphi2(1) = (  phi2(i+3,j,k) - 9*phi2(i+2,j,k) + 45*phi2(i+1,j,k) &
                            - phi2(i-3,j,k) + 9*phi2(i-2,j,k) - 45*phi2(i-1,j,k) ) * odx60

             d1_lphi2(2) = (  phi2(i,j+3,k) - 9*phi2(i,j+2,k) + 45*phi2(i,j+1,k) &
                            - phi2(i,j-3,k) + 9*phi2(i,j-2,k) - 45*phi2(i,j-1,k) ) * ody60

             d1_lphi2(3) = (  phi2(i,j,k+3) - 9*phi2(i,j,k+2) + 45*phi2(i,j,k+1) &
                            - phi2(i,j,k-3) + 9*phi2(i,j,k-2) - 45*phi2(i,j,k-1) ) * odz60

           else
             call CCTK_WARN(0, "derivs_order not yet implemented.")
           end if
           
           if (use_jacobian) then
              call GibbonsMaeda_d1_Vector_apply_jacobian(d1_lA, jac)
              call GibbonsMaeda_d1_Scalar_apply_jacobian(d1_lphi1, jac)
              call GibbonsMaeda_d1_Scalar_apply_jacobian(d1_lphi2, jac)
           end if

           !------------ Levi-Civita tensor ----------
           eps_lc_u        =  0
           eps_lc_u(1,2,3) =  1
           eps_lc_u(2,3,1) =  1
           eps_lc_u(3,1,2) =  1
           eps_lc_u(3,2,1) = -1
           eps_lc_u(2,1,3) = -1
           eps_lc_u(1,3,2) = -1
           eps_lc_u = eps_lc_u / sqrt(detgg)
           eps_lc_d = eps_lc_u * detgg
           !------------------------------------------


           ! magnetic field B (here used as an auxiliary variable)
           lB = 0
           do a = 1, 3
             do b = 1, 3
               do m = 1, 3
                 lB(a) = lB(a) + eps_lc_u(a,b,m) * d1_lA(m,b)
               end do
             end do
           end do

           Bd(1) = gg(1,1) * lB(1) + gg(1,2) * lB(2) + gg(1,3) * lB(3)
           Bd(2) = gg(2,1) * lB(1) + gg(2,2) * lB(2) + gg(2,3) * lB(3)
           Bd(3) = gg(3,1) * lB(1) + gg(3,2) * lB(2) + gg(3,3) * lB(3)
           !-------------------------------------------

           !------------ Matter terms -----------------
           !
           ! mu = 0, 1, 2, 3; i,a = 1,2,3
           !
           ! n_mu = (-alph, 0, 0, 0)
           ! n^mu = (1, -betax, -betay, -betaz)/alph
           !
           ! rho = n^mu n^nu T_{mu nu}
           !     = (T_{00} - 2 beta^i T_{i0} + beta^i beta^j T_{ij})/alph^2
           !
           ! j_a = -h_a^mu n^nu T_{mu nu}
           !     = -(T_{a 0} - beta^j T_{a j})/alph
           !
           ! S_{a b} = h_{a mu} h_{b nu} T^{mu nu} = T_{a b}


           ! srcE_p = rho = n^mu n^nu T_{mu nu}^(phi)
           ! srcE_F = rho = n^mu n^nu T_{mu nu}^F
           ! srcE = rho = n^mu n^nu T_{mu nu} = srcE_p + srcE_F * Exp(-2a phi)
           ! a is coupling constant
           srcE_p = lKphi1*lKphi1
           srcE_F = 0
           do a = 1, 3
              do b = 1, 3
                 srcE_p = srcE_p + ( d1_lphi1(a) * d1_lphi1(b) ) * gu(a,b)
                 srcE_F = srcE_F + ( lE(a) * lE(b) + lB(a) * lB(b) ) * gg(a,b)
              end do
           end do
           srcE = ( srcE_p + EXP(-2 * coupling_constant * lphi1) * srcE_F ) / pi8

           !srcjdi = j_a
           srcjdi_p = 0
           srcjdi_F = 0
           do a = 1, 3
              srcjdi_p(a) = srcjdi_p(a) + d1_lphi1(a) * lKphi1
              do b = 1, 3
                  srcjdi_F(a) = srcjdi_F(a) + lE(b) * ( d1_lA(b,a) - d1_lA(a,b) )
              end do
           end do
           srcjdi = ( srcjdi_p + EXP(-2 * coupling_constant * lphi1) * srcjdi_F ) / pi4


           ! srcSij = S_{a b}

           aux_p = lKphi1 * lKphi1
           aux_F = 0
           do a = 1, 3
              do b = 1, 3
                 aux_p = aux_p - d1_lphi1(a) * d1_lphi1(b) * gu(a,b)
                 aux_F = aux_F + ( lE(a) * lE(b) + lB(a) * lB(b) ) * gg(a,b)
              end do
           end do

           srcSij_p = 0.5 * aux_p * gg
           srcSij_F = 0.5 * aux_F * gg
           do a = 1, 3
              do b = 1, 3
                 srcSij_p(a,b) = srcSij_p(a,b) + d1_lphi1(a) * d1_lphi1(b)
                 srcSij_F(a,b) = srcSij_F(a,b) - ( Ed(a) * Ed(b) + Bd(a) * Bd(b) )
              end do
           end do
           srcSij = ( srcSij_p + EXP(-2 * coupling_constant * lphi1) * srcSij_F ) / pi4
           
           !------------------------------------------
           ! transform to spherical coordinates
           
           xx(1)     = x(i,j,k)
    	   xx(2)     = y(i,j,k)
    	   xx(3)     = z(i,j,k)

           rr = sqrt( xx(1)**2 + xx(2)**2 + xx(3)**2 )
    	   if( rr < eps_r ) rr = eps_r
    	   
    	   jr_p = 0
    	   jr_F = 0
           do a = 1, 3
              jr_p = jr_p + srcjdi_p(a) * xx(a) / rr
              jr_F = jr_F + srcjdi_F(a) * xx(a) / rr
           end do
           
           Sir_p = 0
           Sir_F = 0
           do a = 1, 3
              do b = 1,3
                 Sir_p(a) = Sir_p(a) + srcSij_p(a,b) * xx(b) / rr
                 Sir_F(a) = Sir_F(a) + srcSij_F(a,b) * xx(b) / rr
              end do
           end do
           
           !store it in the Spherical coordinates variables
           jrPhi_gf(i,j,k)  = jr_p
           jrF_gf(i,j,k)    = jr_F
           
           SxrPhi_gf(i,j,k) = Sir_p(1)
           SyrPhi_gf(i,j,k) = Sir_p(2)
           SzrPhi_gf(i,j,k) = Sir_p(3)
           
           SxrF_gf(i,j,k)   = Sir_F(1)
           SyrF_gf(i,j,k)   = Sir_F(2)
           SzrF_gf(i,j,k)   = Sir_F(3)
    	   

           !------------------------------------------


           ! now to fill in the stress-energy tensor. note that we use Tab(4,4)
           ! for T_{0 0}
           !
           ! T_{a b} = S_{a b}
           !
           ! T_{0 0} = alph^2 rho - 2 alph beta^a j_a + beta^a beta^b S_{a b}
           !
           ! T_{0 a} = -alph j_a + beta^b S_{a b}

           Tab(1:3,1:3) = srcSij(1:3,1:3)

           Tab(1:3,4) = -alph * srcjdi(1:3)
           do b = 1, 3
              Tab(1:3,4) = Tab(1:3,4) + beta(b) * srcSij(1:3,b)
           end do
           Tab(4,1:3) = Tab(1:3,4)

           Tab(4,4) = alph**2 * srcE
           do a = 1, 3
              Tab(4,4) = Tab(4,4) - 2 * alph * beta(a) * srcjdi(a)
              do b = 1, 3
                 Tab(4,4) = Tab(4,4) + beta(a) * beta(b) * srcSij(a,b)
              end do
           end do

           ! and finally store it in the Tmunu variables
           eTtt(i,j,k) = eTtt(i,j,k) + Tab(4,4)
           eTtx(i,j,k) = eTtx(i,j,k) + Tab(4,1)
           eTty(i,j,k) = eTty(i,j,k) + Tab(4,2)
           eTtz(i,j,k) = eTtz(i,j,k) + Tab(4,3)
           eTxx(i,j,k) = eTxx(i,j,k) + Tab(1,1)
           eTxy(i,j,k) = eTxy(i,j,k) + Tab(1,2)
           eTxz(i,j,k) = eTxz(i,j,k) + Tab(1,3)
           eTyy(i,j,k) = eTyy(i,j,k) + Tab(2,2)
           eTyz(i,j,k) = eTyz(i,j,k) + Tab(2,3)
           eTzz(i,j,k) = eTzz(i,j,k) + Tab(3,3)

        end do
     end do
  end do
  !$OMP END PARALLEL DO

end subroutine GibbonsMaeda_calc_Tmunu
