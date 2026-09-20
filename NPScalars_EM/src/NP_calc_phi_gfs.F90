! NPEM
! NP_calc_phi_gfs.F90 : Actual calculation of the NP scalar grid functions
!
!=============================================================================

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine NPEM_calcPhiGF( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL gd(3,3), detgd
  CCTK_REAL Ed(3), Bd(3), Bu(3)
  CCTK_REAL d1_lA(3,3)
  CCTK_REAL dx(3), xx(3)
  CCTK_REAL dx12, dy12, dz12, dxsq12, dysq12, dzsq12, dxdy144, dxdz144,    &
            dydz144
  CCTK_REAL odx60, ody60, odz60, odxsq180, odysq180, odzsq180,             &
            odxdy3600, odxdz3600, odydz3600
  CCTK_REAL dx2, dy2, dz2, dxsq, dysq, dzsq, dxdy4, dxdz4, dydz4
  CCTK_REAL u_vec(3), v_vec(3), w_vec(3), dotp1, dotp2
  CCTK_REAL eps_lc_u(3,3,3)
  CCTK_REAL xdvar(3,3)

  CCTK_INT  i, j, k, m, n, p, a, b, c, d
  
  ! jacobian
  integer                  istat
  logical                  use_jacobian
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: lJ11, lJ12, lJ13
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: lJ21, lJ22, lJ23
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: lJ31, lJ32, lJ33
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: ldJ111, ldJ112, ldJ113, ldJ122, ldJ123, ldJ133
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: ldJ211, ldJ212, ldJ213, ldJ222, ldJ223, ldJ233
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: ldJ311, ldJ312, ldJ313, ldJ322, ldJ323, ldJ333

  CCTK_POINTER             lJ11_ptr, lJ12_ptr, lJ13_ptr
  CCTK_POINTER             lJ21_ptr, lJ22_ptr, lJ23_ptr
  CCTK_POINTER             lJ31_ptr, lJ32_ptr, lJ33_ptr
  CCTK_POINTER             ldJ111_ptr, ldJ112_ptr, ldJ113_ptr, ldJ122_ptr, ldJ123_ptr, ldJ133_ptr
  CCTK_POINTER             ldJ211_ptr, ldJ212_ptr, ldJ213_ptr, ldJ222_ptr, ldJ223_ptr, ldJ233_ptr
  CCTK_POINTER             ldJ311_ptr, ldJ312_ptr, ldJ313_ptr, ldJ322_ptr, ldJ323_ptr, ldJ333_ptr

  CCTK_REAL                jac(3,3)

  pointer (lJ11_ptr, lJ11), (lJ12_ptr, lJ12), (lJ13_ptr, lJ13)
  pointer (lJ21_ptr, lJ21), (lJ22_ptr, lJ22), (lJ23_ptr, lJ23)
  pointer (lJ31_ptr, lJ31), (lJ32_ptr, lJ32), (lJ33_ptr, lJ33)

  pointer (ldJ111_ptr, ldJ111), (ldJ112_ptr, ldJ112), (ldJ113_ptr, ldJ113), (ldJ122_ptr, ldJ122), (ldJ123_ptr, ldJ123), (ldJ133_ptr, ldJ133)
  pointer (ldJ211_ptr, ldJ211), (ldJ212_ptr, ldJ212), (ldJ213_ptr, ldJ213), (ldJ222_ptr, ldJ222), (ldJ223_ptr, ldJ223), (ldJ233_ptr, ldJ233)
  pointer (ldJ311_ptr, ldJ311), (ldJ312_ptr, ldJ312), (ldJ313_ptr, ldJ313), (ldJ322_ptr, ldJ322), (ldJ323_ptr, ldJ323), (ldJ333_ptr, ldJ333)

  ! TODO: can this be active but with a cartesian mapping choice?
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

     call CCTK_VarDataPtr(ldJ111_ptr, cctkGH, 0, "Coordinates::dJ111")
     call CCTK_VarDataPtr(ldJ112_ptr, cctkGH, 0, "Coordinates::dJ112")
     call CCTK_VarDataPtr(ldJ113_ptr, cctkGH, 0, "Coordinates::dJ113")
     call CCTK_VarDataPtr(ldJ122_ptr, cctkGH, 0, "Coordinates::dJ122")
     call CCTK_VarDataPtr(ldJ123_ptr, cctkGH, 0, "Coordinates::dJ123")
     call CCTK_VarDataPtr(ldJ133_ptr, cctkGH, 0, "Coordinates::dJ133")

     call CCTK_VarDataPtr(ldJ211_ptr, cctkGH, 0, "Coordinates::dJ211")
     call CCTK_VarDataPtr(ldJ212_ptr, cctkGH, 0, "Coordinates::dJ212")
     call CCTK_VarDataPtr(ldJ213_ptr, cctkGH, 0, "Coordinates::dJ213")
     call CCTK_VarDataPtr(ldJ222_ptr, cctkGH, 0, "Coordinates::dJ222")
     call CCTK_VarDataPtr(ldJ223_ptr, cctkGH, 0, "Coordinates::dJ223")
     call CCTK_VarDataPtr(ldJ233_ptr, cctkGH, 0, "Coordinates::dJ233")

     call CCTK_VarDataPtr(ldJ311_ptr, cctkGH, 0, "Coordinates::dJ311")
     call CCTK_VarDataPtr(ldJ312_ptr, cctkGH, 0, "Coordinates::dJ312")
     call CCTK_VarDataPtr(ldJ313_ptr, cctkGH, 0, "Coordinates::dJ313")
     call CCTK_VarDataPtr(ldJ322_ptr, cctkGH, 0, "Coordinates::dJ322")
     call CCTK_VarDataPtr(ldJ323_ptr, cctkGH, 0, "Coordinates::dJ323")
     call CCTK_VarDataPtr(ldJ333_ptr, cctkGH, 0, "Coordinates::dJ333")
  end if


  dx(:)   = CCTK_DELTA_SPACE(:)

  !--- coefficients for sixth order
  odx60 = 1 / (60 * CCTK_DELTA_SPACE(1))
  ody60 = 1 / (60 * CCTK_DELTA_SPACE(2))
  odz60 = 1 / (60 * CCTK_DELTA_SPACE(3))

  odxsq180 = 1 / (180*CCTK_DELTA_SPACE(1)**2)
  odysq180 = 1 / (180*CCTK_DELTA_SPACE(2)**2)
  odzsq180 = 1 / (180*CCTK_DELTA_SPACE(3)**2)

  odxdy3600 = 1 / (3600*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(2))
  odxdz3600 = 1 / (3600*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(3))
  odydz3600 = 1 / (3600*CCTK_DELTA_SPACE(2)*CCTK_DELTA_SPACE(3))

  !--- coefficients for fourth order
  dx12    = 12 * dx(1)
  dy12    = 12 * dx(2)
  dz12    = 12 * dx(3)
  dxsq12  = 12 * dx(1)**2
  dysq12  = 12 * dx(2)**2
  dzsq12  = 12 * dx(3)**2
  dxdy144 = dx12 * dy12
  dxdz144 = dx12 * dz12
  dydz144 = dy12 * dz12

  !--- coefficients for second order
  dx2   = 2 * dx(1)
  dy2   = 2 * dx(2)
  dz2   = 2 * dx(3)
  dxsq  = dx(1)**2
  dysq  = dx(2)**2
  dzsq  = dx(3)**2
  dxdy4 = 4*dx(1) * dx(2)
  dxdz4 = 4*dx(1) * dx(3)
  dydz4 = 4*dx(2) * dx(3)

  !=== initialize source terms (gfs) as zero ===
  if (calc_phi /= 0) then
     phi0re = 0.0
     phi0im = 0.0
     phi1re = 0.0
     phi1im = 0.0
     phi2re = 0.0
     phi2im = 0.0
  end if

  do k = 1+cctk_nghostzones(3), cctk_lsh(3)-cctk_nghostzones(3)
  do j = 1+cctk_nghostzones(2), cctk_lsh(2)-cctk_nghostzones(2)
  do i = 1+cctk_nghostzones(1), cctk_lsh(1)-cctk_nghostzones(1)
  
  !-------------------------------------------
    
    ! hes(i,j,k) = jac(i,j),k
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
       jac      = 0.0
       jac(1,1) = 1.0
       jac(2,2) = 1.0
       jac(3,3) = 1.0
    end if
    ! write(*,*) 'J = ', jac

    !-------------------------------------------

    !--------------Get local variables ----------
    gd(1,1) = gxx(i,j,k)
    gd(1,2) = gxy(i,j,k)
    gd(1,3) = gxz(i,j,k)
    gd(2,2) = gyy(i,j,k)
    gd(2,3) = gyz(i,j,k)
    gd(3,3) = gzz(i,j,k)
    gd(2,1) = gd(1,2)
    gd(3,1) = gd(1,3)
    gd(3,2) = gd(2,3)

    ! Proca vars
    Ed(1)   = gd(1,1) * Ex(i,j,k) + gd(1,2) * Ey(i,j,k) + gd(1,3) * Ez(i,j,k)
    Ed(2)   = gd(2,1) * Ex(i,j,k) + gd(2,2) * Ey(i,j,k) + gd(2,3) * Ez(i,j,k)
    Ed(3)   = gd(3,1) * Ex(i,j,k) + gd(3,2) * Ey(i,j,k) + gd(3,3) * Ez(i,j,k)

    !--------------------------------------------


    !-------------- det metric ---------------
    detgd =       gd(1,1) * gd(2,2) * gd(3,3)                                &
            + 2 * gd(1,2) * gd(1,3) * gd(2,3)                                &
            -     gd(1,1) * gd(2,3) ** 2                                     &
            -     gd(2,2) * gd(1,3) ** 2                                     &
            -     gd(3,3) * gd(1,2) ** 2
    !--------------------------------------------
    if ( NP_order == 6 ) then
      ! Sixth order derivatives

      !-------------- Centered 1st derivatives ----
      ! d1_lA(3,3)
      d1_lA(1,1) = (  Ax(i+3,j,k) - 9*Ax(i+2,j,k) + 45*Ax(i+1,j,k) &
                    - Ax(i-3,j,k) + 9*Ax(i-2,j,k) - 45*Ax(i-1,j,k) ) * odx60
      d1_lA(2,1) = (  Ay(i+3,j,k) - 9*Ay(i+2,j,k) + 45*Ay(i+1,j,k) &
                    - Ay(i-3,j,k) + 9*Ay(i-2,j,k) - 45*Ay(i-1,j,k) ) * odx60
      d1_lA(3,1) = (  Az(i+3,j,k) - 9*Az(i+2,j,k) + 45*Az(i+1,j,k) &
                    - Az(i-3,j,k) + 9*Az(i-2,j,k) - 45*Az(i-1,j,k) ) * odx60

      d1_lA(1,2) = (  Ax(i,j+3,k) - 9*Ax(i,j+2,k) + 45*Ax(i,j+1,k) &
                    - Ax(i,j-3,k) + 9*Ax(i,j-2,k) - 45*Ax(i,j-1,k) ) * ody60
      d1_lA(2,2) = (  Ay(i,j+3,k) - 9*Ay(i,j+2,k) + 45*Ay(i,j+1,k) &
                    - Ay(i,j-3,k) + 9*Ay(i,j-2,k) - 45*Ay(i,j-1,k) ) * ody60
      d1_lA(3,2) = (  Az(i,j+3,k) - 9*Az(i,j+2,k) + 45*Az(i,j+1,k) &
                    - Az(i,j-3,k) + 9*Az(i,j-2,k) - 45*Az(i,j-1,k) ) * ody60

      d1_lA(1,3) = (  Ax(i,j,k+3) - 9*Ax(i,j,k+2) + 45*Ax(i,j,k+1) &
                    - Ax(i,j,k-3) + 9*Ax(i,j,k-2) - 45*Ax(i,j,k-1) ) * odz60
      d1_lA(2,3) = (  Ay(i,j,k+3) - 9*Ay(i,j,k+2) + 45*Ay(i,j,k+1) &
                    - Ay(i,j,k-3) + 9*Ay(i,j,k-2) - 45*Ay(i,j,k-1) ) * odz60
      d1_lA(3,3) = (  Az(i,j,k+3) - 9*Az(i,j,k+2) + 45*Az(i,j,k+1) &
                    - Az(i,j,k-3) + 9*Az(i,j,k-2) - 45*Az(i,j,k-1) ) * odz60

      !--------------------------------------------

    else if ( NP_order == 4 ) then

      ! Fourth order derivatives

      !-------------- Centered 1st derivatives ----
      ! d1_lA(3,3)
      d1_lA(1,1) = (   -Ax(i+2,j,k) + 8*Ax(i+1,j,k)               &
                    - 8*Ax(i-1,j,k) +   Ax(i-2,j,k) ) / dx12
      d1_lA(2,1) = (   -Ay(i+2,j,k) + 8*Ay(i+1,j,k)               &
                    - 8*Ay(i-1,j,k) +   Ay(i-2,j,k) ) / dx12
      d1_lA(3,1) = (   -Az(i+2,j,k) + 8*Az(i+1,j,k)               &
                    - 8*Az(i-1,j,k) +   Az(i-2,j,k) ) / dx12

      d1_lA(1,2) = (   -Ax(i,j+2,k) + 8*Ax(i,j+1,k)               &
                    - 8*Ax(i,j-1,k) +   Ax(i,j-2,k) ) / dy12
      d1_lA(2,2) = (   -Ay(i,j+2,k) + 8*Ay(i,j+1,k)               &
                    - 8*Ay(i,j-1,k) +   Ay(i,j-2,k) ) / dy12
      d1_lA(3,2) = (   -Az(i,j+2,k) + 8*Az(i,j+1,k)               &
                    - 8*Az(i,j-1,k) +   Az(i,j-2,k) ) / dy12

      d1_lA(1,3) = (   -Ax(i,j,k+2) + 8*Ax(i,j,k+1)               &
                    - 8*Ax(i,j,k-1) +   Ax(i,j,k-2) ) / dz12
      d1_lA(2,3) = (   -Ay(i,j,k+2) + 8*Ay(i,j,k+1)               &
                    - 8*Ay(i,j,k-1) +   Ay(i,j,k-2) ) / dz12
      d1_lA(3,3) = (   -Az(i,j,k+2) + 8*Az(i,j,k+1)               &
                    - 8*Az(i,j,k-1) +   Az(i,j,k-2) ) / dz12

      !--------------------------------------------

    else

      ! second order derivatives as default

      !-------------- Centered 1st derivatives ----
      
      ! d1_lA(3,3)
      d1_lA(1,1) = (Ax(i+1,j,k) - Ax(i-1,j,k)) / dx2
      d1_lA(2,1) = (Ay(i+1,j,k) - Ay(i-1,j,k)) / dx2
      d1_lA(3,1) = (Az(i+1,j,k) - Az(i-1,j,k)) / dx2

      d1_lA(1,2) = (Ax(i,j+1,k) - Ax(i,j-1,k)) / dy2
      d1_lA(2,2) = (Ay(i,j+1,k) - Ay(i,j-1,k)) / dy2
      d1_lA(3,2) = (Az(i,j+1,k) - Az(i,j-1,k)) / dy2

      d1_lA(1,3) = (Ax(i,j,k+1) - Ax(i,j,k-1)) / dz2
      d1_lA(2,3) = (Ay(i,j,k+1) - Ay(i,j,k-1)) / dz2
      d1_lA(3,3) = (Az(i,j,k+1) - Az(i,j,k-1)) / dz2
      
    end if

    !------------------------------------------
    
    if (use_jacobian) then
      xdvar = 0

      do a = 1, 3
         do b = 1, 3
       	    do c = 1, 3
               xdvar(a,b) = xdvar(a,b) + d1_lA(a,c) * jac(c,b)
            end do
         end do
      end do

      d1_lA  = xdvar
      
    end if

    !------------ Levi-Civita tensor ----------
    eps_lc_u        = 0
    eps_lc_u(1,2,3) = 1
    eps_lc_u(2,3,1) = 1
    eps_lc_u(3,1,2) = 1
    eps_lc_u(3,2,1) = -1
    eps_lc_u(2,1,3) = -1
    eps_lc_u(1,3,2) = -1
    eps_lc_u = eps_lc_u / sqrt(detgd)

    ! eps_lc_d = eps_lc_u * detgd
    !------------------------------------------


    ! magnetic field
    Bu = 0
    do a = 1, 3
      do b = 1, 3
        do m = 1, 3
          Bu(a) = Bu(a) + eps_lc_u(a,b,m) * d1_lA(m,b)
        end do
      end do
    end do

    Bd(1)   = gd(1,1) * Bu(1) + gd(1,2) * Bu(2) + gd(1,3) * Bu(3)
    Bd(2)   = gd(2,1) * Bu(1) + gd(2,2) * Bu(2) + gd(2,3) * Bu(3)
    Bd(3)   = gd(3,1) * Bu(1) + gd(3,2) * Bu(2) + gd(3,3) * Bu(3)


    !------------ Orthonormal basis -----------
    ! Starting vectors
    xx(:) = (/ x(i,j,k), y(i,j,k), z(i,j,k) /)

    ! All points on the z-axis are pathological, since the triad vectors
    ! in the phi and theta direction are not well-defined. Take points
    ! just a little off, say at x = +epsilon.
    !
    ! we always assume z_orientation, in agreement with the Multipole thorn

    if( xx(1)**2 + xx(2)**2 < 1.0d-12 ) xx(1) = xx(1) + 1.0d-10

    u_vec(:) = xx(:)
    v_vec(:) = (/ xx(1)*xx(3), xx(2)*xx(3), -xx(1)**2 - xx(2)**2 /)
    w_vec(:) = (/ -xx(2), xx(1), 0.0d0 /)

    ! Orthonormalization
    dotp1 =   gd(1,1) * u_vec(1) * u_vec(1) + gd(1,2) * u_vec(1) * u_vec(2)&
            + gd(1,3) * u_vec(1) * u_vec(3) + gd(2,1) * u_vec(2) * u_vec(1)&
            + gd(2,2) * u_vec(2) * u_vec(2) + gd(2,3) * u_vec(2) * u_vec(3)&
            + gd(3,1) * u_vec(3) * u_vec(1) + gd(3,2) * u_vec(3) * u_vec(2)&
            + gd(3,3) * u_vec(3) * u_vec(3)
    u_vec = u_vec / sqrt(dotp1)

    dotp1 =   gd(1,1) * u_vec(1) * v_vec(1) + gd(1,2) * u_vec(1) * v_vec(2)&
            + gd(1,3) * u_vec(1) * v_vec(3) + gd(2,1) * u_vec(2) * v_vec(1)&
            + gd(2,2) * u_vec(2) * v_vec(2) + gd(2,3) * u_vec(2) * v_vec(3)&
            + gd(3,1) * u_vec(3) * v_vec(1) + gd(3,2) * u_vec(3) * v_vec(2)&
            + gd(3,3) * u_vec(3) * v_vec(3)
    v_vec = v_vec - dotp1 * u_vec

    dotp1 =   gd(1,1) * v_vec(1) * v_vec(1) + gd(1,2) * v_vec(1) * v_vec(2)&
            + gd(1,3) * v_vec(1) * v_vec(3) + gd(2,1) * v_vec(2) * v_vec(1)&
            + gd(2,2) * v_vec(2) * v_vec(2) + gd(2,3) * v_vec(2) * v_vec(3)&
            + gd(3,1) * v_vec(3) * v_vec(1) + gd(3,2) * v_vec(3) * v_vec(2)&
            + gd(3,3) * v_vec(3) * v_vec(3)
    v_vec = v_vec / sqrt(dotp1)

    dotp1 =   gd(1,1) * u_vec(1) * w_vec(1) + gd(1,2) * u_vec(1) * w_vec(2)&
            + gd(1,3) * u_vec(1) * w_vec(3) + gd(2,1) * u_vec(2) * w_vec(1)&
            + gd(2,2) * u_vec(2) * w_vec(2) + gd(2,3) * u_vec(2) * w_vec(3)&
            + gd(3,1) * u_vec(3) * w_vec(1) + gd(3,2) * u_vec(3) * w_vec(2)&
            + gd(3,3) * u_vec(3) * w_vec(3)

    dotp2 =   gd(1,1) * v_vec(1) * w_vec(1) + gd(1,2) * v_vec(1) * w_vec(2)&
            + gd(1,3) * v_vec(1) * w_vec(3) + gd(2,1) * v_vec(2) * w_vec(1)&
            + gd(2,2) * v_vec(2) * w_vec(2) + gd(2,3) * v_vec(2) * w_vec(3)&
            + gd(3,1) * v_vec(3) * w_vec(1) + gd(3,2) * v_vec(3) * w_vec(2)&
            + gd(3,3) * v_vec(3) * w_vec(3)
    w_vec = w_vec - dotp1 * u_vec - dotp2 * v_vec

    dotp1 =   gd(1,1) * w_vec(1) * w_vec(1) + gd(1,2) * w_vec(1) * w_vec(2)&
            + gd(1,3) * w_vec(1) * w_vec(3) + gd(2,1) * w_vec(2) * w_vec(1)&
            + gd(2,2) * w_vec(2) * w_vec(2) + gd(2,3) * w_vec(2) * w_vec(3)&
            + gd(3,1) * w_vec(3) * w_vec(1) + gd(3,2) * w_vec(3) * w_vec(2)&
            + gd(3,3) * w_vec(3) * w_vec(3)
    w_vec = w_vec / sqrt(dotp1)

    ! ud_vec = matmul( gd, u_vec )
    !------------------------------------------

    !------------ Phi0, Phi2, Phi1 ------------
    if (calc_phi /= 0) then
       do m = 1, 3
         phi0re(i,j,k) = phi0re(i,j,k) - 0.5 * (   Ed(m) * v_vec(m)          &
                                                 - Bd(m) * w_vec(m) )
         phi0im(i,j,k) = phi0im(i,j,k) - 0.5 * (   Ed(m) * w_vec(m)          &
                                                 + Bd(m) * v_vec(m) )
       
         phi2re(i,j,k) = phi2re(i,j,k) + 0.5 * (   Ed(m) * v_vec(m)          &
                                                 + Bd(m) * w_vec(m) )
         phi2im(i,j,k) = phi2im(i,j,k) - 0.5 * (   Ed(m) * w_vec(m)          &
                                                 - Bd(m) * v_vec(m) )
       
         phi1re(i,j,k) = phi1re(i,j,k) + 0.5 * Ed(m) * u_vec(m)
         phi1im(i,j,k) = phi1im(i,j,k) + 0.5 * Bd(m) * u_vec(m)
       end do
    end if
    !------------------------------------------


    ! if( abs(y(i,j,k)) < 1.0d-05 .and. abs(z(i,j,k)) < 1.0d-05 ) then

    !    write(*,*) 'i, j, k, x = ', i, j, k, x(i,j,k)
    !    write(*,*) 'psi4re     = ', psi4re(i,j,k)
    !    write(*,*) 'psi4im     = ', psi4im(i,j,k)
    !    write(*,*) 'lA         = ', lA
    !    write(*,*) 'Ed         = ', Ed
    !    write(*,*) 'Bd         = ', Bd
    !    write(*,*) 'detgd      = ', detgd

    ! end if

  end do
  end do
  end do

end subroutine NPEM_CalcPhiGF

