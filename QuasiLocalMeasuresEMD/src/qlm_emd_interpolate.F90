#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



! TODO:
! instead of interpolating to the symmetry points, copy them



! A convenient shortcut
#define P(x) CCTK_PointerTo(x)



subroutine qlm_emd_interpolate (CCTK_ARGUMENTS, hn)
  use cctk
  use qlm_emd_variables
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: hn
  
  CCTK_INT,  parameter :: izero = 0
  integer,   parameter :: ik = kind(izero)
  integer,   parameter :: sk = kind(interpolator)
  CCTK_REAL, parameter :: one = 1
  
  CCTK_REAL, parameter :: poison_value = -42
  
  integer      :: len_coordsystem
  integer      :: len_interpolator
  integer      :: len_interpolator_options
  
  character    :: fort_coordsystem*100
  character    :: fort_interpolator*100
  character    :: fort_interpolator_options*1000
  
  integer      :: nvars
  
  integer      :: coord_handle
  integer      :: interp_handle
  integer      :: options_table
  
  integer      :: ninputs
  integer      :: noutputs
  
  CCTK_REAL, allocatable :: xcoord(:,:)
  CCTK_REAL, allocatable :: ycoord(:,:)
  CCTK_REAL, allocatable :: zcoord(:,:)
  
  integer      :: ind_gxx, ind_gxy, ind_gxz, ind_gyy, ind_gyz, ind_gzz
  integer      :: ind_kxx, ind_kxy, ind_kxz, ind_kyy, ind_kyz, ind_kzz
  integer      :: ind_alpha
  integer      :: ind_betax, ind_betay, ind_betaz
  integer      :: ind_ttt
  integer      :: ind_ttx, ind_tty, ind_ttz
  integer      :: ind_txx, ind_txy, ind_txz, ind_tyy, ind_tyz, ind_tzz
  integer      :: ind_ex, ind_ey, ind_ez
  
  integer      :: coord_type
  CCTK_POINTER :: coords(3)
  CCTK_INT     :: inputs(29)
  CCTK_INT     :: output_types(101)
  CCTK_POINTER :: outputs(101)
  CCTK_INT     :: operand_indices(101)
  CCTK_INT     :: operation_codes(101)
  integer      :: npoints
  
  character    :: msg*1000
  
  integer      :: ni, nj
  
  integer      :: ierr
  
  
  
  if (veryverbose/=0) then
     call CCTK_INFO ("Interpolating 3d grid functions")
  end if
  
  
  
  if (shift_state==0) then
     call CCTK_WARN (0, "The shift must have storage")
  end if
  
!!$  if (stress_energy_state==0) then
!!$     call CCTK_WARN (0, "The stress-energy tensor must have storage")
!!$  end if
  
  
  
  ! Get coordinate system
  call CCTK_FortranString &
       (len_coordsystem, int(coordsystem,sk), fort_coordsystem)
  call CCTK_CoordSystemHandle (coord_handle, fort_coordsystem)
  if (coord_handle<0) then
     write (msg, '("The coordinate system """, a, """ does not exist")') &
          trim(fort_coordsystem)
     call CCTK_WARN (0, msg)
  end if
  
  ! Get interpolator
  call CCTK_FortranString &
       (len_interpolator, int(interpolator,sk), fort_interpolator)
  call CCTK_InterpHandle (interp_handle, fort_interpolator)
  if (interp_handle<0) then
     write (msg, '("The interpolator """,a,""" does not exist")') &
          trim(fort_interpolator)
     call CCTK_WARN (0, msg)
  end if
  
  ! Get interpolator options
  call CCTK_FortranString &
       (len_interpolator_options, int(interpolator_options,sk), &
       fort_interpolator_options)
  call Util_TableCreateFromString (options_table, fort_interpolator_options)
  if (options_table<0) then
     write (msg, '("The interpolator_options """,a,""" have a wrong syntax")') &
          trim(fort_interpolator_options)
     call CCTK_WARN (0, msg)
  end if
  
  
  
  if (hn > 0) then
     
     ni = qlm_emd_ntheta(hn)
     nj = qlm_emd_nphi(hn)
     
     allocate (xcoord(ni,nj))
     allocate (ycoord(ni,nj))
     allocate (zcoord(ni,nj))
     
     xcoord(:,:) = qlm_emd_x(:ni,:nj,hn)
     ycoord(:,:) = qlm_emd_y(:ni,:nj,hn)
     zcoord(:,:) = qlm_emd_z(:ni,:nj,hn)
     
  end if
  
  
  
  ! TODO: check the excision mask
  
  ! Get variable indices
  call CCTK_VarIndex (ind_gxx  , "ADMBase::gxx"   )
  call CCTK_VarIndex (ind_gxy  , "ADMBase::gxy"   )
  call CCTK_VarIndex (ind_gxz  , "ADMBase::gxz"   )
  call CCTK_VarIndex (ind_gyy  , "ADMBase::gyy"   )
  call CCTK_VarIndex (ind_gyz  , "ADMBase::gyz"   )
  call CCTK_VarIndex (ind_gzz  , "ADMBase::gzz"   )
  call CCTK_VarIndex (ind_kxx  , "ADMBase::kxx"   )
  call CCTK_VarIndex (ind_kxy  , "ADMBase::kxy"   )
  call CCTK_VarIndex (ind_kxz  , "ADMBase::kxz"   )
  call CCTK_VarIndex (ind_kyy  , "ADMBase::kyy"   )
  call CCTK_VarIndex (ind_kyz  , "ADMBase::kyz"   )
  call CCTK_VarIndex (ind_kzz  , "ADMBase::kzz"   )
  call CCTK_VarIndex (ind_alpha, "ADMBase::alp"   )
  call CCTK_VarIndex (ind_betax, "ADMBase::betax" )
  call CCTK_VarIndex (ind_betay, "ADMBase::betay" )
  call CCTK_VarIndex (ind_betaz, "ADMBase::betaz" )
  if (stress_energy_state /= 0) then
     call CCTK_VarIndex (ind_ttt  , "TmunuBase::eTtt")
     call CCTK_VarIndex (ind_ttx  , "TmunuBase::eTtx")
     call CCTK_VarIndex (ind_tty  , "TmunuBase::eTty")
     call CCTK_VarIndex (ind_ttz  , "TmunuBase::eTtz")
     call CCTK_VarIndex (ind_txx  , "TmunuBase::eTxx")
     call CCTK_VarIndex (ind_txy  , "TmunuBase::eTxy")
     call CCTK_VarIndex (ind_txz  , "TmunuBase::eTxz")
     call CCTK_VarIndex (ind_tyy  , "TmunuBase::eTyy")
     call CCTK_VarIndex (ind_tyz  , "TmunuBase::eTyz")
     call CCTK_VarIndex (ind_tzz  , "TmunuBase::eTzz")
  else
     ind_ttt = -1
     ind_ttx = -1
     ind_tty = -1
     ind_ttz = -1
     ind_txx = -1
     ind_txy = -1
     ind_txz = -1
     ind_tyy = -1
     ind_tyz = -1
     ind_tzz = -1
  end if
  if (calc_charge /= 0) then
     call CCTK_VarIndex (ind_ex  , "ProcaBase::Ex")
     call CCTK_VarIndex (ind_ey  , "ProcaBase::Ey")
     call CCTK_VarIndex (ind_ez  , "ProcaBase::Ez")
  else
     ind_ex = -1
     ind_ey = -1
     ind_ez = -1
  end if
  
  
  
  ! Set up the interpolator arguments
  coord_type = CCTK_VARIABLE_REAL
  if (hn > 0) then
     npoints = ni * nj
     coords(:) = (/ P(xcoord), P(ycoord), P(zcoord) /)
  else
     npoints = 0
     coords(:) = CCTK_NullPointer()
  end if
  
  inputs = (/ &
       ind_gxx, ind_gxy, ind_gxz, ind_gyy, ind_gyz, ind_gzz, &
       ind_kxx, ind_kxy, ind_kxz, ind_kyy, ind_kyz, ind_kzz, &
       ind_alpha, &
       ind_betax, ind_betay, ind_betaz, &
       ind_ttt, &
       ind_ttx, ind_tty, ind_ttz, &
       ind_txx, ind_txy, ind_txz, ind_tyy, ind_tyz, ind_tzz, &
       ind_ex, ind_ey, ind_ez /)
  
  call CCTK_NumVars (nvars)
  if (nvars < 0) call CCTK_WARN (0, "internal error")
  if (any(inputs /= -1 .and. (inputs < 0 .or. inputs >= nvars))) then
     call CCTK_WARN (0, "internal error")
  end if
  
  operand_indices = (/ &
       00, 01, 02, 03, 04, 05, & ! g_ij
       00, 01, 02, 03, 04, 05, & ! g_ij,k
       00, 01, 02, 03, 04, 05, &
       00, 01, 02, 03, 04, 05, &
       00, 01, 02, 03, 04, 05, & ! g_ij,kl
       00, 01, 02, 03, 04, 05, &
       00, 01, 02, 03, 04, 05, &
       00, 01, 02, 03, 04, 05, &
       00, 01, 02, 03, 04, 05, &
       00, 01, 02, 03, 04, 05, &
       06, 07, 08, 09, 10, 11, & ! K_ij
       06, 07, 08, 09, 10, 11, & ! K_ij,k
       06, 07, 08, 09, 10, 11, &
       06, 07, 08, 09, 10, 11, &
       12, &                     ! alp
       13, 14, 15, &             ! beta^i
       16, &                     ! T_tt
       17, 18, 19, &             ! T_ti
       20, 21, 22, 23, 24, 25, & ! T_ij
       26, 27, 28 /)             ! E^i
  
  operation_codes = (/ &
       0, 0, 0, 0, 0, 0, &      ! g_ij
       1, 1, 1, 1, 1, 1, &      ! g_ij,k
       2, 2, 2, 2, 2, 2, &
       3, 3, 3, 3, 3, 3, &
       11, 11, 11, 11, 11, 11, & ! g_ij,kl
       12, 12, 12, 12, 12, 12, &
       13, 13, 13, 13, 13, 13, &
       22, 22, 22, 22, 22, 22, &
       23, 23, 23, 23, 23, 23, &
       33, 33, 33, 33, 33, 33, &
       0, 0, 0, 0, 0, 0, &      ! K_ij
       1, 1, 1, 1, 1, 1, &      ! K_ij,k
       2, 2, 2, 2, 2, 2, &
       3, 3, 3, 3, 3, 3, &
       0, &                     ! alp
       0, 0, 0, &               ! beta^i
       0, &                     ! T_tt
       0, 0, 0, &               ! T_ti
       0, 0, 0, 0, 0, 0, &      ! T_ij
       0, 0, 0 /)               ! E^i

  output_types(:) = CCTK_VARIABLE_REAL
  if (hn > 0) then
     outputs = (/ &
          P(qlm_emd_gxx), P(qlm_emd_gxy), P(qlm_emd_gxz), P(qlm_emd_gyy), P(qlm_emd_gyz), P(qlm_emd_gzz), &
          P(qlm_emd_dgxxx), P(qlm_emd_dgxyx), P(qlm_emd_dgxzx), P(qlm_emd_dgyyx), P(qlm_emd_dgyzx), P(qlm_emd_dgzzx), &
          P(qlm_emd_dgxxy), P(qlm_emd_dgxyy), P(qlm_emd_dgxzy), P(qlm_emd_dgyyy), P(qlm_emd_dgyzy), P(qlm_emd_dgzzy), &
          P(qlm_emd_dgxxz), P(qlm_emd_dgxyz), P(qlm_emd_dgxzz), P(qlm_emd_dgyyz), P(qlm_emd_dgyzz), P(qlm_emd_dgzzz), &
          P(qlm_emd_ddgxxxx), P(qlm_emd_ddgxyxx), P(qlm_emd_ddgxzxx), P(qlm_emd_ddgyyxx), P(qlm_emd_ddgyzxx), P(qlm_emd_ddgzzxx), &
          P(qlm_emd_ddgxxxy), P(qlm_emd_ddgxyxy), P(qlm_emd_ddgxzxy), P(qlm_emd_ddgyyxy), P(qlm_emd_ddgyzxy), P(qlm_emd_ddgzzxy), &
          P(qlm_emd_ddgxxxz), P(qlm_emd_ddgxyxz), P(qlm_emd_ddgxzxz), P(qlm_emd_ddgyyxz), P(qlm_emd_ddgyzxz), P(qlm_emd_ddgzzxz), &
          P(qlm_emd_ddgxxyy), P(qlm_emd_ddgxyyy), P(qlm_emd_ddgxzyy), P(qlm_emd_ddgyyyy), P(qlm_emd_ddgyzyy), P(qlm_emd_ddgzzyy), &
          P(qlm_emd_ddgxxyz), P(qlm_emd_ddgxyyz), P(qlm_emd_ddgxzyz), P(qlm_emd_ddgyyyz), P(qlm_emd_ddgyzyz), P(qlm_emd_ddgzzyz), &
          P(qlm_emd_ddgxxzz), P(qlm_emd_ddgxyzz), P(qlm_emd_ddgxzzz), P(qlm_emd_ddgyyzz), P(qlm_emd_ddgyzzz), P(qlm_emd_ddgzzzz), &
          P(qlm_emd_kxx), P(qlm_emd_kxy), P(qlm_emd_kxz), P(qlm_emd_kyy), P(qlm_emd_kyz), P(qlm_emd_kzz), &
          P(qlm_emd_dkxxx), P(qlm_emd_dkxyx), P(qlm_emd_dkxzx), P(qlm_emd_dkyyx), P(qlm_emd_dkyzx), P(qlm_emd_dkzzx), &
          P(qlm_emd_dkxxy), P(qlm_emd_dkxyy), P(qlm_emd_dkxzy), P(qlm_emd_dkyyy), P(qlm_emd_dkyzy), P(qlm_emd_dkzzy), &
          P(qlm_emd_dkxxz), P(qlm_emd_dkxyz), P(qlm_emd_dkxzz), P(qlm_emd_dkyyz), P(qlm_emd_dkyzz), P(qlm_emd_dkzzz), &
          P(qlm_emd_alpha), &
          P(qlm_emd_betax), P(qlm_emd_betay), P(qlm_emd_betaz), &
          P(qlm_emd_ttt), &
          P(qlm_emd_ttx), P(qlm_emd_tty), P(qlm_emd_ttz), &
          P(qlm_emd_txx), P(qlm_emd_txy), P(qlm_emd_txz), P(qlm_emd_tyy), P(qlm_emd_tyz), P(qlm_emd_tzz), &
          P(qlm_emd_ex), P(qlm_emd_ey), P(qlm_emd_ez) /)
  else
     outputs(:) = CCTK_NullPointer()
  end if
  
  
  
  ninputs = size(inputs)
  noutputs = size(outputs)
  
#if 0
  ! Poison the output variables
  call poison (qlm_emd_gxx    )
  call poison (qlm_emd_gxy    )
  call poison (qlm_emd_gxz    )
  call poison (qlm_emd_gyy    )
  call poison (qlm_emd_gyz    )
  call poison (qlm_emd_gzz    )
  call poison (qlm_emd_dgxxx  )
  call poison (qlm_emd_dgxyx  )
  call poison (qlm_emd_dgxzx  )
  call poison (qlm_emd_dgyyx  )
  call poison (qlm_emd_dgyzx  )
  call poison (qlm_emd_dgzzx  )
  call poison (qlm_emd_dgxxy  )
  call poison (qlm_emd_dgxyy  )
  call poison (qlm_emd_dgxzy  )
  call poison (qlm_emd_dgyyy  )
  call poison (qlm_emd_dgyzy  )
  call poison (qlm_emd_dgzzy  )
  call poison (qlm_emd_dgxxz  )
  call poison (qlm_emd_dgxyz  )
  call poison (qlm_emd_dgxzz  )
  call poison (qlm_emd_dgyyz  )
  call poison (qlm_emd_dgyzz  )
  call poison (qlm_emd_dgzzz  )
  call poison (qlm_emd_ddgxxxx)
  call poison (qlm_emd_ddgxyxx)
  call poison (qlm_emd_ddgxzxx)
  call poison (qlm_emd_ddgyyxx)
  call poison (qlm_emd_ddgyzxx)
  call poison (qlm_emd_ddgzzxx)
  call poison (qlm_emd_ddgxxxy)
  call poison (qlm_emd_ddgxyxy)
  call poison (qlm_emd_ddgxzxy)
  call poison (qlm_emd_ddgyyxy)
  call poison (qlm_emd_ddgyzxy)
  call poison (qlm_emd_ddgzzxy)
  call poison (qlm_emd_ddgxxxz)
  call poison (qlm_emd_ddgxyxz)
  call poison (qlm_emd_ddgxzxz)
  call poison (qlm_emd_ddgyyxz)
  call poison (qlm_emd_ddgyzxz)
  call poison (qlm_emd_ddgzzxz)
  call poison (qlm_emd_ddgxxyy)
  call poison (qlm_emd_ddgxyyy)
  call poison (qlm_emd_ddgxzyy)
  call poison (qlm_emd_ddgyyyy)
  call poison (qlm_emd_ddgyzyy)
  call poison (qlm_emd_ddgzzyy)
  call poison (qlm_emd_ddgxxyz)
  call poison (qlm_emd_ddgxyyz)
  call poison (qlm_emd_ddgxzyz)
  call poison (qlm_emd_ddgyyyz)
  call poison (qlm_emd_ddgyzyz)
  call poison (qlm_emd_ddgzzyz)
  call poison (qlm_emd_ddgxxzz)
  call poison (qlm_emd_ddgxyzz)
  call poison (qlm_emd_ddgxzzz)
  call poison (qlm_emd_ddgyyzz)
  call poison (qlm_emd_ddgyzzz)
  call poison (qlm_emd_ddgzzzz)
  call poison (qlm_emd_kxx    )
  call poison (qlm_emd_kxy    )
  call poison (qlm_emd_kxz    )
  call poison (qlm_emd_kyy    )
  call poison (qlm_emd_kyz    )
  call poison (qlm_emd_kzz    )
  call poison (qlm_emd_dkxxx  )
  call poison (qlm_emd_dkxyx  )
  call poison (qlm_emd_dkxzx  )
  call poison (qlm_emd_dkyyx  )
  call poison (qlm_emd_dkyzx  )
  call poison (qlm_emd_dkzzx  )
  call poison (qlm_emd_dkxxy  )
  call poison (qlm_emd_dkxyy  )
  call poison (qlm_emd_dkxzy  )
  call poison (qlm_emd_dkyyy  )
  call poison (qlm_emd_dkyzy  )
  call poison (qlm_emd_dkzzy  )
  call poison (qlm_emd_dkxxz  )
  call poison (qlm_emd_dkxyz  )
  call poison (qlm_emd_dkxzz  )
  call poison (qlm_emd_dkyyz  )
  call poison (qlm_emd_dkyzz  )
  call poison (qlm_emd_dkzzz  )
  call poison (qlm_emd_alpha  )
  call poison (qlm_emd_betax  )
  call poison (qlm_emd_betay  )
  call poison (qlm_emd_betaz  )
  call poison (qlm_emd_ttt    )
  call poison (qlm_emd_ttx    )
  call poison (qlm_emd_tty    )
  call poison (qlm_emd_ttz    )
  call poison (qlm_emd_txx    )
  call poison (qlm_emd_txy    )
  call poison (qlm_emd_txz    )
  call poison (qlm_emd_tyy    )
  call poison (qlm_emd_tyz    )
  call poison (qlm_emd_tzz    )
  call poison (qlm_emd_ex     )
  call poison (qlm_emd_ey     )
  call poison (qlm_emd_ez     )
#endif
  
call CCTK_INFO ("Finish poison")
  
  ! Call the interpolator
  call Util_TableSetIntArray &
       (ierr, options_table, noutputs, &
       operand_indices, "operand_indices")
  if (ierr /= 0) call CCTK_WARN (0, "internal error")
  call Util_TableSetIntArray &
       (ierr, options_table, noutputs, &
       operation_codes, "operation_codes")
  if (ierr /= 0) call CCTK_WARN (0, "internal error")
  call CCTK_INFO ("Util_TableSetIntArray")

  write (msg, '("   ninputs:         ",g16.6)') ninputs
  call CCTK_INFO (msg)
  write (msg, '("   noutputs:         ",g16.6)') noutputs
  call CCTK_INFO (msg)
  write (msg, '("   outputs first:         ",g16.6)') outputs(1)
  call CCTK_INFO (msg)
  write (msg, '("   outputs last:         ",g16.6)') outputs(noutputs)
  call CCTK_INFO (msg)
  call CCTK_InterpGridArrays &
       (ierr, cctkGH, 3, &
       interp_handle, options_table, coord_handle, &
       npoints, coord_type, coords, &
       ninputs, inputs, &
       noutputs, output_types, outputs)
   call CCTK_INFO ("CCTK_InterpGridArrays")
  
  if (ierr /= 0) then
     if (hn > 0) then
        qlm_emd_calc_error(hn) = 1
     end if
     call CCTK_WARN (1, "Interpolator failed")
     return
  end if
  
  
  
  ! Unpack the variables
  if (hn > 0) then
     
     call unpack (qlm_emd_gxx    , ni, nj)
     call unpack (qlm_emd_gxy    , ni, nj)
     call unpack (qlm_emd_gxz    , ni, nj)
     call unpack (qlm_emd_gyy    , ni, nj)
     call unpack (qlm_emd_gyz    , ni, nj)
     call unpack (qlm_emd_gzz    , ni, nj)
     call unpack (qlm_emd_dgxxx  , ni, nj)
     call unpack (qlm_emd_dgxyx  , ni, nj)
     call unpack (qlm_emd_dgxzx  , ni, nj)
     call unpack (qlm_emd_dgyyx  , ni, nj)
     call unpack (qlm_emd_dgyzx  , ni, nj)
     call unpack (qlm_emd_dgzzx  , ni, nj)
     call unpack (qlm_emd_dgxxy  , ni, nj)
     call unpack (qlm_emd_dgxyy  , ni, nj)
     call unpack (qlm_emd_dgxzy  , ni, nj)
     call unpack (qlm_emd_dgyyy  , ni, nj)
     call unpack (qlm_emd_dgyzy  , ni, nj)
     call unpack (qlm_emd_dgzzy  , ni, nj)
     call unpack (qlm_emd_dgxxz  , ni, nj)
     call unpack (qlm_emd_dgxyz  , ni, nj)
     call unpack (qlm_emd_dgxzz  , ni, nj)
     call unpack (qlm_emd_dgyyz  , ni, nj)
     call unpack (qlm_emd_dgyzz  , ni, nj)
     call unpack (qlm_emd_dgzzz  , ni, nj)
     call unpack (qlm_emd_ddgxxxx, ni, nj)
     call unpack (qlm_emd_ddgxyxx, ni, nj)
     call unpack (qlm_emd_ddgxzxx, ni, nj)
     call unpack (qlm_emd_ddgyyxx, ni, nj)
     call unpack (qlm_emd_ddgyzxx, ni, nj)
     call unpack (qlm_emd_ddgzzxx, ni, nj)
     call unpack (qlm_emd_ddgxxxy, ni, nj)
     call unpack (qlm_emd_ddgxyxy, ni, nj)
     call unpack (qlm_emd_ddgxzxy, ni, nj)
     call unpack (qlm_emd_ddgyyxy, ni, nj)
     call unpack (qlm_emd_ddgyzxy, ni, nj)
     call unpack (qlm_emd_ddgzzxy, ni, nj)
     call unpack (qlm_emd_ddgxxxz, ni, nj)
     call unpack (qlm_emd_ddgxyxz, ni, nj)
     call unpack (qlm_emd_ddgxzxz, ni, nj)
     call unpack (qlm_emd_ddgyyxz, ni, nj)
     call unpack (qlm_emd_ddgyzxz, ni, nj)
     call unpack (qlm_emd_ddgzzxz, ni, nj)
     call unpack (qlm_emd_ddgxxyy, ni, nj)
     call unpack (qlm_emd_ddgxyyy, ni, nj)
     call unpack (qlm_emd_ddgxzyy, ni, nj)
     call unpack (qlm_emd_ddgyyyy, ni, nj)
     call unpack (qlm_emd_ddgyzyy, ni, nj)
     call unpack (qlm_emd_ddgzzyy, ni, nj)
     call unpack (qlm_emd_ddgxxyz, ni, nj)
     call unpack (qlm_emd_ddgxyyz, ni, nj)
     call unpack (qlm_emd_ddgxzyz, ni, nj)
     call unpack (qlm_emd_ddgyyyz, ni, nj)
     call unpack (qlm_emd_ddgyzyz, ni, nj)
     call unpack (qlm_emd_ddgzzyz, ni, nj)
     call unpack (qlm_emd_ddgxxzz, ni, nj)
     call unpack (qlm_emd_ddgxyzz, ni, nj)
     call unpack (qlm_emd_ddgxzzz, ni, nj)
     call unpack (qlm_emd_ddgyyzz, ni, nj)
     call unpack (qlm_emd_ddgyzzz, ni, nj)
     call unpack (qlm_emd_ddgzzzz, ni, nj)
     call unpack (qlm_emd_kxx    , ni, nj)
     call unpack (qlm_emd_kxy    , ni, nj)
     call unpack (qlm_emd_kxz    , ni, nj)
     call unpack (qlm_emd_kyy    , ni, nj)
     call unpack (qlm_emd_kyz    , ni, nj)
     call unpack (qlm_emd_kzz    , ni, nj)
     call unpack (qlm_emd_dkxxx  , ni, nj)
     call unpack (qlm_emd_dkxyx  , ni, nj)
     call unpack (qlm_emd_dkxzx  , ni, nj)
     call unpack (qlm_emd_dkyyx  , ni, nj)
     call unpack (qlm_emd_dkyzx  , ni, nj)
     call unpack (qlm_emd_dkzzx  , ni, nj)
     call unpack (qlm_emd_dkxxy  , ni, nj)
     call unpack (qlm_emd_dkxyy  , ni, nj)
     call unpack (qlm_emd_dkxzy  , ni, nj)
     call unpack (qlm_emd_dkyyy  , ni, nj)
     call unpack (qlm_emd_dkyzy  , ni, nj)
     call unpack (qlm_emd_dkzzy  , ni, nj)
     call unpack (qlm_emd_dkxxz  , ni, nj)
     call unpack (qlm_emd_dkxyz  , ni, nj)
     call unpack (qlm_emd_dkxzz  , ni, nj)
     call unpack (qlm_emd_dkyyz  , ni, nj)
     call unpack (qlm_emd_dkyzz  , ni, nj)
     call unpack (qlm_emd_dkzzz  , ni, nj)
     call unpack (qlm_emd_alpha  , ni, nj)
     call unpack (qlm_emd_betax  , ni, nj)
     call unpack (qlm_emd_betay  , ni, nj)
     call unpack (qlm_emd_betaz  , ni, nj)
     if (stress_energy_state /= 0) then
        call unpack (qlm_emd_ttt    , ni, nj)
        call unpack (qlm_emd_ttx    , ni, nj)
        call unpack (qlm_emd_tty    , ni, nj)
        call unpack (qlm_emd_ttz    , ni, nj)
        call unpack (qlm_emd_txx    , ni, nj)
        call unpack (qlm_emd_txy    , ni, nj)
        call unpack (qlm_emd_txz    , ni, nj)
        call unpack (qlm_emd_tyy    , ni, nj)
        call unpack (qlm_emd_tyz    , ni, nj)
        call unpack (qlm_emd_tzz    , ni, nj)
     else
        qlm_emd_ttt = 0
        qlm_emd_ttx = 0
        qlm_emd_tty = 0
        qlm_emd_ttz = 0
        qlm_emd_txx = 0
        qlm_emd_txy = 0
        qlm_emd_txz = 0
        qlm_emd_tyy = 0
        qlm_emd_tyz = 0
        qlm_emd_tzz = 0
     end if
     if (calc_charge /= 0) then
        call unpack (qlm_emd_ex    , ni, nj)
        call unpack (qlm_emd_ey    , ni, nj)
        call unpack (qlm_emd_ez    , ni, nj)
     else
        qlm_emd_ex = 0
        qlm_emd_ey = 0
        qlm_emd_ez = 0
     end if
     
   call CCTK_INFO ("Finish unpack")
     
     
#if 0
     ! Check for poison
     call poison_check (qlm_emd_gxx    , "qlm_emd_gxx    ")
     call poison_check (qlm_emd_gxy    , "qlm_emd_gxy    ")
     call poison_check (qlm_emd_gxz    , "qlm_emd_gxz    ")
     call poison_check (qlm_emd_gyy    , "qlm_emd_gyy    ")
     call poison_check (qlm_emd_gyz    , "qlm_emd_gyz    ")
     call poison_check (qlm_emd_gzz    , "qlm_emd_gzz    ")
     call poison_check (qlm_emd_dgxxx  , "qlm_emd_dgxxx  ")
     call poison_check (qlm_emd_dgxyx  , "qlm_emd_dgxyx  ")
     call poison_check (qlm_emd_dgxzx  , "qlm_emd_dgxzx  ")
     call poison_check (qlm_emd_dgyyx  , "qlm_emd_dgyyx  ")
     call poison_check (qlm_emd_dgyzx  , "qlm_emd_dgyzx  ")
     call poison_check (qlm_emd_dgzzx  , "qlm_emd_dgzzx  ")
     call poison_check (qlm_emd_dgxxy  , "qlm_emd_dgxxy  ")
     call poison_check (qlm_emd_dgxyy  , "qlm_emd_dgxyy  ")
     call poison_check (qlm_emd_dgxzy  , "qlm_emd_dgxzy  ")
     call poison_check (qlm_emd_dgyyy  , "qlm_emd_dgyyy  ")
     call poison_check (qlm_emd_dgyzy  , "qlm_emd_dgyzy  ")
     call poison_check (qlm_emd_dgzzy  , "qlm_emd_dgzzy  ")
     call poison_check (qlm_emd_dgxxz  , "qlm_emd_dgxxz  ")
     call poison_check (qlm_emd_dgxyz  , "qlm_emd_dgxyz  ")
     call poison_check (qlm_emd_dgxzz  , "qlm_emd_dgxzz  ")
     call poison_check (qlm_emd_dgyyz  , "qlm_emd_dgyyz  ")
     call poison_check (qlm_emd_dgyzz  , "qlm_emd_dgyzz  ")
     call poison_check (qlm_emd_dgzzz  , "qlm_emd_dgzzz  ")
     call poison_check (qlm_emd_ddgxxxx, "qlm_emd_ddgxxxx")
     call poison_check (qlm_emd_ddgxyxx, "qlm_emd_ddgxyxx")
     call poison_check (qlm_emd_ddgxzxx, "qlm_emd_ddgxzxx")
     call poison_check (qlm_emd_ddgyyxx, "qlm_emd_ddgyyxx")
     call poison_check (qlm_emd_ddgyzxx, "qlm_emd_ddgyzxx")
     call poison_check (qlm_emd_ddgzzxx, "qlm_emd_ddgzzxx")
     call poison_check (qlm_emd_ddgxxxy, "qlm_emd_ddgxxxy")
     call poison_check (qlm_emd_ddgxyxy, "qlm_emd_ddgxyxy")
     call poison_check (qlm_emd_ddgxzxy, "qlm_emd_ddgxzxy")
     call poison_check (qlm_emd_ddgyyxy, "qlm_emd_ddgyyxy")
     call poison_check (qlm_emd_ddgyzxy, "qlm_emd_ddgyzxy")
     call poison_check (qlm_emd_ddgzzxy, "qlm_emd_ddgzzxy")
     call poison_check (qlm_emd_ddgxxxz, "qlm_emd_ddgxxxz")
     call poison_check (qlm_emd_ddgxyxz, "qlm_emd_ddgxyxz")
     call poison_check (qlm_emd_ddgxzxz, "qlm_emd_ddgxzxz")
     call poison_check (qlm_emd_ddgyyxz, "qlm_emd_ddgyyxz")
     call poison_check (qlm_emd_ddgyzxz, "qlm_emd_ddgyzxz")
     call poison_check (qlm_emd_ddgzzxz, "qlm_emd_ddgzzxz")
     call poison_check (qlm_emd_ddgxxyy, "qlm_emd_ddgxxyy")
     call poison_check (qlm_emd_ddgxyyy, "qlm_emd_ddgxyyy")
     call poison_check (qlm_emd_ddgxzyy, "qlm_emd_ddgxzyy")
     call poison_check (qlm_emd_ddgyyyy, "qlm_emd_ddgyyyy")
     call poison_check (qlm_emd_ddgyzyy, "qlm_emd_ddgyzyy")
     call poison_check (qlm_emd_ddgzzyy, "qlm_emd_ddgzzyy")
     call poison_check (qlm_emd_ddgxxyz, "qlm_emd_ddgxxyz")
     call poison_check (qlm_emd_ddgxyyz, "qlm_emd_ddgxyyz")
     call poison_check (qlm_emd_ddgxzyz, "qlm_emd_ddgxzyz")
     call poison_check (qlm_emd_ddgyyyz, "qlm_emd_ddgyyyz")
     call poison_check (qlm_emd_ddgyzyz, "qlm_emd_ddgyzyz")
     call poison_check (qlm_emd_ddgzzyz, "qlm_emd_ddgzzyz")
     call poison_check (qlm_emd_ddgxxzz, "qlm_emd_ddgxxzz")
     call poison_check (qlm_emd_ddgxyzz, "qlm_emd_ddgxyzz")
     call poison_check (qlm_emd_ddgxzzz, "qlm_emd_ddgxzzz")
     call poison_check (qlm_emd_ddgyyzz, "qlm_emd_ddgyyzz")
     call poison_check (qlm_emd_ddgyzzz, "qlm_emd_ddgyzzz")
     call poison_check (qlm_emd_ddgzzzz, "qlm_emd_ddgzzzz")
     call poison_check (qlm_emd_kxx    , "qlm_emd_kxx    ")
     call poison_check (qlm_emd_kxy    , "qlm_emd_kxy    ")
     call poison_check (qlm_emd_kxz    , "qlm_emd_kxz    ")
     call poison_check (qlm_emd_kyy    , "qlm_emd_kyy    ")
     call poison_check (qlm_emd_kyz    , "qlm_emd_kyz    ")
     call poison_check (qlm_emd_kzz    , "qlm_emd_kzz    ")
     call poison_check (qlm_emd_dkxxx  , "qlm_emd_dkxxx  ")
     call poison_check (qlm_emd_dkxyx  , "qlm_emd_dkxyx  ")
     call poison_check (qlm_emd_dkxzx  , "qlm_emd_dkxzx  ")
     call poison_check (qlm_emd_dkyyx  , "qlm_emd_dkyyx  ")
     call poison_check (qlm_emd_dkyzx  , "qlm_emd_dkyzx  ")
     call poison_check (qlm_emd_dkzzx  , "qlm_emd_dkzzx  ")
     call poison_check (qlm_emd_dkxxy  , "qlm_emd_dkxxy  ")
     call poison_check (qlm_emd_dkxyy  , "qlm_emd_dkxyy  ")
     call poison_check (qlm_emd_dkxzy  , "qlm_emd_dkxzy  ")
     call poison_check (qlm_emd_dkyyy  , "qlm_emd_dkyyy  ")
     call poison_check (qlm_emd_dkyzy  , "qlm_emd_dkyzy  ")
     call poison_check (qlm_emd_dkzzy  , "qlm_emd_dkzzy  ")
     call poison_check (qlm_emd_dkxxz  , "qlm_emd_dkxxz  ")
     call poison_check (qlm_emd_dkxyz  , "qlm_emd_dkxyz  ")
     call poison_check (qlm_emd_dkxzz  , "qlm_emd_dkxzz  ")
     call poison_check (qlm_emd_dkyyz  , "qlm_emd_dkyyz  ")
     call poison_check (qlm_emd_dkyzz  , "qlm_emd_dkyzz  ")
     call poison_check (qlm_emd_dkzzz  , "qlm_emd_dkzzz  ")
     call poison_check (qlm_emd_alpha  , "qlm_emd_alpha  ")
     call poison_check (qlm_emd_betax  , "qlm_emd_betax  ")
     call poison_check (qlm_emd_betay  , "qlm_emd_betay  ")
     call poison_check (qlm_emd_betaz  , "qlm_emd_betaz  ")
     call poison_check (qlm_emd_ttt    , "qlm_emd_ttt    ")
     call poison_check (qlm_emd_ttx    , "qlm_emd_ttx    ")
     call poison_check (qlm_emd_tty    , "qlm_emd_tty    ")
     call poison_check (qlm_emd_ttz    , "qlm_emd_ttz    ")
     call poison_check (qlm_emd_txx    , "qlm_emd_txx    ")
     call poison_check (qlm_emd_txy    , "qlm_emd_txy    ")
     call poison_check (qlm_emd_txz    , "qlm_emd_txz    ")
     call poison_check (qlm_emd_tyy    , "qlm_emd_tyy    ")
     call poison_check (qlm_emd_tyz    , "qlm_emd_tyz    ")
     call poison_check (qlm_emd_tzz    , "qlm_emd_tzz    ")
     call poison_check (qlm_emd_ex     , "qlm_emd_ex     ")
     call poison_check (qlm_emd_ey     , "qlm_emd_ey     ")
     call poison_check (qlm_emd_ez     , "qlm_emd_ez     ")
#endif
     
  end if
  
  call CCTK_INFO ("Finish poison_check")
  
  ! Free interpolator options
  call Util_TableDestroy (ierr, options_table)
  
  
  
  if (hn > 0) then
     
     qlm_emd_have_valid_data(hn) = 1
     
     deallocate (xcoord)
     deallocate (ycoord)
     deallocate (zcoord)
     
  end if

  call CCTK_INFO ("Finish interpolate")
  
  
contains
  
  subroutine pack (arr, ni, nj)
    integer,   intent(in)    :: ni, nj
    CCTK_REAL, intent(inout) :: arr(:,:)
    CCTK_REAL :: tmp(ni,nj)
    tmp(:,:) = arr(:ni, :nj)
    call copy (arr, tmp, size(tmp))
  end subroutine pack
  
  subroutine unpack (arr, ni, nj)
    integer,   intent(in)    :: ni, nj
    CCTK_REAL, intent(inout) :: arr(:,:)
    CCTK_REAL :: tmp(ni,nj)
    call copy (tmp, arr, size(tmp))
    arr(:ni, :nj) = tmp(:,:)
    arr(ni+1:, :nj) = 0
    arr(:, nj+1:) = 0
  end subroutine unpack
  
  subroutine copy (a, b, n)
    integer,   intent(in)  :: n
    CCTK_REAL, intent(out) :: a(n)
    CCTK_REAL, intent(in)  :: b(n)
    a = b
  end subroutine copy
  
  subroutine poison (arr)
    CCTK_REAL, intent(out) :: arr(:,:)
    arr = poison_value
  end subroutine poison
  
  subroutine poison_check (arr, name)
    CCTK_REAL,    intent(in) :: arr(:,:)
    character(*), intent(in) :: name
    character*1000 :: msg
!!$    integer        :: i, j
    if (any(arr==poison_value)) then
       write (msg, '("Poison found in ",a)') trim(name)
       call CCTK_WARN (CCTK_WARN_ALERT, msg)
!!$       do j=1,size(arr,2)
!!$          do i=1,size(arr,1)
!!$             print '(2i6)', i,j
!!$          end do
!!$       end do
    end if
  end subroutine poison_check
  
end subroutine qlm_emd_interpolate
