!!  added by Xiantao Fan
!!  Copyright (C) 2009-2017  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!
!!  lesgo is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  lesgo is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with lesgo.  If not, see <http://www.gnu.org/licenses/>.
!!
!=================added by Xiantao Fan to calculate the aerodynamics of trees
!*******************************************************************************
module stat_treesdefs
!*******************************************************************************
use types, only : rprec
use param, only : nx, ny, nz, lh

save
public

type point_t
    integer :: istart, jstart, kstart, coord
    real(rprec) :: xdiff, ydiff, zdiff
    integer :: fid
end type point_t

type plane_t
    integer :: istart
    real(rprec) :: ldiff
end type plane_t

type zplane_t
    integer :: istart, coord
    real(rprec) :: ldiff
end type zplane_t

type rs_t
    real(rprec) :: up2, vp2, wp2, upvp, upwp, vpwp
end type rs_t

!add by Mingwei for pressure related term avg
type ps_t
    real(rprec) :: pp2, ppup, ppvp, ppwp
end type ps_t

type spectra_t
    real(rprec), dimension(:), allocatable :: power
    integer :: istart, coord
    real(rprec) :: ldiff
end type spectra_t

real(rprec) :: spectra_total_time
real(rprec) :: tavg_total_time
#ifdef PPOUTPUT_EXTRA
real(rprec) :: tavg_total_time_sgs
#endif
! Time between calls of tavg_compute, built by summing dt
real(rprec) :: tavg_dt
! Switch for determining if time averaging has been initialized
logical :: tavg_initialized = .false.

!  Sums performed over time
type tavg_t
    real(rprec) :: u, v, w, u_w, v_w, w_uv
    real(rprec) :: u2, v2, w2, uv, uw, vw
    ! real(rprec) :: dudz, dvdz
    real(rprec) :: txx, tyy, tzz, txy, txz, tyz
    real(rprec) :: fx, fy, fz
    real(rprec) :: cs_opt2
    real(rprec) :: p                !add by Mingwei
    real(rprec) :: p2,pu,pv,pw      !add by Mingwei
end type tavg_t

!  Sums performed over time (for subgrid variables)
#ifdef PPOUTPUT_EXTRA
type tavg_sgs_t
    real(rprec) :: Nu_t
end type tavg_sgs_t
#endif

! Types for including trees as drag force
#ifdef PPTREES

! Indicator function calculator
type tree_ind_func_t
    real(rprec), dimension(:), allocatable :: r
    real(rprec), dimension(:), allocatable :: R23
    real(rprec) :: sqrt6overdelta, t_half
contains
    procedure, public :: init
    procedure, public :: val
end type tree_ind_func_t

! Single trees
type tree_t
    real(rprec) :: xloc, yloc, h_all, Rt_all,Lt_all, zm_all, am_all
    ! term used for volume correction
    real(rprec) :: vol_c
    ! efecitive coefficient for trees in three directions
    real(rprec) :: p_x, p_y, p_z

    ! number of nodes associated with each tree
    integer :: num_nodes
    ! location of tree center (local k)
    integer :: icp, jcp, kcp
    ! true if the center is in the processor
    logical :: center_in_proc
    ! drag coefficient
    real(rprec) :: Cd
    ! running time-average of mean disk velocity
    real(rprec) :: u_d, u_d_T
    ! normal force on tree
    real(rprec) :: f_n_x, f_n_y, f_n_z
    ! (nx,ny,nz) of unit normal for each tree
    real(rprec), dimension(3) :: nhat
    ! indicator function - weighting of each node
    real(rprec), dimension(500000) :: ind
    ! object to calculate indicator function
    type(tree_ind_func_t) :: tree_ind_func
    ! (i,j,k) of each included node
    integer, dimension(500000,3) :: nodes
    ! search area for nearby nodes
    integer, dimension(6) :: nodes_max
end type tree_t

! A collection of trees
type tree_farm_tr
    type(tree_t), pointer, dimension(:) :: tree
end type tree_farm_tr

! The wind farm
type(tree_farm_tr) :: tree_farm

#endif

! Create types for outputting data (instantaneous or averaged)
type(point_t), allocatable, dimension(:) :: point
type(plane_t), allocatable, dimension(:) :: xplane, yplane
type(zplane_t), allocatable, dimension(:) :: zplane

type(tavg_t), allocatable, dimension(:,:,:) :: tavg
type(tavg_t), allocatable, dimension(:) :: tavg_zplane

#ifdef PPOUTPUT_EXTRA
type(tavg_sgs_t), allocatable, dimension(:,:,:) :: tavg_sgs
#endif

type(rs_t), allocatable, dimension(:,:,:) :: rs
type(rs_t), allocatable, dimension(:) :: rs_zplane, cnpy_zplane

! add by Mingwei
type(ps_t), allocatable, dimension(:,:,:) :: ps


! Overloaded operators for tavg and rs types
interface operator (.ADD.)
    module procedure tavg_add, tavg_scalar_add, rs_add
end interface

interface operator (.SUB.)
    module procedure tavg_sub, rs_sub
end interface

interface operator (.DIV.)
#ifdef PPOUTPUT_EXTRA
    module procedure tavg_scalar_div, rs_scalar_div, tavg_sgs_scalar_div
#else
    module procedure tavg_scalar_div, rs_scalar_div
#endif
end interface

interface operator (.MUL.)
    module procedure tavg_mul, tavg_scalar_mul
end interface

interface type_set
#ifdef PPOUTPUT_EXTRA
    module procedure tavg_set, rs_set, tavg_sgs_set
#else
    module procedure tavg_set, rs_set
#endif
end interface

interface type_zero_bogus
    module procedure tavg_zero_bogus_2D, tavg_zero_bogus_3D
end interface

contains

#ifdef PPTREES
!*******************************************************************************
function val(this, r, x) result(Rval)
!*******************************************************************************
use functions, only : linear_interp
implicit none
class(tree_ind_func_t), intent(in) :: this
real(rprec), intent(in) :: r, x
real(rprec) :: R1, R23, Rval

R23 = linear_interp(this%r, this%R23, r)
R1 = erf(this%sqrt6overdelta*(this%t_half + x)) +                              &
    erf(this%sqrt6overdelta*(this%t_half - x))
Rval = 0.5 * R1 * R23

end function val

!*******************************************************************************
subroutine init(this, delta2, Rt_all,Lt_all, h_all, N)
!*******************************************************************************
use param, only : write_endian, path, pi
use functions, only : bilinear_interp
implicit none
include'fftw3.f'

class(tree_ind_func_t), intent(inout) :: this
real(rprec), intent(in) :: delta2, Rt_all,Lt_all, h_all
integer, intent(in) :: N

real(rprec) :: L, d, R
integer, dimension(:), allocatable :: ind
real(rprec), dimension(:), allocatable :: yz
real(rprec), dimension(:,:), allocatable :: g, f, h
real(rprec), dimension(:), allocatable :: xi
real(rprec) :: dr, Lr
integer :: i, j

integer*8 plan
complex(rprec), dimension(:,:), allocatable :: ghat, fhat, hhat

L = 4 * Rt_all
d = L / N
R = 0.5 * Rt_all;

allocate(yz(N))
allocate(ind(N))
allocate(g(N, N))
allocate(h(N, N))
allocate(f(N, N))
allocate(ghat(N/2+1, N))
allocate(hhat(N/2+1, N))
allocate(fhat(N/2+1, N))

! Calculate constants
this%t_half = 0.5 * h_all
this%sqrt6overdelta = sqrt(6._rprec) / sqrt(delta2)

! Calculate yz and indices to sort the result
do i = 1, N/2
    yz(i) = d*(i-0.5)
    ind(i) = N/2+i
end do
do i = N/2+1, N
    yz(i) = -L + d*(i-0.5)
    ind(i) = i-N/2
end do

! Calculate g and f
do j = 1, N
    do i = 1, N
        g(i,j) = exp(-6*(yz(i)**2+yz(j)**2)/delta2)
        if (sqrt(yz(i)**2 + yz(j)**2) < R) then
            h(i,j) = 1.0
        else
            h(i,j) = 0.0
        end if
    end do
end do

! Do the convolution f = g*h in fourier space
call dfftw_plan_dft_r2c_2d(plan, N, N, g, ghat, FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan, g, ghat)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_r2c_2d(plan, N, N, h, hhat, FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan, h, hhat)
call dfftw_destroy_plan(plan)

fhat = ghat*hhat

! Compute the inverse fft of fhat
call dfftw_plan_dft_c2r_2d(plan, N, N, fhat, f, FFTW_ESTIMATE)
call dfftw_execute_dft_c2r(plan, fhat, f)
call dfftw_destroy_plan(plan)

! Normalize
f = f / N**2 * d**2

! Sort the results
f = f(ind,ind)
yz = yz(ind);

! Interpolate onto the lookup table
allocate(xi(N))
if (allocated(this%r) ) then
    deallocate(this%r)
end if
allocate( this%r(N) )
allocate( this%R23(N) )

Lr = R + 2 * sqrt(delta2)
dr = Lr / (N - 1)
do i = 1,N
    this%r(i) = (i-1)*dr
    xi(i) = 0
end do
this%R23 = bilinear_interp(yz, yz, f, xi, this%r)
this%R23 = this%R23 / this%R23(1)

end subroutine init
#endif

!///////////////////////////////////////////////////////////////////////////////
!/// TAVG operators
!///////////////////////////////////////////////////////////////////////////////

!*******************************************************************************
function tavg_add( a, b) result(c)
!*******************************************************************************
implicit none
type(tavg_t), intent(in) :: a, b
type(tavg_t) :: c

c % u = a % u + b % u
c % v = a % v + b % v
c % w = a % w + b % w
c % u_w  = a % u_w  + b % u_w
c % v_w  = a % v_w  + b % v_w
c % w_uv = a % w_uv + b % w_uv
c % u2 = a % u2 + b % u2
c % v2 = a % v2 + b % v2
c % w2 = a % w2 + b % w2
c % uv = a % uv + b % uv
c % uw = a % uw + b % uw
c % vw = a % vw + b % vw
!c % dudz = a % dudz + b % dudz
!c % dvdz = a % dvdz + b % dvdz
c % txx = a % txx + b % txx
c % tyy = a % tyy + b % tyy
c % tzz = a % tzz + b % tzz
c % txy = a % txy + b % txy
c % txz = a % txz + b % txz
c % tyz = a % tyz + b % tyz
c % fx = a % fx + b % fx
c % fy = a % fy + b % fy
c % fz = a % fz + b % fz
c % cs_opt2 = a % cs_opt2 + b % cs_opt2

! add by Mingwei
c % p = a % p + b % p
c % p2 = a % p2 + b % p2
c % pu = a % pu + b % pu
c % pv  = a % pv  + b % pv
c % pw  = a % pw  + b % pw

end function tavg_add

!*******************************************************************************
function tavg_sub( a, b) result(c)
!*******************************************************************************
implicit none
type(tavg_t), intent(in) :: a, b
type(tavg_t) :: c

c % u = a % u - b % u
c % v = a % v - b % v
c % w = a % w - b % w
c % u_w  = a % u_w  - b % u_w
c % v_w  = a % v_w  - b % v_w
c % w_uv = a % w_uv - b % w_uv
c % u2 = a % u2 - b % u2
c % v2 = a % v2 - b % v2
c % w2 = a % w2 - b % w2
c % uv = a % uv - b % uv
c % uw = a % uw - b % uw
c % vw = a % vw - b % vw
!c % dudz = a % dudz - b % dudz
!c % dvdz = a % dvdz - b % dvdz
c % txx = a % txx - b % txx
c % tyy = a % tyy - b % tyy
c % tzz = a % tzz - b % tzz
c % txy = a % txy - b % txy
c % txz = a % txz - b % txz
c % tyz = a % tyz - b % tyz
c % fx = a % fx - b % fx
c % fy = a % fy - b % fy
c % fz = a % fz - b % fz
c % cs_opt2 = a % cs_opt2 - b % cs_opt2

!add by Mingwei
c % p = a % p - b % p
c % p2 = a % p2 - b % p2
c % pu = a % pu - b % pu
c % pv  = a % pv  - b % pv
c % pw  = a % pw  - b % pw

end function tavg_sub

!*******************************************************************************
function tavg_scalar_add( a, b ) result(c)
!*******************************************************************************
use types, only : rprec
implicit none

type(tavg_t), intent(in) :: a
real(rprec), intent(in) :: b
type(tavg_t) :: c

c % u = a % u + b
c % v = a % v + b
c % w = a % w + b
c % u_w  = a % u_w  + b
c % v_w  = a % v_w  + b
c % w_uv = a % w_uv + b
c % u2 = a % u2 + b
c % v2 = a % v2 + b
c % w2 = a % w2 + b
c % uv = a % uv + b
c % uw = a % uw + b
c % vw = a % vw + b
!c % dudz = a % dudz + b
!c % dvdz = a % dvdz + b
c % txx = a % txx + b
c % tzz = a % tzz + b
c % tyy = a % tyy + b
c % txy = a % txy + b
c % txz = a % txz + b
c % tyz = a % tyz + b
c % fx = a % fx + b
c % fy = a % fy + b
c % fz = a % fz + b
c % cs_opt2 = a % cs_opt2 + b

!add by Mingwei
c % p = a % p + b
c % p2 = a % p2 + b 
c % pu = a % pu + b 
c % pv  = a % pv  + b 
c % pw  = a % pw  + b

end function tavg_scalar_add

!*******************************************************************************
subroutine tavg_zero_bogus_2D( c )
!*******************************************************************************
use types, only : rprec
implicit none

type(tavg_t), dimension(:,:), intent(inout) :: c

c % txx = 0._rprec
c % tyy = 0._rprec
c % tzz = 0._rprec
c % txy = 0._rprec
c % txz = 0._rprec
c % tyz = 0._rprec
c % fx = 0._rprec
c % fy = 0._rprec
c % fz = 0._rprec

end subroutine tavg_zero_bogus_2D

!*******************************************************************************
subroutine tavg_zero_bogus_3D( c )
!*******************************************************************************
use types, only : rprec
implicit none

type(tavg_t), dimension(:,:,:), intent(inout) :: c

c % txx = 0._rprec
c % tyy = 0._rprec
c % tzz = 0._rprec
c % txy = 0._rprec
c % txz = 0._rprec
c % tyz = 0._rprec
c % fx = 0._rprec
c % fy = 0._rprec
c % fz = 0._rprec

end subroutine tavg_zero_bogus_3D

!*******************************************************************************
function tavg_scalar_div( a, b ) result(c)
!*******************************************************************************
use types, only : rprec
implicit none

type(tavg_t), intent(in) :: a
real(rprec), intent(in) :: b
type(tavg_t) :: c

c % u = a % u / b
c % v = a % v / b
c % w = a % w / b
c % u_w  = a % u_w  / b
c % v_w  = a % v_w  / b
c % w_uv = a % w_uv / b
c % u2 = a % u2 / b
c % v2 = a % v2 / b
c % w2 = a % w2 / b
c % uv = a % uv / b
c % uw = a % uw / b
c % vw = a % vw / b
!c % dudz = a % dudz / b
!c % dvdz = a % dvdz / b
c % txx = a % txx / b
c % tyy = a % tyy / b
c % tzz = a % tzz / b
c % txy = a % txy / b
c % txz = a % txz / b
c % tyz = a % tyz / b
c % fx = a % fx / b
c % fy = a % fy / b
c % fz = a % fz / b
c % cs_opt2 = a % cs_opt2 / b

! add by Mingwei
c % p = a % p / b
c % p2 = a % p2 / b
c % pu = a % pu / b
c % pv  = a % pv  / b
c % pw  = a % pw  / b

end function tavg_scalar_div

!*******************************************************************************
function tavg_mul( a, b) result(c)
!*******************************************************************************
implicit none
type(tavg_t), intent(in) :: a, b
type(tavg_t) :: c

c % u = a % u * b % u
c % v = a % v * b % v
c % w = a % w * b % w
c % u_w  = a % u_w  * b % u_w
c % v_w  = a % v_w  * b % v_w
c % w_uv = a % w_uv * b % w_uv
c % u2 = a % u2 * b % u2
c % v2 = a % v2 * b % v2
c % w2 = a % w2 * b % w2
c % uv = a % uv * b % uv
c % uw = a % uw * b % uw
c % vw = a % vw * b % vw
!c % dudz = a % dudz * b % dudz
!c % dvdz = a % dvdz * b % dvdz
c % txx = a % txx * b % txx
c % tyy = a % tyy * b % tyy
c % tzz = a % tzz * b % tzz
c % txy = a % txy * b % txy
c % txz = a % txz * b % txz
c % tyz = a % tyz * b % tyz
c % fx = a % fx * b % fx
c % fy = a % fy * b % fy
c % fz = a % fz * b % fz
c % cs_opt2 = a % cs_opt2 * b % cs_opt2

! add by Mingwei
c % p = a % p * b % p
c % p2 = a % p2 * b % p2
c % pu = a % pu * b % pu
c % pv  = a % pv * b % pv
c % pw  = a % pw  * b % pw

end function tavg_mul

!*******************************************************************************
function tavg_scalar_mul( a, b ) result(c)
!*******************************************************************************
use types, only : rprec
implicit none

type(tavg_t), intent(in) :: a
real(rprec), intent(in) :: b
type(tavg_t) :: c

c % u = a % u * b
c % v = a % v * b
c % w = a % w * b
c % u_w  = a % u_w  * b
c % v_w  = a % v_w  * b
c % w_uv = a % w_uv * b
c % u2 = a % u2 * b
c % v2 = a % v2 * b
c % w2 = a % w2 * b
c % uv = a % uv * b
c % uw = a % uw * b
c % vw = a % vw * b
!c % dudz = a % dudz * b
!c % dvdz = a % dvdz * b
c % txx = a % txx * b
c % tyy = a % tyy * b
c % tzz = a % tzz * b
c % txy = a % txy * b
c % txz = a % txz * b
c % tyz = a % tyz * b
c % fx = a % fx * b
c % fy = a % fy * b
c % fz = a % fz * b
c % cs_opt2 = a % cs_opt2 * b

! add by Mingwei
c % p = a % p * b
c % p2 = a % p2 * b 
c % pu = a % pu * b 
c % pv  = a % pv * b 
c % pw  = a % pw  * b 

end function tavg_scalar_mul

#ifdef PPOUTPUT_EXTRA
!*******************************************************************************
function tavg_sgs_scalar_div( a, b ) result(c)
!*******************************************************************************
use types, only : rprec
implicit none

type(tavg_sgs_t), intent(in) :: a
real(rprec), intent(in) :: b
type(tavg_sgs_t) :: c

!c % Tn = a % Tn / b
c % Nu_t = a % Nu_t / b
!c % F_LM = a % F_LM / b
!c % F_MM = a % F_MM / b
!c % F_QN = a % F_QN / b
!c % F_NN = a % F_NN / b
!c % ee_now = a % ee_now / b
!#ifdef PPDYN_TN
!c % F_ee2 = a % F_ee2 / b
!c % F_deedt2 = a % F_deedt2 / b
!#endif

end function tavg_sgs_scalar_div
#endif

!*******************************************************************************
function tavg_interp_to_uv_grid( a ) result(c)
!*******************************************************************************
use param, only: lbz
use functions, only : interp_to_uv_grid
implicit none

type(tavg_t), dimension(:,:,lbz:), intent(in) :: a
type(tavg_t), allocatable, dimension(:,:,:) :: c

integer :: ubx, uby, ubz

ubx = ubound(a,1)
uby = ubound(a,2)
ubz = ubound(a,3)

allocate(c(ubx,uby,lbz:ubz))

c = a

c % fz = interp_to_uv_grid(a % fz, lbz )
c % w  = interp_to_uv_grid(a % w,lbz)
c % w2 = interp_to_uv_grid(a % w2,lbz)

end function tavg_interp_to_uv_grid

!*******************************************************************************
function tavg_interp_to_w_grid( a ) result(c)
!*******************************************************************************
use param, only: lbz
use functions, only : interp_to_w_grid
implicit none

type(tavg_t), dimension(:,:,lbz:), intent(in) :: a
type(tavg_t), allocatable, dimension(:,:,:) :: c

integer :: ubx, uby, ubz

ubx = ubound(a,1)
uby = ubound(a,2)
ubz = ubound(a,3)

allocate(c(ubx,uby,lbz:ubz))

c = a

c % txx =  interp_to_w_grid( a % txx, lbz )
c % tyy =  interp_to_w_grid( a % tyy, lbz )
c % tzz =  interp_to_w_grid( a % tzz, lbz )
c % txy =  interp_to_w_grid( a % txy, lbz )

c % fx = interp_to_w_grid( a % fx, lbz )
c % fy = interp_to_w_grid( a % fy, lbz )

end function tavg_interp_to_w_grid

!///////////////////////////////////////////////////////////////////////////////
!/// RS operators
!///////////////////////////////////////////////////////////////////////////////

!*******************************************************************************
function rs_add( a, b) result(c)
!*******************************************************************************
implicit none

type(rs_t), intent(in) :: a, b
type(rs_t) :: c

c % up2 = a % up2 + b % up2
c % vp2 = a % vp2 + b % vp2
c % wp2 = a % wp2 + b % wp2
c % upvp = a % upvp + b % upvp
c % upwp = a % upwp + b % upwp
c % vpwp = a % vpwp + b % vpwp

end function rs_add

!*******************************************************************************
function rs_sub( a, b) result(c)
!*******************************************************************************
implicit none

type(rs_t), intent(in) :: a, b
type(rs_t) :: c

c % up2 = a % up2 - b % up2
c % vp2 = a % vp2 - b % vp2
c % wp2 = a % wp2 - b % wp2
c % upvp = a % upvp - b % upvp
c % upwp = a % upwp - b % upwp
c % vpwp = a % vpwp - b % vpwp

end function rs_sub

!*******************************************************************************
function rs_scalar_div( a, b) result(c)
!*******************************************************************************
implicit none

type(rs_t), intent(in) :: a
real(rprec), intent(in) :: b
type(rs_t) :: c

c % up2 = a % up2 / b
c % vp2 = a % vp2 / b
c % wp2 = a % wp2 / b
c % upvp = a % upvp / b
c % upwp = a % upwp / b
c % vpwp = a % vpwp / b

end function rs_scalar_div

!///////////////////////////////////////////////////////////////////////////////
!/// Spectral RS operators
!///////////////////////////////////////////////////////////////////////////////

!*******************************************************************************
function rs_compute( a , lbz2) result(c)
!*******************************************************************************
implicit none
integer, intent(in) :: lbz2
type(tavg_t), dimension(:,:,lbz2:), intent(in) :: a
type(rs_t), allocatable, dimension(:,:,:) :: c

integer :: ubx, uby, ubz

ubx=ubound(a,1)
uby=ubound(a,2)
ubz=ubound(a,3)

allocate(c(ubx,uby,lbz2:ubz))

c % up2 = a % u2 - a % u * a % u
c % vp2 = a % v2 - a % v * a % v
c % wp2 = a % w2 - a % w * a % w
c % upvp = a % uv - a % u * a % v
!! using u_w and v_w below instead of u and v ensures that the Reynolds
!! stresses are on the same grid as the squared velocities (i.e., w-grid)
c % upwp = a % uw - a % u_w * a % w   !!pj
c % vpwp = a % vw - a % v_w * a % w   !!pj

end function rs_compute

!///////////////////////////////////////////////////////////////////////////////
!/// Spectral PS operators
!///////////////////////////////////////////////////////////////////////////////
!add by Mingwei for pressure term
!*******************************************************************************
function ps_compute( a , lbz2) result(c)
!*******************************************************************************
implicit none
integer, intent(in) :: lbz2
type(tavg_t), dimension(:,:,lbz2:), intent(in) :: a
type(ps_t), allocatable, dimension(:,:,:) :: c

integer :: ubx, uby, ubz

ubx=ubound(a,1)
uby=ubound(a,2)
ubz=ubound(a,3)

allocate(c(ubx,uby,lbz2:ubz))

c % pp2 = a % p2 - a % p * a % p
c % ppup = a % pu - a % p * a % u
c % ppvp = a % pv - a % p * a % v
c % ppwp = a % pw - a % p * a % w_uv  !the average value w_uv_f is the same

end function ps_compute

!*******************************************************************************
function cnpy_tavg_mul( a ) result(c)
!*******************************************************************************
!
! This performs one set of multiplication for the canopy stresses
!
implicit none

type(tavg_t), intent(in) :: a
type(rs_t) :: c

c % up2 = a % u * a % u
c % vp2 = a % v * a % v
c % wp2 = a % w * a % w
c % upvp = a % u * a % v
c % upwp = a % u * a % w
c % vpwp = a % v * a % w

end function cnpy_tavg_mul

!///////////////////////////////////////////////////////////////////////////////
!/// Spectral TAVG operators
!///////////////////////////////////////////////////////////////////////////////

!*******************************************************************************
subroutine tavg_set( c, a )  !not include p related term
!*******************************************************************************
use types, only : rprec
implicit none
real(rprec), intent(in) :: a
type(tavg_t), intent(out) :: c

c % u = a
c % v = a
c % w = a
c % u_w  = a
c % v_w  = a
c % w_uv = a
c % u2 = a
c % v2 = a
c % w2 = a
c % uv = a
c % uw = a
c % vw = a
!c % dudz = a
!c % dvdz = a
c % txx = a
c % tyy = a
c % tzz = a
c % txy = a
c % txz = a
c % tyz = a
c % fx = a
c % fy = a
c % fz = a
c % cs_opt2 = a

!add by Mingwei for type set
c % p = a
c % p2 = a
c % pu = a
c % pv  = a
c % pw  = a

end subroutine tavg_set

#ifdef PPOUTPUT_EXTRA
!*******************************************************************************
subroutine tavg_sgs_set( c, a )
!*******************************************************************************
use types, only : rprec
implicit none
real(rprec), intent(in) :: a
type(tavg_sgs_t), intent(out) :: c

!c % Tn =  a
c % Nu_t =  a
!c % F_LM =  a
!c % F_MM =  a
!c % F_QN =  a
!c % F_NN =  a
!c % ee_now = a
!#ifdef PPDYN_TN
!c % F_ee2 = a
!c % F_deedt2 = a
!#endif

end subroutine tavg_sgs_set
#endif

!///////////////////////////////////////////////////////////////////////////////
!/// Spectral RS subroutines
!///////////////////////////////////////////////////////////////////////////////

!*******************************************************************************
subroutine rs_set( c, a )
!*******************************************************************************
use types, only : rprec
implicit none
real(rprec), intent(in) :: a
type(rs_t), intent(out) :: c

c % up2 = a
c % vp2 = a
c % wp2 = a
c % upvp = a
c % upwp = a
c % vpwp = a

end subroutine rs_set

end module stat_treesdefs

