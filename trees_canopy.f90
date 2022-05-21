!!  added by Xiantao Fan
!!  Copyright (C) 2009-2013  Johns Hopkins University
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
!**********************************************************************
module trees_canopy

use types, only : rprec
use param
use grid_m
use messages 
use string_util
use stat_treesdefs, only : tree_farm

#ifdef PPMPI
use mpi_defs, only : MPI_SYNC_DOWNUP, mpi_sync_real_array 
#endif

implicit none

save
private

!define parameters of trees: Cd,pi,am_all,zm,Rt_all,Lt_all, h_all,n_all(input);a(z)(calculate)
public :: trees_init, trees_forcing,  trees_finalize
public :: trees_checkpoint  
character (*), parameter :: mod_name = 'trees_canopy'

! The following values are read from the input file
! number of trees in the x-direction
integer, public :: num_xx 
! number of trees in the y-direction
integer, public :: num_yy            
! baseline drag coefficient
real(rprec), public :: Cd      
! tree height in meters
real(rprec), public :: h_all
! height in meters of maximum tress leaf-aera-density
real(rprec), public :: zm_all  
! maximum tress leaf-aera-density in meters
real(rprec), public :: am_all      

! the width of trees
real(rprec), public :: Rt_all  
real(rprec), public :: Lt_all

! effective coefficient of LAD in x-direction
real(rprec), public :: p_x
! effective coefficient of LAD in y-direction
real(rprec), public :: p_y
! effective coefficient of LAD in z-direction
real(rprec), public :: p_z

! orientation of trees
integer, public :: orientation1      
! stagger percentage from baseline
real(rprec), public :: stag_perc1 

! Read parameters from input_trees/param.dat
logical, public :: read_param1        

! disk-avg time scale in seconds (default 600) 
real(rprec), public :: T_avg_dim1
! filter size as multiple of grid spacing
real(rprec), public :: alpha1
! indicator function only includes values above this threshold
real(rprec), public :: filter_cutoff1   
! Number of timesteps between the output
integer, public :: tbase1            

! The following are derived from the values above
integer :: nloc             ! total number of trees
real(rprec) :: sx           ! spacing in the x-direction, multiple of diameter
real(rprec) :: sy           ! spacing in the y-direction



! Input files
character(*), parameter :: input_folder = 'input_trees/'
character(*), parameter :: param_dat = path // input_folder // 'param.dat'

! Output files

character(*), parameter :: output_folder = 'tree/'
character(*), parameter :: vel_top_dat = path // output_folder // 'vel_top.dat'
character(*), parameter :: u_d_T_dat = path // output_folder // 'u_d_T.dat'
character(*), parameter :: node_dat = path // output_folder //'node_dat.dat'
integer, dimension(:), allocatable :: forcing_fid    



! epsilon used for velocity time-averaging
real(rprec) :: eps 

! Commonly used indices
integer :: i,j,k,i2,j2,k2,l,s
integer :: k_start, k_end


contains
    
    !*******************************************************************************
subroutine trees_init()
!*******************************************************************************
!
! This subroutine creates the 'tree' folder and starts the tree forcing 
! output files. It also creates the indicator function (Gaussian-filtered from 
! binary locations - in or out) and sets values for tree 
! (node locations, etc)
!
use open_file_fid_mod
implicit none

real(rprec), pointer, dimension(:) :: x,y,z
character (*), parameter :: sub_name = mod_name // '.trees_init'
integer :: fid
real(rprec) :: T_avg_dim_file
logical ::  exst
character (100) :: string1


! Allocate and initialize
nloc = num_xx*num_yy
nullify(tree_farm%tree)
allocate(tree_farm%tree(nloc))
allocate(forcing_fid(nloc))


! Create tree directory
call system("mkdir -vp tree") 


! Spacing between trees (as multiple of mean tree crown size)
sx = L_x / (num_xx * Rt_all )
sy = L_y / (num_yy * Lt_all )

! Place the trees and specify some parameters
call place_trees

! Resize thickness to capture at least on plane of gridpoints
! and set baseline values for size
do k = 1, nloc
    tree_farm%tree(k)%h_all = max(tree_farm%tree(k)%h_all,dz*1.01)
    tree_farm%tree(k)%Rt_all = max(tree_farm%tree(k)%Rt_all,dx*1.01)
    tree_farm%tree(k)%Lt_all = max(tree_farm%tree(k)%Lt_all,dy*1.01)
end do

! Specify starting and ending indices for the processor
#ifdef PPMPI
k_start = 1+coord*(nz-1)
k_end = nz-1+coord*(nz-1)
#else
k_start = 1
k_end = nz
#endif


! Find tree nodes - including filtered ind, n_hat, num_nodes, and nodes for 
! each tree. Each processor finds trees in its domain
call trees_nodes

! Read the time-averaged  velocities from file if available
if (coord == 0) then
    inquire (file=u_d_T_dat, exist=exst)
    if (exst) then
        write(*,*) 'Reading from file ', trim(u_d_T_dat)
        fid = open_file_fid( u_d_T_dat, 'rewind', 'formatted' )
        do i=1,nloc
            read(fid,*) tree_farm%tree(i)%u_d_T    
        end do    
        read(fid,*) T_avg_dim_file
        if (T_avg_dim_file /= T_avg_dim1) then
            write(*,*) 'Time-averaging window does not match value in ',   &
                       trim(u_d_T_dat)
        end if
        close (fid)
    else  
        write (*, *) 'File ', trim(u_d_T_dat), ' not found'
        write (*, *) 'Assuming u_d_T = -1. for all trees'
        do k=1,nloc
            tree_farm%tree(k)%u_d_T = -1.
        end do
    end if                                    
end if

! Generate top of domain file
if (coord .eq. nproc-1) then
    fid = open_file_fid( vel_top_dat, 'rewind', 'formatted' )
    close(fid)
end if

! Generate the files for the tree forcing output
if(coord==0) then
    do s=1,nloc
        call string_splice( string1, path // 'tree/tree_', s, '.dat' )
        forcing_fid(s) = open_file_fid( string1, 'append', 'formatted' )
    end do
end if


end subroutine trees_init

!*******************************************************************************
subroutine trees_nodes
!*******************************************************************************
!
! This subroutine locates nodes for each tree and builds the arrays: ind, 
! n_hat, num_nodes, and nodes
!
implicit none

character (*), parameter :: sub_name = mod_name // '.trees_nodes'

real(rprec), pointer :: p_xloc => null(), p_yloc => null(), p_height => null()
real(rprec), pointer :: p_Rt => null(), p_Lt => null()

integer :: icp, jcp, kcp
integer :: imax, jmax, kmax, k3
integer :: min_i, max_i, min_j, max_j, min_k, max_k
integer :: count_i, count_n

real(rprec) :: filt
real(rprec), dimension(:), allocatable :: sumA

allocate(sumA(nloc))
sumA = 0

! z_tot for total domain (since z is local to the processor)

do s=1,nloc
    
    count_n = 0    !used for counting nodes for each tree
    count_i = 1    !index count - used for writing to array "nodes"

    p_xloc => tree_farm%tree(s)%xloc     
    p_yloc => tree_farm%tree(s)%yloc  
    p_height => tree_farm%tree(s)%h_all 
    p_Rt => tree_farm%tree(s)%Rt_all
    p_Lt => tree_farm%tree(s)%Lt_all

    !identify "search area"
    imax = int(0.5*p_Rt/dx +0.5)
    jmax = int(0.5*p_Lt/dy )
    kmax = int((0.5*p_height/dz))



    !determine nearest (i,j,k) to tree center
    icp = nint(p_xloc/dx)
    jcp = nint(p_yloc/dy)
    kcp = nint(0.5*p_height/dz )

    !determine limits for checking i,j,k
    !due to spectral BCs, i and j may be < 1 or > nx,ny
    !the mod function accounts for this when these values are used
    min_i = icp-imax
    max_i = icp+imax
    min_j = jcp-jmax
    max_j = jcp+jmax
    min_k = max((kcp-kmax),1)
    max_k = min((kcp+kmax),nz_tot)
    tree_farm%tree(s)%nodes_max(1) = min_i
    tree_farm%tree(s)%nodes_max(2) = max_i
    tree_farm%tree(s)%nodes_max(3) = min_j
    tree_farm%tree(s)%nodes_max(4) = max_j
    tree_farm%tree(s)%nodes_max(5) = min_k
    tree_farm%tree(s)%nodes_max(6) = max_k      

    ! check neighboring grid points
    ! update num_nodes, nodes, and ind for this tree
    ! split domain between processors
    ! z(nz) and z(1) of neighboring coords match so each coord gets 
    ! (local) 1 to nz-1
    tree_farm%tree(s)%ind = 0._rprec
    tree_farm%tree(s)%nodes = 0
    tree_farm%tree(s)%num_nodes = 0
    count_n = 0
    count_i = 1
    
    do k=k_start,k_end  !global k     
        do j=min_j,max_j
            do i=min_i,max_i
             if (i<1) then
                    i2 = mod(i+nx-1,nx)+1
                    
                elseif (i>nx) then
                    i2 = mod(i+nx-1,nx)+1
                    
                else
                    i2 = i
                    
                end if            
                if (j<1) then
                    j2 = mod(j+ny-1,ny)+1
                                   
                elseif (j>ny) then
                    j2 = mod(j+ny-1,ny)+1
                    
                else
                    j2 = j
                    
                end if   
               ! get the filter value
                filt = 1 !tree_farm%tree(s)%tree_ind_func%val(r_Rt_all, r_norm)
                if ( filt > filter_cutoff1 ) then
                    tree_farm%tree(s)%nodes(count_i,1) = i2
                    tree_farm%tree(s)%nodes(count_i,2) = j2
                    tree_farm%tree(s)%nodes(count_i,3) = k-coord*(nz-1)!local
                    count_n = count_n + 1
                    count_i = count_i + 1
                    sumA(s) = sumA(s) + filt * dx * dy * dz
                    k2=k-coord*(nz-1)
                    k3=k                                       
                end if
           end do
       end do
    end do
    tree_farm%tree(s)%num_nodes = count_n
    
   
end do

 open(unit=1,file=node_dat,status='unknown',form='formatted',            &
    action='write',position='append')
    write(1,*) i2, j2, k2, k3
    close(1) 
       
! Cleanup
deallocate(sumA)

end subroutine trees_nodes


!********************************************************************************
subroutine trees_forcing()
!*******************************************************************************
! 
! This subroutine applies the drag forcing by trees
! 
use sim_param, only : u,v,w, fxa_t,fya_t,fza_t
use functions, only : interp_to_uv_grid
implicit none

character (*), parameter :: sub_name = mod_name // '.trees_forcing'



real(rprec) :: ind2, velocity

real(rprec), pointer :: p_Cd => null(), p_px=> null(), p_py=> null(), p_pz=> null()
real(rprec), pointer :: p_am=> null(),p_zm=> null(),p_h=> null()
real(rprec), allocatable, dimension(:,:,:) :: w_uv

integer :: kstart, k3


#ifdef PPMPI
!syncing intermediate w-velocities
call mpi_sync_real_array(w, 0, MPI_SYNC_DOWNUP)     
#endif


!transfer the velocity in x,y,z to MPI.
w_uv = interp_to_uv_grid(w, lbz)

   
!apply forcing to each node
  do s=1,nloc
        p_Cd => tree_farm%tree(s)%Cd
        p_px => tree_farm%tree(s)%p_x
        p_py => tree_farm%tree(s)%p_y
        p_pz => tree_farm%tree(s)%p_z
        p_am => tree_farm%tree(s)%am_all
        p_zm => tree_farm%tree(s)%zm_all
        p_h => tree_farm%tree(s)%h_all       
   do l=1,tree_farm%tree(s)%num_nodes
      i2 = tree_farm%tree(s)%nodes(l,1)
      j2 = tree_farm%tree(s)%nodes(l,2)
      k2 = tree_farm%tree(s)%nodes(l,3)
      
! Make k2 the global index 
#ifdef PPMPI
k3=k2 + coord*(nz-1)
#else
k3=k2
#endif

!================test the node
! set the point of the velocity sqrt(u(i2,j2,k3)**2 + v(i2,j2,k3)**2 + w(i2,j2,k3)**2)
velocity=sqrt(u(i2,j2,k2)**2 + v(i2,j2,k2)**2 + w_uv(i2,j2,k2)**2)


if (k3*dz .LE. p_zm .and. k3 .GE. 0) then 
       fxa_t(i2,j2,k2) = -p_Cd*abs(velocity)*u(i2,j2,k2)*p_am*(((p_h-p_zm)/(p_h-k3*dz))**6) *exp(6*(1-(p_h-p_zm)/(p_h-k3*dz)))*p_px
       fya_t(i2,j2,k2) = -p_Cd*abs(velocity)*v(i2,j2,k2)*p_am*(((p_h-p_zm)/(p_h-k3*dz))**6) *exp(6*(1-(p_h-p_zm)/(p_h-k3*dz)))*p_py
       fza_t(i2,j2,k2)= -p_Cd*abs(velocity)*w_uv(i2,j2,k2)*p_am*(((p_h-p_zm)/(p_h-k3*dz))**6) *exp(6*(1-(p_h-p_zm)/(p_h-k3*dz)))*p_pz
else if (k3*dz .GT. p_zm .and. k3*dz .LT. p_h) then
       fxa_t(i2,j2,k2) = -p_Cd*abs(velocity)*u(i2,j2,k2)*p_am*(((p_h-p_zm)/(p_h-k3*dz))**0.5) *exp(0.5*(1-(p_h-p_zm)/(p_h-k3*dz)))*p_px
       fya_t(i2,j2,k2)= -p_Cd*abs(velocity)*v(i2,j2,k2)*p_am*(((p_h-p_zm)/(p_h-k3*dz))**0.5) *exp(0.5*(1-(p_h-p_zm)/(p_h-k3*dz)))*p_py
       fza_t(i2,j2,k2)= -p_Cd*abs(velocity)*w_uv(i2,j2,k2)*p_am*(((p_h-p_zm)/(p_h-k3*dz))**0.5) *exp(0.5*(1-(p_h-p_zm)/(p_h-k3*dz)))*p_pz
else if (k3*dz .GE. p_h) then
       fxa_t(i2,j2,k2) = 0._rprec
       fya_t(i2,j2,k2) = 0._rprec
       fza_t(i2,j2,k2)= 0._rprec
end if

end do
end do



!spatially average velocity at the top of the domain and write to file
if (coord .eq. nproc-1) then
    open(unit=1,file=vel_top_dat,status='unknown',form='formatted',            &
    action='write',position='append')
    write(1,*) total_time, sum(u(:,:,nz-1))/(nx*ny)
    close(1)
end if

! Cleanup


end subroutine trees_forcing


!*******************************************************************************
subroutine trees_finalize ()
!*******************************************************************************
implicit none

character (*), parameter :: sub_name = mod_name // '.trees_finalize'

!write tree-averaged velocity to file along with T_avg_dim
!useful if simulation has multiple runs   >> may not make a large difference
call trees_checkpoint ()
    
!deallocate
deallocate(tree_farm%tree) 
    
end subroutine trees_finalize

!*******************************************************************************
subroutine trees_checkpoint ()
!*******************************************************************************
! 
!
!
use open_file_fid_mod
implicit none

character (*), parameter :: sub_name = mod_name // '.trees_checkpoint'
integer :: fid

!write disk-averaged velocity to file along with T_avg_dim
!useful if simulation has multiple runs   >> may not make a large difference
if (coord == 0) then  
    fid = open_file_fid( u_d_T_dat, 'rewind', 'formatted' )
    do i=1,nloc
        write(fid,*) tree_farm%tree(i)%u_d_T
    end do           
    write(fid,*) T_avg_dim1
    close (fid)
end if
    
end subroutine trees_checkpoint




!********************************************************************************
subroutine place_trees
!*******************************************************************************
!
! This subroutine places the trees on the domain. It also sets the values for
! each individual tree. After the subroutine is called, the following values 
! are set for each tree : xloc, yloc, drag coefficient cd, height h_all, heighet of LAD zm_all, maximum LAD am_all,  crow size of trees Rt_all and Lt_all
! efecitive coefficient p_x, p_y, p_z.
!
use param, only: pi, z_i
use open_file_fid_mod
use messages
implicit none

character(*), parameter :: sub_name = mod_name // '.place_trees'

logical :: exst
integer :: fid

! Read parameters from file if needed
if (read_param1) then
    ! Check if file exists and open
    inquire (file = param_dat, exist = exst)
    if (.not. exst) then
        call error (sub_name, 'file ' // param_dat // 'does not exist')
    end if

    ! Check that there are enough lines from which to read data
    nloc = count_lines(param_dat)
    if (nloc < num_xx*num_yy) then
        nloc = num_xx*num_yy
        call error(sub_name, param_dat // 'must have num_xx*num_yy lines')
    else if (nloc > num_xx*num_yy) then
        call warn(sub_name, param_dat // ' has more than num_xx*num_yy lines. '  &
                  // 'Only reading first num_xx*num_yy lines')
    end if
    
    ! Read from parameters file, which should be in this format:
    !xloc [meters], yloc[meters], drag coefficient cd, height h_all[meters], heighet of LAD zm_all[meters], maximum LAD am_all,  width of trees Rt_all and Lt_all
    ! efecitive coefficient p_x, p_y, p_z. 
  
    write(*,*) "Reading from", param_dat
    fid = open_file_fid(param_dat, 'rewind', 'formatted')
    do k = 1, nloc
        read(fid,*) tree_farm%tree(k)%xloc, tree_farm%tree(k)%yloc,      &
            tree_farm%tree(k)%Cd, tree_farm%tree(k)%h_all,               &
            tree_farm%tree(k)%zm_all, tree_farm%tree(k)%am_all,          &
            tree_farm%tree(k)%Rt_all, tree_farm%tree(k)%Lt_all, tree_farm%tree(k)%p_x,             &
            tree_farm%tree(k)%p_y, tree_farm%tree(k)%p_z 
    end do
    close(fid)
    

end if

end subroutine place_trees



!********************************************************************************
function count_lines(fname) result(N)
!*******************************************************************************
!
! This function counts the number of lines in a file
!
use open_file_fid_mod
use messages
use param, only : CHAR_BUFF_LENGTH
implicit none
character(*), intent(in) :: fname
logical :: exst
integer :: fid, ios
integer :: N

character(*), parameter :: sub_name = mod_name // '.count_lines'

! Check if file exists and open
inquire (file = trim(fname), exist = exst)
if (.not. exst) then
    call error (sub_name, 'file ' // trim(fname) // 'does not exist')
end if
fid = open_file_fid(trim(fname), 'rewind', 'formatted')

! count number of lines and close
ios = 0
N = 0
do 
    read(fid, *, IOstat = ios)
    if (ios /= 0) exit
    N = N + 1
end do

! Close file
close(fid)

end function count_lines

end module trees_canopy


