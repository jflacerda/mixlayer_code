program mixing-layer

! Program to generate baseflow of a 2D compressible mixing layer
! This code is a extension of code from Leonardo da Costa Salemi (LCP/INPE)
! and Márcio Teixeira de Mendonça (CTA/IAE)

!     This code solves the compressible binary shear layer equations transformed
!    according to Lees-Dorodnitsyn Transformation (Anderson, Hypersonic 
!     and High Temperature Gas Dynamics).
! 
!     The code solves the equations for Chi,Pr,Le =1 and Chi,Pr,Le = constant.
!
!     Velocity, Enthalpy and Mass Fraction profiles: 
!
!     f' - velocity - u/u1     
!     g  - enthalpy - He/He1

! Jonatas F Lacerda    jonatasflacerda.hotmail.com
! ICMC - USP

! date: 19.03.2013

implicit none
include '../var.f90'
include 'comm.mshbase'
include 'comm.varbase'

 ! Mesh generation
 call mesh_baseflow
 
 ! Read input data from coupled.in
 open(unit=1,file='coupled.in',status='UNKNOWN')
  read(1,*)
  read(1,*)
  read(1,*)eps
  read(1,*)h1,hmin
  read(1,*)
  read(1,*)
  read(1,*)L,R
  read(1,*)
  read(1,*)
  read(1,*)P,Q
  read(1,*)
  read(1,*)
  read(1,*)W,Z
  read(1,*)
  read(1,*)
  read(1,*)press
  read(1,*)T(kmax),u(kmax)      ! Fast side data
  read(1,*)T(1),u(1)            ! Slow side data
  read(1,*)
  read(1,*)
  read(1,*)ngas1
  read(1,*)
  read(1,*)
  read(1,*)ngas2
  read(1,*)
  read(1,*)
  read(1,*)const
  read(1,*)
  read(1,*)
  read(1,*)xpos
 close(1)
	  
 stop

end program mixing-layer

!!###############################################################################

subroutine mesh_baseflow

! generates non-uniform mesh in x and y direction for the mixing layer
! in y direction is created the upper part of the mixing layer grid
! and this grid is copied to the lower side and later both parts
! are united into 'y' variable

implicit none
include '../var.f90'
include 'comm.mshbase'

integer(4):: icount, icountmax
integer(4):: j, jmid, jmaxver, jmaxaux, i, istf1, istf2, ind
real(8):: a3, eta0, eta, stxf, ep, posx, dx
real(8),allocatable,dimension(:):: yaux1, yaux2
logical:: lstop, lwrite

! logical variable to print mesh
lwrite = .false.

!******************************************************
  ! 1) y direction calculation
  lstop = .false.
  jmaxaux =  jmax
  icountmax = 40
  icount = 1
  write(*,*) 'jmax =', jmax
  do while (.not.lstop)
    ! calculating middle value of y-index: jmid
    jmid = (jmaxaux+1)/2
    
    ! calculating index verification
    jmaxver = (2*jmid)-1
    
    !verifying if jmid is correct
    if (jmaxver.ne.jmaxaux) then
      jmaxaux = jmaxaux + 1
      icount = icount + 1
	  if (icount.gt.icountmax) then
        write(*,*) 'not worked in', icountmax,' iterations'
        lstop=.true.
      endif
    else
      allocate (yaux1(jmid))
      allocate (yaux2(jmaxaux))	  
      lstop=.true.
    endif
  enddo
  
  !coeficients for mesh calculation (According to Babucke, 2009)
  a3 = 0.000502236538416527d0
  eta0 =  0.d0

  ! calculating the upper part yaux1
  eta = 0.d0
  yaux1(1) = eta
  do j = 2, jmid
    yaux1(j) = a3 * (eta - eta0)**3 + (eta - eta0)
	eta = eta + dn
  enddo
  
  ! copying yaux1 to the lower part yaux2
  yaux2(1) = yaux1(jmid)
  do j = 1,jmid-1
    yaux2(j+1) = yaux1(jmid-j)
  enddo
  
  ! rewrites yaux1 on the upper side of yaux2
  do j = jmid+1, jmaxaux
    yaux2(j) = yaux1(j-jmid+1)
  enddo
  
  ! calculates y on the lower side
  do j = 1, jmid
    y(j) = -yaux2(j)
  enddo

  ! calculates y on the upper side
  do j = jmid+1, jmax
    y(j) = yaux2(j)
  enddo
!******************************************************
  ! 2) x direction calculation

  !coeficients for mesh calculation (According to Babucke, 2009)
  stx = 1.d0
  stxf = 1.009d0
  istf1 = 2100
  istf2 = istf1 + 30
  dx = 0.157d0

  posx = 30.d0
  do i = 1,imax
    ! assign x position  
    x(i) = posx
	
	! identify the 'i' index where the mesh in x direction transitions
	! from regular to stretched mesh
    if (i.ge.istf1 .and. i.le.istf2) then
      ep  = dble(i-istf1)/dble(istf2-istf1)
      stx = (stxf-1.d0)*((6.d0*ep-15.d0)*ep+10.d0)*ep**3+1.d0
    else
      if (i.gt.istf2) stx = stxf  
    endif

    posx = posx+dx
    dx = dx*stx
  enddo

!******************************************************
  if(lwrite) then
    !write a mesh file
    open(unit=1, file='meshbase.out',status='unknown')
     ind = 1
     write(1,*) 'variables = "x","y"'
     write(1,*) 'zone i=',imax/ind,', j=',jmax/ind
     do j = 1, jmax, ind
       do i = 1, imax, ind
         write(1,10) x(i),y(j)
       enddo
     enddo
     10 format(1x,2(e14.7,x))
    close (1)
  endif

return
end subroutine mesh_baseflow

!!###############################################################################
! subroutine escreve(itime)
! implicit none
! include 'var.f90'
! include 'comm.par'
! include 'comm.var'
! include 'comm.msh'
! 
! character(len=18)      :: nome
! integer(4), intent(in) :: itime
! integer(4):: i, j
! real(8),dimension(ptsx,ptsy)     :: rho, rhou, rhov, ezaot, u, v, et
!
! equivalence (q(1,1,1),rho(1,1)), (q(1,1,2),rhou(1,1)), (q(1,1,3),rhov(1,1)), (q(1,1,4),ezaot(1,1))
! 
! write(nome,'(a,i0.2,a)')'saida_',myrank,'_.bin'
! 
! if (lprintform) then
!   u  = rhou  / rho
!   v  = rhov  / rho
!   et = ezaot / rho
!   open (1,file=nome,status='unknown')
!     write(1,*) 'variables = "x","y","rho","u","v","et","erro1","erro2","erro3","erro4"'
!     write(1,*) 'zone i=',ptsx,', j=',ptsy
!     10 format(x,10(e24.10,x))
!     do j = 1, ptsy
!       do i = 1,ptsx
!         write(1,10) x(i), y(j), rho(i,j), u(i,j), v(i,j), et(i,j), &
!                     erro(i,j,1), erro(i,j,2), erro(i,j,3), erro(i,j,4)
!       enddo
!     enddo
!   close (1)
! else
!   open(1,file=nome,form='unformatted')
!    write(1) itime
!    write(1) q, f, g
!   close (unit=1)
! endif
! 
! return
! end subroutine
!!###############################################################################
