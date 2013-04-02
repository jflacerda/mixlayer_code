program baseflow

! Program to generate baseflow of a 2D compressible mixing layer
! This code is a extension of code from Leonardo da Costa Salemi (LCP/INPE)
! and Márcio Teixeira de Mendonça (CTA/IAE)
!
!     This code soluves the compressible binary shear layer equations transformed
!    according to Lees-Dorodnitsyn Transformation (Anderson, Hypersonic 
!     and High Temperature Gas Dynamics).
! 
!     The code soluves the equations for Chi,Pr,Le =1 and Chi,Pr,Le = constant.
!
!     Velocity, Enthalpy and Mass Fraction profiles: 
!
!     f' - velocity - u/u1     
!     g  - enthalpy - He/He1

! Jonatas F Lacerda    jonatasflacerda.hotmail.com
! ICMC - USP

! date: 04.02.2013

implicit none
include '../var.f90'
include 'comm.mshbase'
include 'comm.varbase'

 !........................ Mesh generation
 call mesh_baseflow
 
 !........................ Calculation of eta variable
 call calc_eta
 
 !........................ Integration of equations
 call integra
 

 stop

end program baseflow

!*******************************************************************************

subroutine mesh_baseflow

! generates non-uniform mesh in x and y direction for the mixing layer
! in y direction is created the upper part of the mixing layer grid
! and this grid is copied to the lower side and later both parts
! are united into 'y' variable

implicit none
include '../var.f90'
include 'comm.mshbase'

integer(4):: icount, icountmax
!integer(4):: j, jmid, jmaxver, jmaxaux, i, istf1, istf2, ind
integer(4):: j, jmaxver, i, istf1, istf2, ind
real(8):: a3, ksi0, ksi, stxf, ep, posx, dx, dy
real(8),allocatable,dimension(:):: yaux1, yaux2
logical:: lstop, lwrite

! logical variable to print mesh
lwrite = .false.

!******************************************************
  !1) y direction calculation - uniform mesh
  !dy = 0.15d0
  !
  !lstop = .false.
  !jmaxaux =  jmax
  !icountmax = 40
  !icount = 1
  !do while (.not.lstop)
  !  ! calculating middle value of y-index: jmid
  !  jmid = (jmaxaux+1)/2
  !  
  !  ! calculating index verification
  !  jmaxver = (2*jmid)-1
  !  
  !  !verifying if jmid is correct
  !  if (jmaxver.ne.jmaxaux) then
  !    jmaxaux = jmaxaux + 1
  !    icount = icount + 1
  !    if (icount.gt.icountmax) then
  !      write(*,*) 'not worked in', icountmax,' iterations'
  !      lstop=.true.
  !    endif
  !  else
  !    write(*,*) 'jmax = ', jmax, 'jmaxaux= ', jmaxaux
  !    allocate (yaux1(jmid))
  !    allocate (yaux2(jmaxaux))	  
  !    lstop=.true.
  !  endif
  !enddo
  !
  !! calculating the upper part yaux1
  !ksi = 0.d0
  !do j = 1, jmid
  !  yaux1(j) = ksi
  !  ksi = ksi + dy
  !enddo
  !
  !! copying yaux1 to the lower part yaux2
  !yaux2(1) = yaux1(jmid)
  !do j = 1,jmid-1
  !  yaux2(j+1) = yaux1(jmid-j)
  !enddo
  !
  !! rewrites yaux1 on the upper side of yaux2
  !do j = jmid+1, jmaxaux
  !  yaux2(j) = yaux1(j-jmid+1)
  !enddo
  !
  !! calculates y on the lower side
  !do j = 1, jmid-1
  !  y(j) = -yaux2(j)
  !enddo
  !
  !! calculates y on the upper side
  !do j = jmid, jmax
  !  y(j) = yaux2(j)
  !enddo
  
  
  ! 1) y direction calculation - non-uniform mesh
  lstop = .false.
  jmaxaux =  jmax
  icountmax = 40
  icount = 1
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
      write(*,*) 'jmax = ', jmax, 'jmaxaux= ', jmaxaux
      allocate (yaux1(jmid))
      allocate (yaux2(jmaxaux))	  
      lstop=.true.
    endif
  enddo
  
  !coeficients for mesh calculation (According to Babucke, 2009)
  a3 = 0.000502236538416527d0
  ksi0 =  0.d0
  
  ! calculating the upper part yaux1
  ksi = 0.d0
  yaux1(1) = ksi
  do j = 1, jmid
    yaux1(j) = a3 * (ksi - ksi0)**3 + (ksi - ksi0)
  ksi = ksi + dn
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
  do j = 1, jmid-1
    y(j) = -yaux2(j)
  enddo
  
  ! calculates y on the upper side
  do j = jmid, jmax
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

!*******************************************************************************
subroutine calc_eta

 ! calculates eta variable as function of y and Reynolds number
 
 implicit none
 include '../var.f90'
 include 'comm.mshbase'
 
 integer(4):: j
 real(8)   :: xpos
 
 xpos = dabs(x(1))

 !open (unit=1,file='eta.dat',status='unknown')
  do j = 1,jmax
    eta(j) = y(j)*dsqrt(Re/xpos)
    !write(1,*) j, eta(j)
  enddo
 !close(1)

return
end subroutine calc_eta

!*******************************************************************************

!*******************************************************************************
!subroutine calc_eta
!
!! calculates eta variable as function of y and Reynolds number
!
!implicit none
!include '../var.f90'
!include 'comm.mshbase'
!
!integer(4):: j
!real(8):: ymax, etamax
!
!ymax = dabs(y(jmax))
!etamax = 10.d0
!
!!open (unit=1,file='eta.dat',status='unknown')
! do j = 1,jmax
!   eta(j) = y(j)*etamax/ymax
!   !write(1,*) j, eta(j)
! enddo
!!close(1)
!
!return
!end subroutine calc_eta

!*******************************************************************************
subroutine integra
implicit none
include '../var.f90'
include 'comm.varbase'
include 'comm.mshbase'

integer(4)          :: j, k, pt
real(8),dimension(5):: f_eta
real(8)             :: h
 
 ! integration of the lower part
 
 ! initial guesses
 f_eta(1) = 0.d0
 f_eta(2) = 0.76503626724495066d0 
 f_eta(3) = 0.17283431626225954d0
 f_eta(4) = 0.93109227001435113d0
 f_eta(5) = 5.06285541989586352d-3
 
 solu(:,jmid) = f_eta
 
 do j = jmid - 1, 1, -1
   h = eta(j) - eta(j+1)

   pt = int ( dabs(h) / 1.d-4)
   h = h / dfloat(pt)

   do k = 1, pt
     call rk4(f_eta, h)
   enddo

   solu(:,j) = f_eta

 enddo
 ! end of integration of the lower part

 
 ! integration of the upper
 
 ! initial guesses
 f_eta(1) = 0.d0
 f_eta(2) = 0.76503626724495066d0 
 f_eta(3) = 0.17283431626225954d0
 f_eta(4) = 0.93109227001435113d0
 f_eta(5) = 5.06285541989586352d-3
 
 do j = jmid + 1, jmax
   h = eta(j) - eta(j-1)

   pt = int ( dabs(h) / 1.d-4)
   h = h / dfloat(pt)

   do k = 1, pt
     call rk4(f_eta, h)
   enddo

   solu(:,j) = f_eta
 enddo
 ! end of integration of the upper part

 open (3, file='perfil.out',status='unknown')
   write(3,*) 'j f f(1) f(2) g g(1) y eta'
   do j = 1, jmax
     write(3,*) j, solu(:,j), y(j), eta(j)
   enddo
 close(3)
 
return
end subroutine integra
!*******************************************************************************
subroutine rk4(f_eta, h)

implicit none
real(8), intent(inout) :: f_eta(5)
real(8), intent(in)    :: h
real(8)                :: hh, h6, dfd_eta(5), f_etat(5), dfd_eta2(5)

   hh = h * 0.5d0
   h6 = h / 6.d0

   call derivs(f_eta, dfd_eta)
   f_etat   = f_eta + hh * dfd_eta

   call derivs(f_etat, dfd_eta2) 
   f_etat   = f_eta + hh * dfd_eta2
   dfd_eta = dfd_eta + 2.d0 * dfd_eta2

   call derivs(f_etat, dfd_eta2)
   f_etat   = f_eta + h * dfd_eta2
   dfd_eta = dfd_eta + 2.d0 * dfd_eta2

   call derivs(f_etat, dfd_eta2)
   f_eta   = f_eta + h6 * ( dfd_eta + dfd_eta2 )

 end subroutine rk4
!*******************************************************************************
subroutine derivs (f_eta, dfd_eta)
  
!  This subroutine sets the system of first order differential 
!  equations for the compressible shear layer 
!  
!  ( Chi*f'' )' + f*f" = 0

!  Momentum equations
!  ( Chi*f'' )' + f*f" = 0 
!   y(1)  = f
!   y(2)  = f'
!   y(3)  = f''
!   dy(3) = f'''

!  Energy equations
!   (( Chi/Pr)*g' )' + f*g' + Chi*(u_e^2/h_e)*( f'' )^2 = 0
!   y(4)  = g
!   y(5)  = g'
!   dy(5) = g''

! Chi = 2.0

implicit none
real(8), intent(in)    :: f_eta(5)
real(8), intent(out)   :: dfd_eta(5)
real(8)                :: Pr, U1, h1, Chi

 Pr = 0.66465d0
 U1 = 155.83d0
 h1 = -9444.04955d0
 Chi = 1.d0
 
 dfd_eta(1) =   f_eta(2)
 dfd_eta(2) =   f_eta(3)
 dfd_eta(3) = - ( f_eta(1) * f_eta(3) )/ Chi
 dfd_eta(4) =   f_eta(5)
 dfd_eta(5) = - (Pr * f_eta(1) * f_eta(5))/Chi - (Pr * U1 * U1 * f_eta(3) * f_eta(3) / h1)

 
 !dfd_eta(1) =   f_eta(2)
 !dfd_eta(2) =   f_eta(3)
 !dfd_eta(3) = - 0.5d0 * f_eta(1) * f_eta(3)
 !dfd_eta(4) =   f_eta(5)
 !dfd_eta(5) = - 0.5d0 * Pr * f_eta(1) * f_eta(5) - Pr * U1 * U1 * f_eta(3) * f_eta(3) / h1


end subroutine derivs
!*******************************************************************************

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
!!*******************************************************************************
