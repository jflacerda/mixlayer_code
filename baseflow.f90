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
 write(*,*) 'generating mesh...'
 call mesh_baseflow
 
 !........................ Calculation of eta variable
 write(*,*) 'calculating eta...' 
 call calc_eta
 
 !........................ Integration of equations
 write(*,*) 'integrating similar equations...' 
 call integra
 
 !........................ Calculation of flow variables
 write(*,*) 'calculating flow variables...' 
 call calc_flowvar 

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
integer(4):: j, jmaxver, i, istf1, istf2, ind
real(8):: a3, ksi0, ksi, stxf, ep, posx, dx, dy, ymax
real(8),allocatable,dimension(:):: yaux1, yaux2
logical:: lstop, lwrite

! logical variable to print mesh
lwrite = .false.

!******************************************************
  !1) y direction calculation - uniform mesh

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
  !! calculate dy
  !ymax = 10.d0
  !dy = ymax/dfloat(jmid-1)
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
!open (unit=1,file='eta.out',status='unknown')
! do j = 1,jmax
!   eta(j) = y(j)*etamax/ymax
!   write(1,*) j, y(j), eta(j)
! enddo
!close(1)
!return
!end subroutine calc_eta
!*******************************************************************************
subroutine calc_eta

 ! calculates eta variable as function of y and Reynolds number
 
 implicit none
 include '../var.f90'
 include 'comm.mshbase'
 
 integer(4):: j
 real(8)   :: xpos
 
 xpos = dabs(x(1))

 !open (unit=1,file='eta.out',status='unknown')
  do j = 1,jmax
    eta(j) = y(j)*dsqrt(Re/xpos)
    !write(1,*) j, eta(j)
  enddo
 !close(1)

return
end subroutine calc_eta

!*******************************************************************************

!*******************************************************************************
subroutine integra
! in this subroutine are searched: f'(0), f''(0), g(0), g'(0) that satisfies
! similar equations up to +inf and -inf simultaneously

! The search process is:
! (a) guess f''(0) and g'(0) and search f'(0+) and g(0+) through integration 
!     that satisfies eqs at +inf
! (b) guess f''(0) and g'(0) and search f'(0-) and g(0-) through integration 
!     that satisfies eqs at -inf
! (c) repeat (a) and (b) modifying now f''(0) and g'(0) until:
!     (|f'(0+) - f'(0-)'| < tol) .and. (|g(0+) - g(0-)| < tol)

implicit none
include '../var.f90'
include 'comm.varbase'
include 'comm.mshbase'

real(8):: R1, Q1, R2, Q2, tol, fR, fQ, flR, flQ
real(8):: L1, P1, L2, P2, swap, Rnew, Rold, Qnew, Qold
real(8):: dxR, dxQ
integer(4):: icond, icondR, icondQ, iconditer, iter, itmax, j

parameter(tol   = 1.d-8)
parameter(itmax = 40)

open(unit=1,file='secant.out',status='unknown')

 write(1,*)
 write(1,*)' The initial boundary conditions at eta(0) are:'
 write(1,*)
 write(1,*)'f_i (0) = ',Lguess,'f_ii (0) = ',Rguess
 write(1,*)'g   (0) = ',Pguess,'g_i  (0) = ',Qguess
 write(1,*) 

! 1) Initial search values -----------------
     ! first guess
     R1 = Rguess
     Q1 = Qguess
     
     ! integration of equations to + inf
     call integra_plus(R1,Q1,L1,P1)
	 
     ! integration of equations to - inf
     call integra_minus(R1,Q1,L2,P2)
	 
     ! compare difference between L1 and L2
     flR = L1 - L2
     flQ = P1 - P2
     
     ! second guess   
     R2 = R1 + 0.01d0
     Q2 = Q1 + 0.01d0
     
     ! integration of equations to + inf
     call integra_plus(R2,Q2,L1,P1)
     
     ! integration of equations to - inf
     call integra_minus(R2,Q2,L2,P2)
     
     ! compare difference between L1 and L2
     fR = L1 - L2
     fQ = P1 - P2
     
     ! comparing first and second guesses
     ! for R
     if (dabs(flR).lt.dabs(fR)) then ! pick the bound with smaller value
       Rnew   = R1
       Rold   = R2
       swap   = flR
       flR    = fR
       fR     = swap
     else
       Rnew   = R2
       Rold   = R1
     endif
     
     ! for Q
     if (dabs(flQ).lt.dabs(fQ)) then ! pick the bound with smaller value
       Qnew = Q1
       Qold   = Q2
       swap   = flQ
       flQ    = fQ
       fQ     = swap
     else
       Qnew = Q2
       Qold   = Q1
     endif
! ----- End of Initial search values ----------

! 2) Iteration process ------------------------
     iter = 1
     
     icondR = 1 ! condition for R
     icondQ = 1 ! condition for Q
     iconditer = 1 ! condition for max iteration
     icond = (icondR + icondQ) * iconditer
	 
     write(1,*)'************************************************'
     write(1,*)
     write(1,*)'subroutine secant - iterations in R,Q'
     write(1,*)
     write(1,*)'************************************************'	 

	 it_process: do while (icond .ne. 0)
        ! Alteration of R and Q
        dxR = (Rold-Rnew)*fR/(fR-flR) !R
        dxQ = (Qold-Qnew)*fQ/(fQ-flQ) !Q
        
        Rold = Rnew
        Qold = Qnew
        
        flR = fR
        flQ = fQ
        
        Rnew = Rnew + dxR !R
        Qnew = Qnew + dxQ !Q
        
        ! integration of equations to + inf
        call integra_plus(Rnew,Qnew,L1,P1)
        
        ! integration of equations to - inf
        call integra_minus(Rnew,Qnew,L2,P2)
        
        ! compare difference between L1 and L2
        fR = L1 - L2
        fQ = P1 - P2
		
		! Printing results
        write(1,*)'************************************************'
        write(1,*)
        write(1,*)'iterat =',iter
        write(1,10)'R =',Rnew,'L(+inf)=',L1,'L(-inf)=',L2,'f_i (+inf) - f_i (-inf)= ',fR
        write(1,10)'Q =',Qnew,'P(+inf)=',P1,'P(-inf)=',P2,'  g (+inf) - g   (-inf)= ',fQ		
        write(1,*)
        write(1,*)'************************************************'
        10 format(a3,2x,F15.10,2x,a8,2x,F15.10,2x,a8,2x,F15.10,2x,a25,2x,F15.10)
		
		! Evaluating conditions
		if (dabs(fR).lt.tol) icondR = 0 ! R condition
        if (dabs(fQ).lt.tol) icondQ = 0 ! Q condition
        if (iter.gt.itmax) iconditer = 0 ! iter max condition		
        icond = (icondR + icondQ) * iconditer
		
		iter = iter + 1	

     enddo it_process
     
	 if (iter.gt.itmax) then
       write(1,*)	 
       write(1,*)'Maximum iterations exceeded'
       write(1,*)
     else 
       write(1,*)
       write(1,*)'Secant method for total integration converged in',iter-1,' iterations'
       write(1,*)
     endif
	 
	 ! Writes final boundary conditions on secant.out
     write(1,*)
     write(1,*)' The boundary conditions at eta(0) are:'
     write(1,*)
     write(1,*)'f_i (0) = ',L2,'f_ii (0) = ',Rnew
     write(1,*)'g   (0) = ',P2,'g_i  (0) = ',Qnew
! ---------------------------------------------
close(1)
 

 open (3, file='perfil.out',status='unknown')
   write(3,*) 'j f f(1) f(2) g g(1) y eta'
   do j = 1, jmax
     write(3,*) j, solu(:,j), y(j), eta(j)
   enddo
 close(3)
 
return
end subroutine integra
!*******************************************************************************

!*******************************************************************************
subroutine integra_plus(R,Q,L,P)
! in this subroutine are performed integration to +inf
! searched: f'(0) and g(0) that satisfies similar equations

! The search process is:

implicit none
include '../var.f90'
include 'comm.varbase'
include 'comm.mshbase'

real(8):: R, Q, tol, fL, fP, flL, flP
real(8):: L, L1, L2, P, P1, P2, swap, Lold, Pold
real(8):: dxL, dxP
integer(4):: icond, icondL, icondP, iconditer, iter, itmax

parameter(tol   = 1.d-8)
parameter(itmax = 40) 

! 1) Initial search values -----------------
     ! first guess
     L1 = Lguess
     P1 = Pguess
   
     ! integration of equations to + inf
     call int_plus(R,Q,L1,P1)
	 flL = (solu(2,jmax)-1.d0)
     flP = (solu(4,jmax)-1.d0)
	 
	 ! second guess	 
     L2 = L1 + dabs(solu(2,jmax)-1.d0)
     P2 = P1 + dabs(solu(4,jmax)-1.d0)
	 
     13 format (a3,x,E20.13)
     write(1,*) '**************************'	 
	 write(1,*) 'subroutine integra_plus'	 
	 write(1,13) 'R= ', R
	 write(1,13) 'Q= ', Q
	 write(1,13) 'L1=', L1
	 write(1,13) 'P1=', P1	 
	 write(1,13) 'L2=', L2
	 write(1,13) 'P2=', P2	 
     write(1,*) '**************************'
	 write(1,*)	 

	 ! integration of equations to + inf
     call int_plus(R,Q,L2,P2)
     fL = (solu(2,jmax)-1.d0)
     fP = (solu(4,jmax)-1.d0)
	 
     ! comparing first and second guesses
     ! for L
     if (dabs(flL).lt.dabs(fL)) then ! pick the bound with smaller value
       L   = L1
       Lold   = L2
       swap   = flL
       flL    = fL
       fL     = swap
     else
       L   = L2
       Lold   = L1
     endif
     
     ! for P
     if (dabs(flP).lt.dabs(fP)) then ! pick the bound with smaller value
       P   = P1
       Pold   = P2
       swap   = flP
       flP    = fP
       fP     = swap
     else
       P = P2
       Pold   = P1
     endif
! ----- End of Initial search values ----------

! 2) Iteration process ------------------------
     iter = 1
     
     icondL = 1 ! condition for L
     icondP = 1 ! condition for P
     iconditer = 1 ! condition for max iteration
     icond = (icondL + icondP) * iconditer
	 
     !open(unit=2,file='integra_plus.out',status='unknown')
	 
     write(1,*)'************************************************'
     write(1,*)
     write(1,*)'subroutine integra_plus - iterations in L,P'
     write(1,*)
     write(1,*)'************************************************'	 

	 it_process_plus: do while (icond .ne. 0)
        ! Alteration of L and P
        dxL = (Lold-L)*fL/(fL-flL) !L
        dxP = (Pold-P)*fP/(fP-flP) !P
        
        Lold = L
        Pold = P
        
        flL = fL
        flP = fP
        
        L = L + dxL !L
        P = P + dxP !P
        
        ! integration of equations to + inf
        call int_plus(R,Q,L,P)
        fL = (solu(2,jmax)-1.d0)
        fP = (solu(4,jmax)-1.d0)		
        
		! Printing results
        write(1,*)'************************************************'
        write(1,*)
        write(1,*)'iterat =',iter
        write(1,11)'L =',L,'R =',R,'f_i (+inf) - 1 = ',fL
        write(1,11)'P =',P,'Q =',Q,'g   (+inf) - 1 = ',fP
        write(1,*)'************************************************'
        11 format(a3,2x,F15.10,2x,a3,2x,F15.10,2x,a17,2x,F15.10)
		
		! Evaluating conditions
        if (dabs(fL).lt.tol) icondL = 0 ! L condition
        if (dabs(fP).lt.tol) icondP = 0 ! P condition
        if (iter.gt.itmax) iconditer = 0 ! iter max condition		
        icond = (icondL + icondP) * iconditer
		
		iter = iter + 1	

     enddo it_process_plus
     
	 if (iter.gt.itmax) then
       write(1,*)	 
       write(1,*)'Maximum iterations exceeded'
       write(1,*)
     else    
       write(1,*)
       write(1,*)'Secant method for integra_plus converged in',iter-1,' iterations'
       write(1,*)
     endif
	 
	 ! Writes final boundary conditions on integra_plus.out
     write(1,*)
     write(1,*)' The boundary conditions at eta(0) are:'
     write(1,*)
     write(1,*)'f_i (0) = ',L,'f_ii (0) = ',R
     write(1,*)'g   (0) = ',P,'g_i  (0) = ',Q
     
     !close(2)
! ---------------------------------------------
 
return
end subroutine integra_plus
!*******************************************************************************
subroutine int_plus(R,Q,L,P)
implicit none
! This subroutine integrates the upper part
include '../var.f90'
include 'comm.varbase'
include 'comm.mshbase'

integer(4)          :: j, k, pt
real(8),dimension(5):: f_eta
real(8)             :: h, R, Q, L, P

 ! initial guesses
 f_eta(1) = 0.d0
 f_eta(2) = L
 f_eta(3) = R
 f_eta(4) = P
 f_eta(5) = Q
 
 solu(:,jmid) = f_eta
 
 do j = jmid + 1, jmax
   h = eta(j) - eta(j-1)
 
   pt = int ( dabs(h) / 1.d-4)
   h = h / dfloat(pt)
 
   do k = 1, pt
     call rk4(f_eta, h)
   enddo
 
   solu(:,j) = f_eta
 enddo

return
end subroutine int_plus
!*******************************************************************************
subroutine integra_minus(R,Q,L,P)
! in this subroutine are performed integration to -inf
! searched: f'(0) and g(0) that satisfies similar equations

! The search process is:

implicit none
include '../var.f90'
include 'comm.varbase'
include 'comm.mshbase'

real(8):: R, Q, tol, fL, fP, flL, flP
real(8):: L, L1, L2, P, P1, P2, swap, Lold, Pold
real(8):: dxL, dxP
integer(4):: icond, icondL, icondP, iconditer, iter, itmax

parameter(tol   = 1.d-8)
parameter(itmax = 40) 

! 1) Initial search values -----------------
     ! first guess
     L1 = Lguess
     P1 = Pguess
	 
     ! integration of equations to - inf
     call int_minus(R,Q,L1,P1)
	 flL = (solu(2,1)-Uratio)
     flP = (solu(4,1)-Heratio)
	 
	 ! second guess	 
     L2 = L1 + dabs(solu(2,1)- Uratio)
     P2 = P1 + dabs(solu(4,1)- Heratio)
	 
     14 format (a3,x,E20.13)
     write(1,*) '**************************'	 
	 write(1,*) 'subroutine integra_minus'	 
	 write(1,14) 'R= ', R
	 write(1,14) 'Q= ', Q
	 write(1,14) 'L1=', L1
	 write(1,14) 'P1=', P1	 
	 write(1,14) 'L2=', L2
	 write(1,14) 'P2=', P2	 
     write(1,*) '**************************'
	 write(1,*)	

	 ! integration of equations to - inf
     call int_minus(R,Q,L2,P2)
     fL = (solu(2,1)-Uratio)
     fP = (solu(4,1)-Heratio)
	 
     ! comparing first and second guesses
     ! for L
     if (dabs(flL).lt.dabs(fL)) then ! pick the bound with smaller value
       L   = L1
       Lold   = L2
       swap   = flL
       flL    = fL
       fL     = swap
     else
       L   = L2
       Lold   = L1
     endif
     
     ! for P
     if (dabs(flP).lt.dabs(fP)) then ! pick the bound with smaller value
       P   = P1
       Pold   = P2
       swap   = flP
       flP    = fP
       fP     = swap
     else
       P = P2
       Pold   = P1
     endif
! ----- End of Initial search values ----------

! 2) Iteration process ------------------------
     iter = 1
     
     icondL = 1 ! condition for L
     icondP = 1 ! condition for P
     iconditer = 1 ! condition for max iteration
     icond = (icondL + icondP) * iconditer
	 
     !open(unit=3,file='integra_minus.out',status='unknown')
	 
     write(1,*)'************************************************'
     write(1,*)
     write(1,*)'subroutine integra_minus - iterations in L,P'
     write(1,*)
     write(1,*)'************************************************'	 
	 
	 it_process_minus: do while (icond .ne. 0)
        ! Alteration of L and P
        dxL = (Lold-L)*fL/(fL-flL) !L
        dxP = (Pold-P)*fP/(fP-flP) !P
        
        Lold = L
        Pold = P
        
        flL = fL
        flP = fP
        
        L = L + dxL !L
        P = P + dxP !P
        
        ! integration of equations to - inf
        call int_minus(R,Q,L,P)
        fL = (solu(2,1)-Uratio)
        fP = (solu(4,1)-Heratio)		
        
        ! Printing results
        write(1,*)'************************************************'
        write(1,*)
        write(1,*)'iterat =',iter
        write(1,12)'L =',L,'R =',R,'f_i (-inf) - Uratio = ',fL
        write(1,12)'P =',P,'Q =',Q,'g   (-inf) - Heratio= ',fP
        write(1,*)'************************************************'
        12 format(a3,2x,F15.10,2x,a3,2x,F15.10,2x,a27,2x,F15.10)
		
		! Evaluating conditions
        if (dabs(fL).lt.tol) icondL = 0 ! L condition
        if (dabs(fP).lt.tol) icondP = 0 ! P condition
        if (iter.gt.itmax) iconditer = 0 ! iter max condition
        icond = (icondL + icondP) * iconditer
		
		iter = iter + 1	

     enddo it_process_minus
     
	 if (iter.gt.itmax) then
       write(1,*)	 
       write(1,*)'Maximum iterations exceeded'
       write(1,*)
     else    
       write(1,*)
       write(1,*)'Secant method for integra_minus converged in',iter-1,' iterations'
       write(1,*)
     endif
	 
	 ! Writes final boundary conditions on integra_minus.out
     write(1,*)
     write(1,*)' The boundary conditions at eta(0) are:'
     write(1,*)
     write(1,*)'f_i (0) = ',L,'f_ii (0) = ',R
     write(1,*)'g   (0) = ',P,'g_i  (0) = ',Q
     
     !close(3)
! ---------------------------------------------
 
return
end subroutine integra_minus
!*******************************************************************************
subroutine int_minus(R,Q,L,P)
implicit none
! This subroutine integrates the lower part
include '../var.f90'
include 'comm.varbase'
include 'comm.mshbase'

integer(4)          :: j, k, pt
real(8),dimension(5):: f_eta
real(8)             :: h, R, Q, L, P

 ! initial guesses
 f_eta(1) = 0.d0
 f_eta(2) = L
 f_eta(3) = R
 f_eta(4) = P
 f_eta(5) = Q
 
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
 
return
end subroutine int_minus
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
include '../var.f90'
real(8), intent(in)    :: f_eta(5)
real(8), intent(out)   :: dfd_eta(5)
real(8)                :: Pr1, U1, h1, Chi

 Pr1 = Pr
 U1  = Umax1
 h1  = Hemax1
 Chi = 2.d0
 
 dfd_eta(1) =   f_eta(2)
 dfd_eta(2) =   f_eta(3)
 dfd_eta(3) = - ( f_eta(1) * f_eta(3) )/ Chi
 dfd_eta(4) =   f_eta(5)
 dfd_eta(5) = - (Pr1 * f_eta(1) * f_eta(5))/Chi - (Pr1 * U1 * U1 * f_eta(3) * f_eta(3) / h1)

 
 !dfd_eta(1) =   f_eta(2)
 !dfd_eta(2) =   f_eta(3)
 !dfd_eta(3) = - 0.5d0 * f_eta(1) * f_eta(3)
 !dfd_eta(4) =   f_eta(5)
 !dfd_eta(5) = - 0.5d0 * Pr1 * f_eta(1) * f_eta(5) - Pr1 * U1 * U1 * f_eta(3) * f_eta(3) / h1


end subroutine derivs
!*******************************************************************************
subroutine calc_flowvar

 ! calculates non-dimensionles flow variables
 ! u   - x velocity
 ! v   - y velocity
 ! rho - density
 ! T   - Temperature
 
 implicit none
 include '../var.f90'
 include 'comm.mshbase'
 include 'comm.varbase'
 
 integer(4):: j
 
 ! u: x velocity
 u = solu(2,:)
 
 ! v: y velocity
 v = eta*solu(2,:) - solu(1,:)
 
 ! T: temperature
 T = solu(4,:)
 
 ! rho: density
 rho = 1.d0/T
 
 open (1, file='results.out',status='unknown')
   write(1,*) 'j u v rho T y eta'
   do j = 1, jmax
     write(1,*) j, u(j), v(j), rho(j), T(j), y(j), eta(j)
   enddo
 close(1)

return
end subroutine calc_flowvar

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
