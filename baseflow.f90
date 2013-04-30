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
!     Velocities, Density and Temperature profiles:
!
!     f' - velocity - u/u1
!     g  - enthalpy - He/He1

! Jonatas F Lacerda    jonatasflacerda.hotmail.com
! ICMC - USP

! date: 14.04.2013

implicit none
include '../var.f90'
include 'comm.mshbase'
include 'comm.varbase'

 !........................ Mesh generation
 write(*,*) 'generating mesh...'
 call mesh_baseflow

 !........................ Calculation of eta variable
 write(*,*) 'calculating eta...'
 call calc_eta(x(1))

 !........................ Finding the BC for EDO's
 write(*,*) 'finding boundary condition for EDOs at x_0...'
 call bc_finder
 
 !........................ Calculation of flow variables
 write(*,*) 'calculating flow variables...'
 call calc_flowvar
 
 !........................ Calculation of vorticity thickness
 write(*,*) 'calculating vorticity thickness...'
 call calc_vorticitythickness
 
 !........................ Integration of equations
 write(*,*) 'Integrating the equations in all domain'
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
integer(4):: j, jmaxver, i, istf1, istf2, ind
real(8):: a3, ksi0, ksi, stxf, ep, posx, dx, dy, ymax
real(8),allocatable,dimension(:):: yaux1, yaux2
logical:: lstop, lwrite

! logical variable to print mesh
lwrite = .false.

!******************************************************
  !1) y direction calculation - uniform mesh

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

   ! calculate dy
   ymax = 10.d0
   dy = ymax/dfloat(jmid-1)

   ! calculating the upper part yaux1
   ksi = 0.d0
   do j = 1, jmid
     yaux1(j) = ksi
     ksi = ksi + dy
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


  ! 1) y direction calculation - non-uniform mesh
! lstop = .false.
! jmaxaux =  jmax
! icountmax = 40
! icount = 1
! do while (.not.lstop)
!   ! calculating middle value of y-index: jmid
!   jmid = (jmaxaux+1)/2
!
!   ! calculating index verification
!   jmaxver = (2*jmid)-1
!
!   !verifying if jmid is correct
!   if (jmaxver.ne.jmaxaux) then
!     jmaxaux = jmaxaux + 1
!     icount = icount + 1
!     if (icount.gt.icountmax) then
!       write(*,*) 'not worked in', icountmax,' iterations'
!       lstop=.true.
!     endif
!   else
!     allocate (yaux1(jmid))
!     allocate (yaux2(jmaxaux))
!     lstop=.true.
!   endif
! enddo
!
! !coeficients for mesh calculation (According to Babucke, 2009)
! a3 = 0.000502236538416527d0
! ksi0 =  0.d0
!
! ! calculating the upper part yaux1
! ksi = 0.d0
! yaux1(1) = ksi
! do j = 1, jmid
!   yaux1(j) = a3 * (ksi - ksi0)**3 + (ksi - ksi0)
! ksi = ksi + dn
! enddo
!
! ! copying yaux1 to the lower part yaux2
! yaux2(1) = yaux1(jmid)
! do j = 1,jmid-1
!   yaux2(j+1) = yaux1(jmid-j)
! enddo
!
! ! rewrites yaux1 on the upper side of yaux2
! do j = jmid+1, jmaxaux
!   yaux2(j) = yaux1(j-jmid+1)
! enddo
!
! ! calculates y on the lower side
! do j = 1, jmid-1
!   y(j) = -yaux2(j)
! enddo
!
! ! calculates y on the upper side
! do j = jmid, jmax
!   y(j) = yaux2(j)
! enddo

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
! subroutine calc_eta
!
! ! calculates eta variable as function of y and Reynolds number
!
! implicit none
! include '../var.f90'
! include 'comm.mshbase'
!
! integer(4):: j
! real(8):: ymax, etamax
!
! ymax = dabs(y(jmax))
! etamax = 10.d0
!
! open (unit=1,file='eta.out',status='unknown')
!  do j = 1,jmax
!    eta(j) = y(j)*etamax/ymax
!    write(1,*) j, y(j), eta(j)
!  enddo
! close(1)
! return
! end subroutine calc_eta
!*******************************************************************************
subroutine calc_eta(xpos)

 ! calculates eta variable as function of y and Reynolds number

 implicit none
 include '../var.f90'
 include 'comm.mshbase'

 integer(4):: j
 real(8)   :: xpos

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
subroutine bc_finder
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

real(8):: fpp1, gp1, fpp2, gp2, tol, ffpp, fgp, flfpp, flgp
real(8):: fp1, g1, fp2, g2, swap, fppnew, fppold, gpnew, gpold
real(8):: dxfpp, dxgp
integer(4):: icond, icondfpp, icondgp, iconditer, iter, itmax, j

parameter(tol   = 1.d-8)
parameter(itmax = 40)

open(unit=1,file='secant.out',status='unknown')

 write(1,*)
 write(1,*)' The initial boundary conditions at eta(0) are:'
 write(1,*)
 write(1,*)'f_i (0) = ',fp_guess,'f_ii (0) = ',fpp_guess
 write(1,*)'g   (0) = ',g_guess,'g_i  (0) = ',gp_guess
 write(1,*)

! 1) Initial search values -----------------
     ! first guess
     fpp1 = fpp_guess
     gp1  = gp_guess

     ! integration of equations to + inf
     call integra_plus(fpp1,gp1,fp1,g1)

     ! integration of equations to - inf
     call integra_minus(fpp1,gp1,fp2,g2)

     ! compare difference between fp1 and fp2
     flfpp = fp1 - fp2
     flgp  = g1 - g2

     ! second guess
     fpp2 = fpp1 + 0.01d0
     gp2  = gp1 + 0.01d0

     ! integration of equations to + inf
     call integra_plus(fpp2,gp2,fp1,g1)

     ! integration of equations to - inf
     call integra_minus(fpp2,gp2,fp2,g2)

     ! compare difference between fp1 and fp2
     ffpp = fp1 - fp2
     fgp  = g1 - g2

     ! comparing first and second guesses
     ! for fpp
     if (dabs(flfpp).lt.dabs(ffpp)) then ! pick the bound with smaller value
       fppnew = fpp1
       fppold = fpp2
       swap   = flfpp
       flfpp  = ffpp
       ffpp   = swap
     else
       fppnew = fpp2
       fppold = fpp1
     endif

     ! for gp
     if (dabs(flgp).lt.dabs(fgp)) then ! pick the bound with smaller value
       gpnew = gp1
       gpold = gp2
       swap  = flgp
       flgp  = fgp
       fgp   = swap
     else
       gpnew = gp2
       gpold = gp1
     endif
! ----- End of Initial search values ----------

! 2) Iteration process ------------------------
     iter = 1

     icondfpp  = 1 ! condition for fpp
     icondgp   = 1 ! condition for gp
     iconditer = 1 ! condition for max iteration
     icond     = (icondfpp + icondgp) * iconditer

     write(1,*)'************************************************'
     write(1,*)
     write(1,*)'subroutine secant - iterations in fpp,gp'
     write(1,*)
     write(1,*)'************************************************'

     it_process: do while (icond .ne. 0)
        ! Alteration of fpp and gp
        dxfpp = (fppold-fppnew)*ffpp/(ffpp-flfpp) !fpp
        dxgp  = (gpold-gpnew)*fgp/(fgp-flgp) !gp

        fppold = fppnew
        gpold  = gpnew

        flfpp = ffpp
        flgp  = fgp

        fppnew = fppnew + dxfpp !fpp
        gpnew  = gpnew  + dxgp  !gp

        ! integration of equations to + inf
        call integra_plus(fppnew,gpnew,fp1,g1)

        ! integration of equations to - inf
        call integra_minus(fppnew,gpnew,fp2,g2)

        ! compare difference between fp1 and fp2
        ffpp = fp1 - fp2
        fgp  = g1 - g2

        ! Printing results
        write(1,*)'************************************************'
        write(1,*)
        write(1,*)'iterat =',iter
        write(1,10)'fpp =',fppnew,'fp(+inf)=',fp1,'fp(-inf)=',fp2,'f_i (+inf) - f_i (-inf)= ',ffpp
        write(1,10)'gp =',gpnew,'g(+inf)=',g1,'g(-inf)=',g2,'  g (+inf) - g   (-inf)= ',fgp
        write(1,*)
        write(1,*)'************************************************'
        10 format(a3,2x,F15.10,2x,a8,2x,F15.10,2x,a8,2x,F15.10,2x,a25,2x,F15.10)

        ! Evaluating conditions
        if (dabs(ffpp).lt.tol) icondfpp = 0 ! fpp condition
        if (dabs(fgp).lt.tol) icondgp = 0   ! gp condition
        if (iter.gt.itmax) iconditer = 0    ! iter max condition
        icond = (icondfpp + icondgp) * iconditer

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
     fp_f  = fp2
     fpp_f = fppnew
     g_f   = g2
     gp_f  = gpnew
     write(1,*)
     write(1,*)' The boundary conditions at eta(0) are:'
     write(1,*)
     write(1,*)'f_i (0) = ',fp_f,'f_ii (0) = ',fpp_f
     write(1,*)'g   (0) = ',g_f,'g_i  (0) = ',gp_f
! ---------------------------------------------
close(1)


 open (3, file='perfil.out',status='unknown')
   write(3,*) 'j f fp fpp g gp y eta'
   do j = 1, jmax
     write(3,*) j, solu(:,j), y(j), eta(j)
   enddo
 close(3)

return
end subroutine bc_finder
!*******************************************************************************

!*******************************************************************************
subroutine integra_plus(fpp,gp,fp,g)
! in this subroutine are performed integration to +inf
! searched: f'(0) and g(0) that satisfies similar equations

! The search process is:

implicit none
include '../var.f90'
include 'comm.varbase'
include 'comm.mshbase'

real(8):: fpp, gp, tol, ffp, fg, flfp, flg
real(8):: fp, fp1, fp2, g, g1, g2, swap, fpold, gold
real(8):: dxfp, dxg
integer(4):: icond, icondfp, icondg, iconditer, iter, itmax

parameter(tol   = 1.d-8)
parameter(itmax = 40)

! 1) Initial search values -----------------
     ! first guess
     fp1 = fp_guess
     g1 = g_guess

     ! integration of equations to + inf
     call int_plus(fpp,gp,fp1,g1)
     flfp = (solu(2,jmax)-1.d0)
     flg = (solu(4,jmax)-1.d0)

     ! second guess
     fp2 = fp1 + dabs(solu(2,jmax)-1.d0)
     g2 = g1 + dabs(solu(4,jmax)-1.d0)

     13 format (a3,x,E20.13)
     write(1,*) '**************************'
     write(1,*) 'subroutine integra_plus'
     write(1,13) 'fpp= ', fpp
     write(1,13) 'gp= ', gp
     write(1,13) 'fp1=', fp1
     write(1,13) 'g1=', g1
     write(1,13) 'fp2=', fp2
     write(1,13) 'g2=', g2
     write(1,*) '**************************'
     write(1,*)

     ! integration of equations to + inf
     call int_plus(fpp,gp,fp2,g2)
     ffp = (solu(2,jmax)-1.d0)
     fg  = (solu(4,jmax)-1.d0)

     ! comparing first and second guesses
     ! for fp
     if (dabs(flfp).lt.dabs(ffp)) then ! pick the bound with smaller value
       fp    = fp1
       fpold = fp2
       swap  = flfp
       flfp  = ffp
       ffp   = swap
     else
       fp   = fp2
       fpold   = fp1
     endif

     ! for g
     if (dabs(flg).lt.dabs(fg)) then ! pick the bound with smaller value
       g    = g1
       gold = g2
       swap = flg
       flg  = fg
       fg   = swap
     else
       g = g2
       gold   = g1
     endif
! ----- End of Initial search values ----------

! 2) Iteration process ------------------------
     iter = 1

     icondfp = 1 ! condition for fp
     icondg = 1 ! condition for g
     iconditer = 1 ! condition for max iteration
     icond = (icondfp + icondg) * iconditer

     !open(unit=2,file='integra_plus.out',status='unknown')

     write(1,*)'************************************************'
     write(1,*)
     write(1,*)'subroutine integra_plus - iterations in fp,g'
     write(1,*)
     write(1,*)'************************************************'

     it_process_plus: do while (icond .ne. 0)
        ! Alteration of fp and g
        dxfp = (fpold-fp)*ffp/(ffp-flfp) !fp
        dxg = (gold-g)*fg/(fg-flg) !g

        fpold = fp
        gold = g

        flfp = ffp
        flg = fg

        fp = fp + dxfp !fp
        g = g + dxg !g

        ! integration of equations to + inf
        call int_plus(fpp,gp,fp,g)
        ffp = (solu(2,jmax)-1.d0)
        fg = (solu(4,jmax)-1.d0)

        ! Printing results
        write(1,*)'************************************************'
        write(1,*)
        write(1,*)'iterat =',iter
        write(1,11)'fp =',fp,'fpp =',fpp,'f_i (+inf) - 1 = ',ffp
        write(1,11)'g =',g,'gp =',gp,'g   (+inf) - 1 = ',fg
        write(1,*)'************************************************'
        11 format(a3,2x,F15.10,2x,a3,2x,F15.10,2x,a17,2x,F15.10)

        ! Evaluating conditions
        if (dabs(ffp).lt.tol) icondfp = 0 ! fp condition
        if (dabs(fg).lt.tol) icondg = 0 ! g condition
        if (iter.gt.itmax) iconditer = 0 ! iter max condition
        icond = (icondfp + icondg) * iconditer

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
     write(1,*)'f_i (0) = ',fp,'f_ii (0) = ',fpp
     write(1,*)'g   (0) = ',g,'g_i  (0) = ',gp

     !close(2)
! ---------------------------------------------

return
end subroutine integra_plus
!*******************************************************************************
subroutine int_plus(fpp,gp,fp,g)
implicit none
! This subroutine integrates the upper part
include '../var.f90'
include 'comm.varbase'
include 'comm.mshbase'

integer(4)          :: j, k, pt
real(8),dimension(5):: f_eta
real(8)             :: h, fpp, gp, fp, g

 ! initial guesses
 f_eta(1) = 0.d0
 f_eta(2) = fp
 f_eta(3) = fpp
 f_eta(4) = g
 f_eta(5) = gp

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
subroutine integra_minus(fpp,gp,fp,g)
! in this subroutine are performed integration to -inf
! searched: f'(0) and g(0) that satisfies similar equations

! The search process is:

implicit none
include '../var.f90'
include 'comm.varbase'
include 'comm.mshbase'

real(8):: fpp, gp, tol, ffp, fg, flfp, flg
real(8):: fp, fp1, fp2, g, g1, g2, swap, fpold, gold
real(8):: dxfp, dxg
integer(4):: icond, icondfp, icondg, iconditer, iter, itmax

parameter(tol   = 1.d-8)
parameter(itmax = 40)

! 1) Initial search values -----------------
     ! first guess
     fp1 = fp_guess
     g1 = g_guess

     ! integration of equations to - inf
     call int_minus(fpp,gp,fp1,g1)
     flfp = (solu(2,1)-Uratio)
     flg  = (solu(4,1)-Heratio)

     ! second guess
     fp2 = fp1 + dabs(solu(2,1)- Uratio)
     g2  = g1 + dabs(solu(4,1)- Heratio)

     14 format (a3,x,E20.13)
     write(1,*) '**************************'
     write(1,*) 'subroutine integra_minus'
     write(1,14) 'fpp= ', fpp
     write(1,14) 'gp= ', gp
     write(1,14) 'fp1=', fp1
     write(1,14) 'g1=', g1
     write(1,14) 'fp2=', fp2
     write(1,14) 'g2=', g2
     write(1,*) '**************************'
     write(1,*)

     ! integration of equations to - inf
     call int_minus(fpp,gp,fp2,g2)
     ffp = (solu(2,1)-Uratio)
     fg = (solu(4,1)-Heratio)

     ! comparing first and second guesses
     ! for fp
     if (dabs(flfp).lt.dabs(ffp)) then ! pick the bound with smaller value
       fp    = fp1
       fpold = fp2
       swap  = flfp
       flfp  = ffp
       ffp   = swap
     else
       fp    = fp2
       fpold = fp1
     endif

     ! for g
     if (dabs(flg).lt.dabs(fg)) then ! pick the bound with smaller value
       g    = g1
       gold = g2
       swap = flg
       flg  = fg
       fg   = swap
     else
       g    = g2
       gold = g1
     endif
! ----- End of Initial search values ----------

! 2) Iteration process ------------------------
     iter = 1

     icondfp   = 1 ! condition for fp
     icondg    = 1 ! condition for g
     iconditer = 1 ! condition for max iteration
     icond     = (icondfp + icondg) * iconditer

     !open(unit=3,file='integra_minus.out',status='unknown')

     write(1,*)'************************************************'
     write(1,*)
     write(1,*)'subroutine integra_minus - iterations in fp,g'
     write(1,*)
     write(1,*)'************************************************'

     it_process_minus: do while (icond .ne. 0)
        ! Alteration of fp and g
        dxfp = (fpold-fp)*ffp/(ffp-flfp) !fp
        dxg  = (gold-g)*fg/(fg-flg) !g

        fpold = fp
        gold  = g

        flfp = ffp
        flg  = fg

        fp = fp + dxfp !fp
        g  = g + dxg !g

        ! integration of equations to - inf
        call int_minus(fpp,gp,fp,g)
        ffp = (solu(2,1)-Uratio)
        fg  = (solu(4,1)-Heratio)

        ! Printing results
        write(1,*)'************************************************'
        write(1,*)
        write(1,*)'iterat =',iter
        write(1,12)'fp =',fp,'fpp =',fpp,'f_i (-inf) - Uratio = ',ffp
        write(1,12)'g =',g,'gp =',gp,'g   (-inf) - Heratio= ',fg
        write(1,*)'************************************************'
        12 format(a3,2x,F15.10,2x,a3,2x,F15.10,2x,a27,2x,F15.10)

        ! Evaluating conditions
        if (dabs(ffp).lt.tol) icondfp = 0 ! fp condition
        if (dabs(fg).lt.tol) icondg = 0 ! g condition
        if (iter.gt.itmax) iconditer = 0 ! iter max condition
        icond = (icondfp + icondg) * iconditer

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
     write(1,*)'f_i (0) = ',fp,'f_ii (0) = ',fpp
     write(1,*)'g   (0) = ',g,'g_i  (0) = ',gp

     !close(3)
! ---------------------------------------------

return
end subroutine integra_minus
!*******************************************************************************
subroutine int_minus(fpp,gp,fp,g)
implicit none
! This subroutine integrates the lower part
include '../var.f90'
include 'comm.varbase'
include 'comm.mshbase'

integer(4)          :: j, k, pt
real(8),dimension(5):: f_eta
real(8)             :: h, fpp, gp, fp, g

 ! initial guesses
 f_eta(1) = 0.d0
 f_eta(2) = fp
 f_eta(3) = fpp
 f_eta(4) = g
 f_eta(5) = gp

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
!   (( Chi/gr)*g' )' + f*g' + Chi*(u_e^2/h_e)*( f'' )^2 = 0
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
!dfd_eta(5) = - (Pr1 * f_eta(1) * f_eta(5))/Chi - 2.d0 * Pr1 * f_eta(3) * f_eta(3)


 !dfd_eta(1) =   f_eta(2)
 !dfd_eta(2) =   f_eta(3)
 !dfd_eta(3) = - 0.5d0 * f_eta(1) * f_eta(3)
 !dfd_eta(4) =   f_eta(5)
 !dfd_eta(5) = - 0.5d0 * Pr1 * f_eta(1) * f_eta(5) - Pr1 * U1 * U1 * f_eta(3) * f_eta(3) / h1


end subroutine derivs
!*******************************************************************************
subroutine integra
 ! integrate equations in all domain

 implicit none
 include '../var.f90'
 include 'comm.mshbase'
 include 'comm.varbase'
 
 integer(4):: i, j
 
  !open (unit=1,file='initdata.out',form='unformatted')
  open (unit=1,file='initdata.out',status='unknown')
   !write(1) x,y  
   60 format('variables = "x","y","u","v","rho","T"')
   61 format(' zone i=',i4,' j=',i4, ' f=point')
   59 format(6(x,e14.7))
   write(1,60)
   write(1,61) jmax, imax

   do i = 1, imax
     !write(*,*) '.......integrating for x(',i,')=',x(i)
     ! calculating eta
     call calc_eta(x(i))
     
     ! integrating positive side
     call int_plus(fpp_f,gp_f,fp_f,g_f)
     
     ! integrating negative side
     call int_minus(fpp_f,gp_f,fp_f,g_f)
     
     ! calculating flow variables
     call calc_flowvar
	 
	 !write(1) i, u, v, rho, T
	 do j= 1, jmax
	   write(1,59) x(i), y(j), u(j), v(j), rho(j), T(j)
     enddo
	 
   enddo
  close (1)


return
end subroutine integra

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
 real(8):: visc, densi, ksii

 ! u: x velocity
 u = solu(2,:)

 ! T: temperature
 T = solu(4,:)

 ! rho: density
 rho = 1.d0/T
 
  ! v: y velocity
  visc  = 1.785d-5
  densi = 1.3519d0
  ksii = densi*Umax1*visc*x(1)
  
  !v = (visc/dsqrt(2.d0*ksii))*T*(eta*solu(2,:) - solu(1,:))
  
  v = 0.5d0*dsqrt(Umax1*mi/x(1))*(eta*solu(2,:) - solu(1,:))
  !v = (1.d0/Umax1)*(0.5d0*dsqrt(Umax1*mi/x(1))*(eta*solu(2,:) - solu(1,:)))
  !v = (1.d0/dsqrt(2.d0*Re))*(eta*solu(2,:) - solu(1,:))

 !open (1, file='results.out',status='unknown')
 !  write(1,*) 'j u v rho T y eta'
 !  do j = 1, jmax
 !    write(1,*) j, u(j), v(j), rho(j), T(j), y(j), eta(j)
 !  enddo
 !close(1)

return
end subroutine calc_flowvar

!*******************************************************************************
subroutine calc_vorticitythickness

 ! calculates non-dimnesional vorticity thickness
 ! delta_vorti = delta_fp / |dfp/dy|max
 
 implicit none
 include '../var.f90'
 include 'comm.mshbase'
 include 'comm.varbase'

 integer(4):: j
 real(8):: deltafpmax, deltafp, dfpdy, delta_vorti

 ! Calculating delta_fp
 deltafp = solu(2,jmax)-solu(2,1)
 
 ! Calculates |dfp/dy|max
 deltafpmax = 0.d0

 do j = 2, jmax
   dfpdy = (solu(2,j) - solu(2,j-1))/(y(j)-y(j-1))
   deltafpmax = dmax1(deltafpmax,dfpdy)
 enddo
 
 ! calculates delta_vorti
 delta_vorti = deltafp / deltafpmax
 
 write(*,*) 'vorticity thickness = ', delta_vorti

return
end subroutine calc_vorticitythickness

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
