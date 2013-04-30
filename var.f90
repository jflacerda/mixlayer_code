integer(4):: nprocx, nprocy, lvc, mp, ncomu, ivar, inter, interderv, interfilt, imaxin, jmaxin, ptsx, ptsy, &
             imax, jmax, iteramax, nmodes, npdamp, npfilter, icomp

logical:: restart, ldisturb, lprintform, lfilter

real(8):: dt, tfinal
real(8):: dksi, dn, stx, sty, xl, yl, posxin, posyin
real(8):: mach, Re, Pr, cv, kappa, kcond, mi
real(8):: Umax1, Umax2, Uratio, Tmax1, Hemax1, Hemax2, Heratio
real(8):: af, bf, cf, df, alphaf
real(8):: omega, t_period
real(8):: damp, pi
real(8):: fp_guess, g_guess, fpp_guess, gp_guess


parameter (restart  = .false.)   ! se reinicia ou não a simulação a partir de uma anterior
parameter (lprintform = .true.) ! se imprime formatado ou nao

parameter (nprocx = 1)  ! numero de sub-dominios na direcao x
parameter (nprocy = 1)  ! numero de sub-dominios na direcao y
parameter (lvc = 1)     ! numero de niveis do ciclo v
parameter (mp = 5)     ! numero de pontos da molecula computacional necessaria para 
                       ! o cálculo das derivadas
parameter (ncomu =  3) ! numero de cominicacoes na volta do pipeline
parameter (ivar  =  4) ! numero de dimensoes da variavel Q
parameter (icomp  = 2500) ! intervalo de passos de tempo para fazer comparação

!Definindo o valor de inter (numero de volumes no overlap de cada malha)
parameter(interderv = 2**(lvc-1)*(mp-2))
parameter(interfilt = 2) ! o +2 ocorre porque o filtro precisa de 2 pontos adicionais
parameter(inter = interderv + interfilt)

!parameter (imaxin = 61)  
!parameter (jmaxin = imaxin)
!parameter (imaxin = 2500, jmaxin=850)  ! Malha Quali
parameter (imaxin = 200, jmaxin=200)  ! Malha regular em y - comparacao Salemi

! Definindo numero de pontos para cada nó
parameter (ptsx = (imaxin+(inter+1)*(nprocx-1))/nprocx)
parameter (ptsy = (jmaxin+(inter+1)*(nprocy-1))/nprocy)

!Total de número de pontos em cada direcao somando todos os nós
parameter (imax = ptsx*nprocx -(inter+1)*(nprocx-1), jmax = ptsy*nprocy -(inter+1)*(nprocy-1))

parameter (dksi     = 0.157d0)              !tamanho do delta na direcao ksi
parameter (dn       = 0.149998305009144d0)  !tamanho do delta na direcao n

! Para MMS ------------------------------------------------------------------
!parameter (pi       = dacos(-1.d0))
!parameter (xl       =  2.d0*pi )
!parameter (yl       =  2.d0*pi )
!parameter (posxin   =  -pi + pi/16.d0)
!parameter (posyin   =  -pi + pi/16.d0 )
!parameter (dksi     = xl/dfloat(imax-1))  !tamanho do delta na direcao ksi
!parameter (dn       = yl/dfloat(jmax-1))  !tamanho do delta na direcao n
!!modificado para o teste de ordem com stretching
!parameter (stx      = 1.d0)
!parameter (sty      = 1.d0)
!----------------------------------------------------------------------------

parameter (mach     = 0.1d0)
parameter (Re       = 250.d0)
parameter (Pr       = 0.71d0)
parameter (Umax1    = 171.6d0)
parameter (Umax2    = 85.8d0)
parameter (Uratio   = Umax2/Umax1)
parameter (Tmax1    = 280.d0)
parameter (Hemax1   = 280.1d3)
parameter (Hemax2   = 280.1d3)
parameter (Heratio  = Hemax1/Hemax2)

parameter (cv       = 0.1d0)
parameter (mi       = 13.204d-6)
parameter (kcond    = 1.d0)
parameter (kappa    = 1.4d0 )
parameter (dt       = 1.d-10 )
parameter (tfinal   = 5.d-9 )
parameter (iteramax = dint(tfinal/dt)+1)

!guesses for baseflow
parameter (fp_guess   = 0.78d0)
parameter (g_guess   = -1.d0)
parameter (fpp_guess   = 0.025d0)
parameter (gp_guess   = 7.d-2)


!constantes para MMS
!real(8):: rho0, u0, v0, et0, omg, eps
!parameter(rho0 = 1.d0)
!parameter(  u0 = 0.1d0)
!parameter(  v0 = 1.0d0)
!parameter( et0 = 0.5d0)
!parameter( omg = 1.0d0)
!parameter( eps = 0.5d0)

!Parametros para o filtro
!Filter constants (Lele C.2.5)
parameter ( lfilter = .false. )
parameter ( npfilter = 5)
parameter ( alphaf = 0.48d0 )
parameter ( af = (11.d0+10.d0*alphaf)/16.d0 )
parameter ( bf = (15.d0+34.d0*alphaf)/64.d0 ) !/2
parameter ( cf = (-3.d0+ 6.d0*alphaf)/32.d0 ) !/2
parameter ( df = ( 1.d0- 2.d0*alphaf)/64.d0 ) !/2

! Parametros para o damping na direcao y
parameter ( npdamp = 40    ) ! numero de pontos no amortecimento
parameter ( damp   = 1.d-2 ) ! fator de amortecimento

!Variáveis para excitacao do LST
parameter (ldisturb = .false.)  ! se calcula ou nao as perturbaçoes dadas pela LST e as características
parameter (nmodes   = 4)        ! número de modos que serao excitados
parameter (omega    = 1.d0)     ! freq de excitacao
parameter (t_period = 1.d0)     ! period
