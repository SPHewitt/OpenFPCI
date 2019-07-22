  SUBROUTINE runbd(node,val,num_var,mat_prop,nod,nr,loaded_nodes,timeStep, 	&
                      g_g_pp,g_num_pp,g_coord_pp,gravlo_pp,Dfield,Ufield,Afield,flag)
  !/****f* parafembd/runbd
  !*  NAME
  !*    SUBROUTINE: runbd
  !*
  !*  SYNOPSIS
  !*    Usage:     runbd_(forceNodes_,fext_OF_,numSchemes_,solidProps_,
  !*                     &numRestrNodes_,&numFixedForceNodes_, &dtim,
  !*                     g_g_pp_OF_,g_num_pp_OF_,g_coord_pp_OF_,gravlo_,
  !*                     ptDtemp_,ptUtemp_,ptAtemp_);
  !*
  !*  FUNCTION
  !*    Reads in the current timesteps displacement, velocity,
  !*    acceleration and external force field. Loads the structure
  !*    and solves the governing equations of the problem
  !*
  !*    Backward difference scheme:
  !*
  !*    Liu TY, Li Q Bin, Zhao C Bin. An efficient time-integration method for
  !*    nonlinear dynamic analysis of solids and structures.
  !*    Sci China Physics, Mech Astron. 2013;56(4):798–804.
  !*
  !*
  !*  INPUTS
  !*    val         (ndim,loaded_nodes)	 - Force vector of loaded nodes
  !*    node        (loaded_nodes)       - Node # of loaded_nodes
  !*    num_var     (a1 b1 theta dTim)   - Numerical Variables
  !*    mat_prop    (e v rho)            - Material Properties
  !*
  !*    nr                               - Number of restrained nodes
  !*    loaded_nodes                     - Number of loaded nodes
  !*    time_step                        - time step
  !*
  !*    g_g_pp      (ntot,nels_pp)       - Distributed equation steering array
  !*    g_num_pp    (nod,nels_pp)        - Distributed element steering array
  !*    g_coord_pp  (nod,nels_pp)        - Distributed nodal coordinates
  !*    store_km_pp (ntot,ntot,nels_pp)  - Distributed stiffness matrix [k]
  !*    store_mm_pp (ntot,ntot,nels_pp)  - Distributed mass matrix [M]
  !*
  !*    diag_precon_pp (neq_pp)          - Distributed diagonal preconditioner
  !*    gravlo_pp      (neq_pp)          - Vector of gravity loads
  !*
  !*  OUTPUT
  !*    Dfield  (ntot,nels_pp)           - Distributed nodal displacements
  !*    Ufield  (ntot,nels_pp)           - Distributed nodal velocities
  !*    Afield  (ntot,nels_pp)           - Distributed nodal accelerations
  !*
  !*  AUTHOR
  !*    S. Hewitt
  !******
  !*/

  USE mpi_wrapper;     USE precision;    USE global_variables;
  USE mp_interface;    USE input;        USE output;
  USE loading;         USE timing;       USE maths;
  USE gather_scatter;  USE steering;     USE new_library;
  USE large_strain;

  IMPLICIT NONE

!------------------------------------------------------------------------------
! 1. Declare variables used in the main program
!------------------------------------------------------------------------------

  INTEGER,PARAMETER         :: nodof=3,ndim=3,nst=6
  REAL(iwp),PARAMETER       :: zero=0.0_iwp,one=1.0_iwp

  INTEGER,INTENT(IN)        :: loaded_nodes,nr,nod

  INTEGER,INTENT(INOUT)     :: g_g_pp(ntot,nels_pp),flag
  INTEGER,INTENT(INOUT)     :: g_num_pp(nod,nels_pp),node(loaded_nodes)


  INTEGER                   :: printres
  INTEGER                   :: nels,nn,nip,nlen
  INTEGER                   :: nf_start, fmt=1, i, j, k, m
  INTEGER                   :: iters, limit, iel
  INTEGER                   :: num_load_steps, iload, igauss
  INTEGER                   :: dimH, inewton, jump, npes_pp
  INTEGER                   :: partitioner=1
  INTEGER                   :: nodes_pp, node_start
  INTEGER                   :: node_end, idx1, idx2

  REAL(iwp),INTENT(IN)      :: num_var(8),mat_prop(3),timeStep,g_coord_pp(nod,ndim,nels_pp)

  REAL(iwp),INTENT(INOUT)   :: gravlo_pp(neq_pp)
  REAL(iwp),INTENT(INOUT)   :: Dfield(ntot,nels_pp),Ufield(ntot,nels_pp)
  REAL(iwp),INTENT(INOUT)   :: Afield(ntot,nels_pp), val(ndim,loaded_nodes)

  REAL(iwp)                 :: up,ray_a,ray_b,radius
  REAL(iwp)                 :: e,v,rho,det,tol, maxdiff, tol2, detF
  REAL(iwp)                 :: energy, energy1, rn0
  REAL(iwp)                 :: a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
  REAL(iwp)                 :: dtim

  REAL(iwp)                 :: xi,eta,zeta,etam,xim,zetam,etap,xip,zetap
  REAL(iwp),SAVE            :: time

  CHARACTER(len=15)         :: element
  CHARACTER(len=50)         :: text
  CHARACTER(len=50)         :: fname_base, fname
  CHARACTER(LEN=50)         :: argv

  LOGICAL :: converged, timewrite=.TRUE.

!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------
  REAL(iwp),SAVE,ALLOCATABLE  :: timest(:),nr_timest(:,:)

  REAL(iwp),SAVE,ALLOCATABLE  :: points(:,:),coord(:,:),weights(:)
  REAL(iwp),SAVE,ALLOCATABLE  :: r_pp(:),xnew_pp(:),bee(:,:)
  REAL(iwp),SAVE,ALLOCATABLE  :: diag_precon_pp(:)
  REAL(iwp),SAVE,ALLOCATABLE  :: diag_precon_tmp(:,:)
  REAL(iwp),SAVE,ALLOCATABLE  :: res_pp(:),fint_pp(:)
  REAL(iwp),SAVE,ALLOCATABLE  :: kmat_elem(:,:), kgeo_elem(:,:)
  REAL(iwp),SAVE,ALLOCATABLE  :: xnewel_pp(:,:), jacF(:,:),auxm(:,:)
  REAL(iwp),SAVE,ALLOCATABLE  :: derivFtran(:,:), derivF(:,:)
  REAL(iwp),SAVE,ALLOCATABLE  :: beeF(:,:), defE(:,:)
  REAL(iwp),SAVE,ALLOCATABLE  :: cmat(:,:,:,:), sigma(:,:), cspa(:,:,:,:)
  REAL(iwp),SAVE,ALLOCATABLE  :: sigma1C(:), storefint_pp(:,:), deeF(:,:)
  REAL(iwp),SAVE,ALLOCATABLE  :: geomH(:,:)
  REAL(iwp),SAVE,ALLOCATABLE  :: piolaS(:,:)
  REAL(iwp),SAVE,ALLOCATABLE  :: fextpiece_pp(:)
  REAL(iwp),SAVE,ALLOCATABLE  :: fext_pp(:), deltax_pp(:)
  REAL(iwp),SAVE,ALLOCATABLE  :: x1_pp(:),d1x1_pp(:),d2x1_pp(:)
  REAL(iwp),SAVE,ALLOCATABLE  :: x0_pp(:),d1x0_pp(:),d2x0_pp(:)
  REAL(iwp),SAVE,ALLOCATABLE  :: vu_pp(:),xu_pp(:)
  REAL(iwp),SAVE,ALLOCATABLE  :: pmul_pp(:,:),utemp_pp(:,:)
  REAL(iwp),SAVE,ALLOCATABLE  :: temp_pp(:,:,:)
  REAL(iwp),SAVE,ALLOCATABLE  :: vel_pp(:),acel_pp(:),eld_pp(:,:)
  REAL(iwp),SAVE,ALLOCATABLE  :: shape_integral_pp(:,:)

  REAL(iwp),SAVE,ALLOCATABLE  :: fun(:),emm(:,:),ecm(:,:)
  REAL(iwp),SAVE,ALLOCATABLE  :: d2x1_ppstar(:),meff_pp(:),ceff_pp(:)
  REAL(iwp),SAVE,ALLOCATABLE  :: storemm_pp(:,:,:),storekm_pp(:,:,:)
  REAL(iwp),SAVE,ALLOCATABLE  :: storecm_pp(:,:,:)


  INTEGER,SAVE,ALLOCATABLE  :: num(:),nr_iters(:,:)
  INTEGER,SAVE,ALLOCATABLE  :: comp(:,:)


!------------------------------------------------------------------------------
! 3. Start Program
!------------------------------------------------------------------------------

  IF(numpe .EQ. 1)WRITE(*,*)"ParaFEM Large Strain Solver: "

  ! Barrier (may not be needed but safe)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

  ! Variables required for writing and partioning
  argv="Case"; nlen=4; partitioner=1

  SELECT CASE(nod)
    CASE(8)
      nip = 8; element = "hexahedron"; dimH = 8
    CASE(4)
      nip = 4; element = "tetrahedron"; dimH = 4
  END SELECT

  num_load_steps = 1         ! Number of load steps
  printres   = 0             ! Write .res file

  ! Set Numerical and Material Values
  ray_a  =  num_var(4)    ! Damping parameter A
  ray_b  =  num_var(5)    ! Damping parameter B
  tol    =  num_var(6)    ! Tolerance of PCG loop
  limit  =  num_var(7)    ! Max number of Interation in PCG
  tol2   =  num_var(8)    ! Tolerance for Newton-Raphson loop

  dtim    =  timeStep

  ! Youndgs Modulus, Poissons Ratio and Density
  e     =  mat_prop(1)
  v     =  mat_prop(2)
  rho   =  mat_prop(3)

  ! Allocate memory required for the time loop
  IF(.NOT.ALLOCATED(timest))THEN
    ! Vectors{eqns}
    ALLOCATE(diag_precon_pp(0:neq_pp))
    ALLOCATE(r_pp(0:neq_pp))
    ALLOCATE(res_pp(0:neq_pp))

    ALLOCATE(fint_pp(0:neq_pp))
    ALLOCATE(fextpiece_pp(0:neq_pp))
    ALLOCATE(fext_pp(0:neq_pp))

    ALLOCATE(xnew_pp(0:neq_pp))
    ALLOCATE(deltax_pp(0:neq_pp))

    ALLOCATE(x0_pp(0:neq_pp))
    ALLOCATE(d1x0_pp(0:neq_pp))
    ALLOCATE(d2x0_pp(0:neq_pp))

    ALLOCATE(d2x1_ppstar(0:neq_pp))

    ALLOCATE(x1_pp(0:neq_pp))
    ALLOCATE(d1x1_pp(0:neq_pp))
    ALLOCATE(d2x1_pp(0:neq_pp))

    ALLOCATE(vu_pp(0:neq_pp))
    ALLOCATE(xu_pp(0:neq_pp))
    ALLOCATE(meff_pp(0:neq_pp))
    ALLOCATE(ceff_pp(0:neq_pp))

    ! Matricies
    ALLOCATE(temp_pp(ntot,ntot,nels_pp))
    ALLOCATE(pmul_pp(ntot,nels_pp))
    ALLOCATE(utemp_pp(ntot,nels_pp))
    ALLOCATE(diag_precon_tmp(ntot,nels_pp))
    ALLOCATE(ecm(ntot,ntot))
    ALLOCATE(emm(ntot,ntot))
    ALLOCATE(fun(nod))
    ALLOCATE(weights(nip))
    ALLOCATE(timest(20))
    ALLOCATE(nr_timest(10,20))
    ALLOCATE(nr_iters(10,1))

    printres=1

   IF(numpe .EQ. 1)WRITE(*,*)"------------------------------ "
   IF(numpe .EQ. 1)WRITE(*,'(A)')"Backward Difference method"
   IF(numpe .EQ. 1)WRITE(*,*)"------------------------------ "
  ENDIF

  !---- Clean Arrays ------

  x0_pp    =  0._iwp;    d1x0_pp =  0._iwp;
  d2x0_pp  =  0._iwp;

  d1x1_pp  =  0._iwp;    d2x1_pp =  0._iwp;
  x1_pp    =  0._iwp;

  vu_pp    =  0._iwp;    xu_pp   =  0._iwp;

  nr_timest = zero;


  IF(.NOT.ALLOCATED(coord))THEN
    ALLOCATE(coord(nod,ndim))
    ALLOCATE(bee(nst,ntot))
    ALLOCATE(num(nod))
    ALLOCATE(storekm_pp(ntot,ntot,nels_pp))
    ALLOCATE(storemm_pp(ntot,ntot,nels_pp))
    ALLOCATE(storecm_pp(ntot,ntot,nels_pp))
    ALLOCATE(kmat_elem(ntot,ntot))
    ALLOCATE(kgeo_elem(ntot,ntot))
    ALLOCATE(xnewel_pp(ntot,nels_pp))
    ALLOCATE(comp(nod,ndim))
    ALLOCATE(jacF(ndim,ndim))
    ALLOCATE(auxm(nod,ndim))
    ALLOCATE(derivFtran(nod,ndim))
    ALLOCATE(derivF(ndim,nod))
    ALLOCATE(beeF(nst,ntot))
    ALLOCATE(defE(ndim,ndim))
    ALLOCATE(piolaS(ndim,ndim))
    ALLOCATE(cmat(ndim,ndim,ndim,ndim))
    ALLOCATE(sigma(ndim,ndim))
    ALLOCATE(cspa(ndim,ndim,ndim,ndim))
    ALLOCATE(sigma1C(nst))
    ALLOCATE(storefint_pp(ntot,nels_pp))
    ALLOCATE(deeF(nst,nst))
    ALLOCATE(geomH(dimH,dimH))
    ALLOCATE(points(ndim,nip))
  ENDIF

  timest        =  zero
  timest(1)     =  elap_time()

!------------------------------------------------------------------------------
! 4. Get integration Gauss points and weights in the element
!------------------------------------------------------------------------------

  CALL GET_GAUSS_POINTS(element,points,weights)

!------------------------------------------------------------------------------
! 5. Set Loads
!------------------------------------------------------------------------------
  timest(2)     =  elap_time()

  CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

  fext_pp = zero
  CALL load(g_g_pp,g_num_pp,node,val,fext_pp(1:))

  CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

  fext_pp(1:) = fext_pp(1:) + gravlo_pp

  fextpiece_pp(1:) = fext_pp(1:)/FLOAT(num_load_steps)

!------------------------------------------------------------------------------
! 6. Set Initial Conditions
!------------------------------------------------------------------------------
  timest(3)     =  elap_time()

  ! Scatter_noadd has no barriers in the subroutine
  CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
  CALL scatter_noadd(Dfield,x0_pp(1:))

  CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
  CALL scatter_noadd(Ufield,d1x0_pp(1:))

  CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
  CALL scatter_noadd(Afield,d2x0_pp(1:))

  CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

!-------------------------------------------------------------------------
! 7. Initialise the solution vector to 0.0
!-------------------------------------------------------------------------
  timest(4)     =  elap_time()

  ! U_n = U
  xnew_pp = x0_pp

  ! Vector comp to compute F (gradient of deformation)
  DO i = 1,nod
    DO j = 1,ndim
      comp(i,j) = (i-1)*ndim + j
    END DO
  END DO

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!----------------------- Start Load LOOP --------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  timest(5)     =  elap_time()

  DO iload = 1,num_load_steps
    converged = .FALSE.

   !fext_pp(1:) = FLOAT(iload)*fextpiece_pp(1:)

!------------------------------------------------------------------------------
!----------------------- Start Newton-Raphson iterations ----------------------
!------------------------------------------------------------------------------
    inewton = 0
    iterations: DO
      inewton = inewton + 1

      timest(7)     =  elap_time()

      storefint_pp = 0._iwp
      xnewel_pp = zero

      CALL GATHER(xnew_pp(1:),xnewel_pp)

      timest(8)     =  elap_time()
      nr_timest(inewton,1)= timest(8)-timest(7)

!-------------------------------------------------------------------------
! 8. Build Matricies (K , M , f_int)
!-------------------------------------------------------------------------
      timest(9)     =  elap_time()

      ! Clean [K], [M] and [C]
      storekm_pp  =  zero
      storemm_pp  =  zero
      storecm_pp  =  zero

      DO iel = 1,nels_pp
        kmat_elem = 0._iwp
        kgeo_elem = 0._iwp
        DO i = 1,nod
          num(i) = g_num_pp(i,iel)
        END DO
        coord = g_coord_pp(:,:,iel)
        auxm(:,1) = xnewel_pp(comp(:,1),iel)
        auxm(:,2) = xnewel_pp(comp(:,2),iel)
        auxm(:,3) = xnewel_pp(comp(:,3),iel)

        emm  	=  0._iwp
        ecm  	=  0._iwp

        DO igauss = 1,nip
          CALL kine3D(igauss,auxm,coord,points,det,detF,beeF,defE,derivF,jacF)
          CALL venantkirchhoff(defE,e,v,piolaS,cmat)
          CALL push2r(detF,jacF,piolaS,sigma)
          CALL push4r(detF,jacF,cmat,cspa)
          CALL voigt2to1(sigma,sigma1C)
          CALL voigt4to2(cspa,deeF)

          ! F_int
          storefint_pp(:,iel) = storefint_pp(:,iel) +         &
             MATMUL(TRANSPOSE(beeF),sigma1C)*det*detF*weights(igauss)
          ! K_l
          kmat_elem(:,:) = kmat_elem(:,:) +     &
             MATMUL(MATMUL(TRANSPOSE(beeF),deeF),beeF)*det*detF*weights(igauss)
          ! K_nl
          geomH = MATMUL(MATMUL(TRANSPOSE(derivF),sigma),derivF)

          DO i = 1,dimH
            DO j = 1,dimH
              kgeo_elem(3*i-2,3*j-2) = kgeo_elem(3*i-2,3*j-2) + &
                  geomH(i,j)*det*detF*weights(igauss)
              kgeo_elem(3*i-1,3*j-1) = kgeo_elem(3*i-1,3*j-1) + &
                  geomH(i,j)*det*detF*weights(igauss)
              kgeo_elem(3*i,3*j) = kgeo_elem(3*i,3*j)         + &
                  geomH(i,j)*det*detF*weights(igauss)
            END DO
          END DO


        ! Mass Matrix
        fun = zero

        ! BUG: shape_Fun in shared/new_library.f90
        ! Line 253
        ! ndim = UBOUND(points,2) -> ndim = UBOUND(points,1)

        !CALL shape_fun(fun,points,igauss)
        ! SHAPE_ FUNCITON ERRORS
        ! fun calculated manually
        ! nod =8, ndim=3
        xi    = points(1,igauss)
        eta   = points(2,igauss)
        zeta  = points(3,igauss)
        etam  = one  - eta
        xim   = one  - xi
        zetam = one  - zeta
        etap  = eta  + one
        xip   = xi   + one
        zetap = zeta + one

       fun = (/0.125_iwp*xim*etam*zetam,0.125_iwp*xim*etam*zetap,          &
               0.125_iwp*xip*etam*zetap,0.125_iwp*xip*etam*zetam,          &
               0.125_iwp*xim*etap*zetam,0.125_iwp*xim*etap*zetap,          &
               0.125_iwp*xip*etap*zetap,0.125_iwp*xip*etap*zetam/)

        CALL ecmat(ecm,fun,ntot,nodof)
        ecm  =  ecm*det*weights(igauss)*rho
        emm  =  emm+ecm

        END DO ! Gauss Points

        ! k = k_l + k_nl
        storekm_pp(:,:,iel) = kmat_elem(:,:) + kgeo_elem(:,:)

        ! M
        storemm_pp(:,:,iel)  =  emm

        ! C
        storecm_pp(:,:,iel)=ray_a*storemm_pp(:,:,iel)+ray_b*storekm_pp(:,:,iel)

      END DO ! nels_pp

      timest(10)     =  elap_time()
      nr_timest(inewton,2)= timest(10)-timest(9)

      ! F_int
      fint_pp(:) = .0_iwp
      CALL SCATTER(fint_pp(1:),storefint_pp)

      timest(11) = elap_time()
      nr_timest(inewton,3)= timest(11)-timest(10)

!-------------------------------------------------------------------------
! 9. Backward difference scheme
!-------------------------------------------------------------------------
     a0  = (dtim*dtim)/(2.0)
     a1  = (dtim)/(2.0)

     ! M_eff
     meff_pp = zero
     meff_pp(1:) = x0_pp(1:)-xnew_pp(1:) + dtim*d1x0_pp(1:)

     temp_pp =  zero
     temp_pp = storemm_pp
     pmul_pp = zero

     ! M*M_eff
     CALL GATHER(meff_pp(1:),pmul_pp) ; utemp_pp=zero

     DO iel=1,nels_pp
       CALL DGEMV('N',ntot,ntot,one,temp_pp(:,:,iel),ntot,                 &
                pmul_pp(:,iel),1,zero,utemp_pp(:,iel),1)
     END DO

     vu_pp = zero
     CALL SCATTER(vu_pp(1:),utemp_pp)

      !C_eff
     !ceff_pp = zero
     !ceff_pp(1:) = a1*(x0_pp(1:)-xnew_pp(1:)) + a4*d1x0_pp(1:) + a5*d2x0_pp(1:)

     !temp_pp  =  zero
     !temp_pp  =  storecm_pp
     !pmul_pp = zero

     ! C*C_eff
     !CALL GATHER(ceff_pp(1:),pmul_pp) ; utemp_pp=zero

     !DO iel=1,nels_pp
     !  CALL DGEMV('N',ntot,ntot,one,temp_pp(:,:,iel),ntot,                 &
     !           pmul_pp(:,iel),1,zero,utemp_pp(:,iel),1)
     !END DO

     !xu_pp = zero
     !CALL SCATTER(xu_pp(1:),utemp_pp)

!-------------------------------------------------------------------------
! 10. Get residual
!-------------------------------------------------------------------------

     ! {r_pp}
     r_pp(1:) = a0*(fext_pp(1:) - fint_pp(1:)) + vu_pp(1:)

     ! Compute maxdiff of residual
     maxdiff =  MAXABSVAL_P(r_pp(1:))

     ! Normalise residual vector and stiffness matrix for pcg
     IF (maxdiff == 0.0) THEN
       EXIT
     END IF

     ! [k]
     storekm_pp =  a0*storekm_pp  + a1*storecm_pp + storemm_pp

     timest(12)     =  elap_time()
     nr_timest(inewton,4)= timest(12)-timest(11)

!-------------------------------------------------------------------------
! 11. Diagonal Preconditioner
!-------------------------------------------------------------------------
      diag_precon_tmp = .0_iwp
      DO iel = 1,nels_pp
        DO k = 1,ntot
          diag_precon_tmp(k,iel)=diag_precon_tmp(k,iel) + storekm_pp(k,k,iel)
        END DO
      END DO

      diag_precon_pp(:) = .0_iwp
      CALL SCATTER(diag_precon_pp(1:),diag_precon_tmp)


      diag_precon_pp(1:) = 1._iwp/diag_precon_pp(1:)
      diag_precon_pp(0)  = .0_iwp

      timest(13)     =  elap_time()
      nr_timest(inewton,5)= timest(13)-timest(12)
!---------------------------------------------------------------------------
!------------------------------- Solve using PCG ---------------------------
!---------------------------------------------------------------------------
      deltax_pp = .0_iwp
      res_pp    = r_pp

      ! Include dynamic Effects
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

      CALL PCG_VER1(inewton,limit,tol,storekm_pp,r_pp(1:), &
                    diag_precon_pp(1:),rn0,deltax_pp(1:),iters)

      nr_iters(inewton,1)=iters

      ! To ensure completion of PCG
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

      !IF(numpe .EQ. 1)WRITE(*,'(2(a,I3))'),"N-R: ",inewton," PCG iters: ",iters

      xnew_pp(1:) = xnew_pp(1:) + deltax_pp(1:)
      xnew_pp(0) = .0_iwp

      timest(14)     =  elap_time()
      nr_timest(inewton,6)= timest(14)-timest(13)

!-------------------------------------------------------------------------
! 12. Check convergence for Newton-Raphson iterations
!-------------------------------------------------------------------------

      energy = ABS(DOT_PRODUCT_P(res_pp(1:),deltax_pp(1:)))
      IF (inewton==1) THEN
       energy1 = energy
      END IF

      !if(numpe .EQ. 1)WRITE(*,*)iload,inewton,energy,energy/energy1

      IF (inewton>1) THEN
        IF ((energy/energy1)<=tol2) THEN
          converged = .TRUE.
        END IF
      END IF

      IF(converged .OR. inewton==10) THEN
        EXIT
      END IF

    END DO iterations

   IF(numpe .EQ. 1)WRITE(*,'(a,I3,a,ES10.3)')," Newton-Raphson Iters: ",inewton,&
                                        ",  Final residual: ", (energy/energy1)

   nr_timest(inewton,7)= elap_time()-timest(14)

  END DO !iload

!-------------------------------------------------------------------------
! 13. Update Velocity and Acceleration
!-------------------------------------------------------------------------

   a0  = (dtim*dtim)/(2.0)
   a1  = (dtim)/(2.0)
   a3  = (2.0)/(dtim*dtim)

   timest(15)     =  elap_time()

   x1_pp=zero; d2x1_ppstar= zero; d1x1_pp= zero; d2x1_pp=zero;

   x1_pp(1:) = xnew_pp(1:)

   d2x1_ppstar(1:) = a3*(x1_pp(1:)-x0_pp(1:) - dtim*d1x0_pp(1:))
   d1x1_pp(1:)     = d1x0_pp(1:) + a1*(d2x0_pp(1:) + d2x1_ppstar(1:))
   d2x1_pp(1:)     = d2x1_ppstar(1:)

!------------------------------------------------------------------------------
! 14. Gather Data from ntot,nels_pp to ndim,nodes_pp
!------------------------------------------------------------------------------
   timest(16)     =  elap_time()

  IF(.NOT.ALLOCATED(eld_pp))THEN
    ALLOCATE(eld_pp(ntot,nels_pp))
    printres=1
  END IF

  ! Displacement
  eld_pp   =  zero
  CALL gather(x1_pp(1:),eld_pp)
  Dfield=eld_pp

  ! Velocity
  eld_pp   =  zero
  CALL gather(d1x1_pp(1:),eld_pp)
  Ufield=eld_pp

  ! Acceleration
  eld_pp   =  zero
  CALL gather(d2x1_pp(1:),eld_pp)
  Afield=eld_pp

  CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

  timest(17)     =  elap_time()

  ! Write at Information Timestep XX
  IF(printres == 1)THEN
    CALL WRITE_LARGESTRAIN(argv,nn,nr,loaded_nodes,timest,nr_timest,inewton,nr_iters)
  ENDIF

  END SUBROUTINE

  !------------------------------------------------------------------------
  !------------------------------------------------------------------------
  !------------------------------------------------------------------------
