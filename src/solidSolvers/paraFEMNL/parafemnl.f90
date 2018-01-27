
  SUBROUTINE initnl(g_coord,rest,nn,nr,g_num_pp,g_g_pp,nn_pp,nn_start)

  USE mpi_wrapper;    USE precision;  USE global_variables; 
  USE mp_interface;   USE input;      USE output; 
  USE loading;        USE timing;     USE maths; 
  USE gather_scatter; USE steering;   USE new_library;
  USE large_strain;

  IMPLICIT NONE

!----------------------------------------------------------------------
! 1. Declare variables
!----------------------------------------------------------------------


  INTEGER,PARAMETER         :: nodof=3,ndim=3,nst=6,nod=8
  REAL(iwp),PARAMETER       :: zero=0.0_iwp

  INTEGER,INTENT(IN) 	    :: nn,nr

  INTEGER,INTENT(INOUT)     :: g_g_pp(ndim*nod,nels_pp),rest(nr,nodof+1)
  INTEGER,INTENT(INOUT)     :: g_num_pp(nod,nels_pp),nn_pp,nn_start

  INTEGER                   :: iel,nip,npes_pp,partitioner,ndof
  INTEGER                   :: nlen,i

  REAL(iwp),INTENT(IN)      :: g_coord(ndim,nn)

  CHARACTER(LEN=50)         :: argv
  CHARACTER(LEN=15)         :: element

  LOGICAL                   :: initialised
 
!----------------------------------------------------------------------
! 3. Input and Initialisation
!---------------------------------------------------------------------- 

  ! Degrees of Freedon per Element
  ndof  =  nod*nodof
  ntot  =  ndof

  ! Variables required for writing and partioning
  argv="Case"; nlen=4; nip=8; element="hexahedron"; partitioner=1

  ! Test for MPI Initialisation
  CALL MPI_INITIALIZED(initialised,ier) 
  
  IF(initialised .EQV. .false.)THEN
    CALL find_pe_procs(numpe,npes)
    numpe  =  1
    npes   =  1
  ELSE
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,numpe,ier)
    ! C++ starts at 0, fortran 1
    numpe  =  numpe+1
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,npes,ier); 
  ENDIF
  
  ! Calculate number of Elements per Core
  ! CALL calc_nels_pp(argv,nels,npes,numpe,partitioner,nels_pp)

  ! nels_pp calculated in Foam-Extend
  ! Calculate iel_start
  CALL setielstart() 
 
!----------------------------------------------------------------------
! 4. Populate g_num_pp
!----------------------------------------------------------------------
  
  ! Poulate the Steering matrix
  !CALL POPULATE_G_NUM_PP(g_num,g_num_pp,npes,nod,nels)

  ! Convert from Foam-Extend to Smith Gritths format
  DO iel=1,nels_pp
    CALL of2sg(element,g_num_pp(:,iel),nod)
  ENDDO

   
!----------------------------------------------------------------------
! 5. Calculate g_g_pp
!----------------------------------------------------------------------

  ! Rearrange the rest Array
  CALL rearrange(rest)

  ! Clean arrays
  g_g_pp  =  zero
  neq	  =  zero
  
  ! Find the global Steering Matrix
  elements_0: DO iel=1,nels_pp
   CALL find_g(g_num_pp(:,iel),g_g_pp(:,iel),rest) ! Stable but slow
  !  CALL find_g3(g_num_pp(:,iel),g_g_pp(:,iel),rest)! Unstable but Fast
  END DO elements_0

  
  ! Build GGL Array
  neq  =  MAXVAL(g_g_pp)
  neq  =  max_p(neq)
  CALL calc_neq_pp
  CALL calc_npes_pp(npes,npes_pp)
  CALL make_ggl(npes_pp,npes,g_g_pp)

  CALL CALC_NN_PP(g_num_pp,nn_pp,nn_start)

  ! output g_g_pp,g_num_pp
 
  END SUBROUTINE

  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
  !--------------------------------------------------------------------

  SUBROUTINE pop_gcoordpp(nn_pp,nn_start,g_coord,g_coord_pp,g_num_pp,nn)

  USE mpi_wrapper;    USE precision;  USE global_variables; 
  USE mp_interface;   USE input;      USE output; 
  USE loading;        USE timing;     USE maths; 
  USE gather_scatter; USE steering;   USE new_library;
  USE large_strain;

  IMPLICIT NONE

  INTEGER,PARAMETER         :: ndim=3, nod=8


  INTEGER,INTENT(IN) 	    :: nn
  INTEGER,INTENT(INOUT)     :: g_num_pp(nod,nels_pp),nn_pp,nn_start

  REAL(iwp),INTENT(INOUT)   :: g_coord_pp(nod,ndim,nels_pp)
  REAL(iwp),INTENT(IN)      :: g_coord(ndim,nn)

  INTEGER                   :: i, k, inpe

   ! Routine hacked from READ_NODES in input.f90
  
    bufsize = ndim*nn

    !CALL MPI_BCAST(g_coord,bufsize,MPI_REAL8,0,MPI_COMM_WORLD,ier)

    nn_pp = UBOUND(g_coord_pp,2) 
    inpe  = nn_start 

    !DO i = 1,nn_pp
    !  g_coord_pp(:,i) = g_coord(:,inpe)
    !  inpe = inpe + 1
    !END DO

  CALL POPULATE_G_COORD_PP2(g_coord,g_coord_pp,g_num_pp,nn,nod,ndim) 

  END SUBROUTINE

  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
  !--------------------------------------------------------------------

  SUBROUTINE runnl(node,val,num_var,mat_prop,nr,loaded_nodes,timeStep,nn_pp,nn_start, 	&
                      g_g_pp,g_num_pp,g_coord_pp,gravlo,Dfield,Ufield,Afield)
  
  USE mpi_wrapper;     USE precision;    USE global_variables; 
  USE mp_interface;    USE input;        USE output; 
  USE loading;         USE timing;       USE maths; 
  USE gather_scatter;  USE steering;     USE new_library; 
  USE large_strain;
  
  IMPLICIT NONE
  
!------------------------------------------------------------------------------ 
! 1. Declare variables used in the main program
!------------------------------------------------------------------------------

  INTEGER,PARAMETER         :: nodof=3,ndim=3,nst=6,nod=8
  REAL(iwp),PARAMETER       :: zero=0.0_iwp,one=1.0_iwp

  INTEGER,INTENT(IN)        :: loaded_nodes,nr

  INTEGER,INTENT(INOUT)     :: g_g_pp(ntot,nels_pp)
  INTEGER,INTENT(INOUT)     :: g_num_pp(nod,nels_pp),node(loaded_nodes)


  INTEGER                   :: nels,nn,nip, nn_pp,nlen
  INTEGER                   :: nf_start, fmt=1, i, j, k, m
  INTEGER                   :: iters, limit, iel, nn_start
  INTEGER                   :: num_load_steps, iload, igauss
  INTEGER                   :: dimH, inewton, jump, npes_pp
  INTEGER                   :: partitioner=1
  INTEGER                   :: nodes_pp, node_start
  INTEGER                   :: node_end, idx1, idx2

  REAL(iwp),INTENT(IN)      :: num_var(4),mat_prop(3),timeStep,g_coord_pp(nod,ndim,nels_pp)

  REAL(iwp),INTENT(INOUT)   :: gravlo(neq_pp)
  REAL(iwp),INTENT(INOUT)   :: Dfield(ntot,nels_pp),Ufield(ntot,nels_pp)
  REAL(iwp),INTENT(INOUT)   :: Afield(ntot,nels_pp), val(ndim,loaded_nodes)

  REAL(iwp)                 :: up,alpha1,beta1,theta
  REAL(iwp)                 :: e,v,rho,det,tol, maxdiff, tol2, detF
  REAL(iwp)                 :: energy, energy1, rn0
  REAL(iwp)                 :: c1,c2,c3,c4
  REAL(iwp)                 :: a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
  REAL(iwp)                 :: dtim,beta,delta,alpha

  REAL(iwp)                 :: xi,eta,zeta,etam,xim,zetam,etap,xip,zetap

  CHARACTER(len=15)         :: element
  CHARACTER(len=50)         :: text
  CHARACTER(len=50)         :: fname_base, fname
  CHARACTER(LEN=50)         :: argv

  INTEGER :: argc, iargc
  INTEGER :: fixed_nodes, numfix_pp, fixdim, writetimes=0
 
  LOGICAL :: converged, timewrite=.TRUE.

!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------
  REAL(iwp),SAVE,ALLOCATABLE  :: timest(:)

  REAL(iwp),SAVE,ALLOCATABLE  :: points(:,:),coord(:,:),weights(:)
  REAL(iwp),SAVE,ALLOCATABLE  :: r_pp(:),xnew_pp(:),bee(:,:)
  REAL(iwp),SAVE,ALLOCATABLE  :: diag_precon_pp(:),load_value(:,:)
  REAL(iwp),SAVE,ALLOCATABLE  :: diag_precon_tmp(:,:), storekm_pp(:,:,:)
  REAL(iwp),SAVE,ALLOCATABLE  :: disp(:,:),g_coord(:,:),val_pp(:)
  REAL(iwp),SAVE,ALLOCATABLE  :: disp_pp(:,:), res_pp(:),fint_pp(:)
  REAL(iwp),SAVE,ALLOCATABLE  :: kmat_elem(:,:), kgeo_elem(:,:)
  REAL(iwp),SAVE,ALLOCATABLE  :: xnewel_pp(:,:), jacF(:,:),auxm(:,:)
  REAL(iwp),SAVE,ALLOCATABLE  :: derivFtran(:,:), derivF(:,:)
  REAL(iwp),SAVE,ALLOCATABLE  :: beeF(:,:), rightCG(:,:), defE(:,:)
  REAL(iwp),SAVE,ALLOCATABLE  :: cmat(:,:,:,:), sigma(:,:), cspa(:,:,:,:)
  REAL(iwp),SAVE,ALLOCATABLE  :: sigma1C(:), storefint_pp(:,:), deeF(:,:)
  REAL(iwp),SAVE,ALLOCATABLE  :: geomH(:,:), fixed_value(:), fixval_pp(:)
  REAL(iwp),SAVE,ALLOCATABLE  :: jacFinv(:,:), piolaS(:,:),fixvalpiece_pp(:)
  REAL(iwp),SAVE,ALLOCATABLE  :: elemdisp(:), fextpiece_pp(:)
  REAL(iwp),SAVE,ALLOCATABLE  :: fext_pp(:), deltax_pp(:)
  REAL(iwp),SAVE,ALLOCATABLE  :: value_shape(:),xnewnodes_pp(:)
  REAL(iwp),SAVE,ALLOCATABLE  :: stress_integral_pp(:,:), stressnodes_pp(:)
  REAL(iwp),SAVE,ALLOCATABLE  :: principal_integral_pp(:,:), princinodes_pp(:)
  REAL(iwp),SAVE,ALLOCATABLE  :: principal(:), reacnodes_pp(:)
  REAL(iwp),SAVE,ALLOCATABLE  :: loads_pp(:),x1_pp(:),d1x1_pp(:)
  REAL(iwp),SAVE,ALLOCATABLE  :: d2x1_pp(:),x0_pp(:),d1x0_pp(:),d2x0_pp(:)
  REAL(iwp),SAVE,ALLOCATABLE  :: vu_pp(:),p_pp(:),d_pp(:),x_pp(:)
  REAL(iwp),SAVE,ALLOCATABLE  :: pmul_pp(:,:),utemp_pp(:,:)
  REAL(iwp),SAVE,ALLOCATABLE  :: temp(:),d1temp(:),d2temp(:),temp_pp(:,:,:)
  REAL(iwp),SAVE,ALLOCATABLE  :: vel_pp(:),acel_pp(:),eld_pp(:,:)
  REAL(iwp),SAVE,ALLOCATABLE  :: fext_o_pp(:),shape_integral_pp(:,:)

  REAL(iwp),SAVE,ALLOCATABLE  :: fun(:),emm(:,:),ecm(:,:)
  REAL(iwp),SAVE,ALLOCATABLE  :: d2x1_ppstar(:),meff_pp(:)
  REAL(iwp),SAVE,ALLOCATABLE  :: storemm_pp(:,:,:),points_r(:,:)

 
  INTEGER,SAVE,ALLOCATABLE  :: num(:),load_node(:),nf_pp(:,:),no_pp(:)             
  INTEGER,SAVE,ALLOCATABLE  :: comp(:,:),fixed_node(:), fixed_dof(:)
  INTEGER,SAVE,ALLOCATABLE  :: fixelem_pp(:), fixdof_pp(:)

 
!------------------------------------------------------------------------------
! 3. Start Program
!------------------------------------------------------------------------------

  IF(numpe .EQ. 1)PRINT*,"ParaFEM: "
  
  ! Barrier (may not be needed but safe)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
  
  ! Set Base paramenter
  argv       =  "Case"       ! Name files write to
  nlen       =  4            ! Length of Name
  nip        =  8            ! Number of Integration Points
  limit      =  1000         ! Max number of Interation in PCG
  tol        =  1.0e-8       ! Tolerance of PCG loop
  element    =  "hexahedron" ! Element Name

  num_load_steps = 1         ! Number of load steps		
  fixed_nodes    = 0         ! Number of restrained degrees of freedom
                             ! with a non-zero applied value

  tol2 = 1.0e-16             ! Tolerance for Newton-Raphson loop

  dimH = 8		     ! NOT SURE WHAT THIS VARIABLE DOES
  
  ! Set Numerical and Material Values 
  alpha1  =  num_var(1)
  beta1   =  num_var(2)
  theta   =  num_var(3)
  dtim    =  timeStep
  !dtim	  =  num_var(4)  

  !IF(numpe .EQ. 1)PRINT*,"Time Step: ",dtim

  c1	    =  (1._iwp-theta)*dtim
  c2	    =  beta1-c1
  c3	    =  alpha1+1._iwp/(theta*dtim);
  c4	    =  beta1+theta*dtim

  ! Youndgs ModulusPoissons Ratio and Density

  e     =  mat_prop(1)
  v     =  mat_prop(2)
  rho   =  mat_prop(3)
 
  ! Allocate memory required for the time loop
  IF(.NOT.ALLOCATED(timest))THEN

    IF(numpe .EQ. 1)PRINT*,"ALLOCATE : 1"

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
    ALLOCATE(d_pp(0:neq_pp))
    ALLOCATE(p_pp(0:neq_pp))
    ALLOCATE(x_pp(0:neq_pp))
    ALLOCATE(meff_pp(0:neq_pp))

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

   !ALLOCATE(loads_pp(0:neq_pp))
   !ALLOCATE(fext_o_pp(neq_pp))
   !fext_o_pp  =  zero
  ENDIF

  !---- Clean Arrays ------

  x0_pp    =  0._iwp;    d1x0_pp =  0._iwp;
  vu_pp    =  0._iwp;    d2x0_pp =  0._iwp; 
  d1x1_pp  =  0._iwp;    d2x1_pp =  0._iwp;
  d_pp	   =  0._iwp;    p_pp    =  0._iwp;    
  x_pp     =  0._iwp;    x1_pp   =  0._iwp;  

      
  IF(.NOT.ALLOCATED(coord))THEN

    IF(numpe .EQ. 1)PRINT*,"ALLOCATE : 2"

    ALLOCATE(coord(nod,ndim))
    ALLOCATE(bee(nst,ntot))
    ALLOCATE(num(nod))
    !ALLOCATE(load_value(ndim,loaded_nodes))
    !ALLOCATE(load_node(loaded_nodes))
    ALLOCATE(nf_pp(nodof,nn_pp))
    ALLOCATE(storekm_pp(ntot,ntot,nels_pp))
    ALLOCATE(storemm_pp(ntot,ntot,nels_pp))
    ALLOCATE(kmat_elem(ntot,ntot))
    ALLOCATE(kgeo_elem(ntot,ntot))
    ALLOCATE(xnewel_pp(ntot,nels_pp))
    ALLOCATE(comp(nod,ndim))
    ALLOCATE(jacF(ndim,ndim))
    ALLOCATE(auxm(nod,ndim))
    ALLOCATE(jacFinv(ndim,ndim))
    ALLOCATE(derivFtran(nod,ndim))
    ALLOCATE(derivF(ndim,nod))
    ALLOCATE(beeF(nst,ntot))
    ALLOCATE(rightCG(ndim,ndim))
    ALLOCATE(defE(ndim,ndim))
    ALLOCATE(piolaS(ndim,ndim))
    ALLOCATE(cmat(ndim,ndim,ndim,ndim))
    ALLOCATE(sigma(ndim,ndim))
    ALLOCATE(cspa(ndim,ndim,ndim,ndim))
    ALLOCATE(sigma1C(nst))
    ALLOCATE(storefint_pp(ntot,nels_pp))
    ALLOCATE(deeF(nst,nst))
    ALLOCATE(geomH(dimH,dimH))
    ALLOCATE(principal(ndim))
    ALLOCATE(elemdisp(ntot))
    ALLOCATE(value_shape(nod))
    ALLOCATE(points(ndim,nip))
    ALLOCATE(points_r(nip,ndim))
  ENDIF

!------------------------------------------------------------------------------
! 1. Get integration Gauss points and weights in the element
!------------------------------------------------------------------------------

  CALL GET_GAUSS_POINTS(element,points,weights)

!------------------------------------------------------------------------------
! 4. Set Loads
!------------------------------------------------------------------------------
   !fext_pp =  0._iwp
   fext_pp = zero

  ! Load fext_pp based on global load vector

  CALL load(g_g_pp,g_num_pp,node,val,fext_pp(1:))
  !CALL load_2((nn_start,g_num_pp,node,val,nf_pp,fext_pp(1:))
 
  CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

  !IF(numpe .EQ. 1)PRINT*,"Gravity Loads: ",SUM_P(gravlo)  
  
  fext_pp(1:) = fext_pp(1:) + gravlo
 
  fextpiece_pp(1:) = fext_pp(1:)/FLOAT(num_load_steps)
  
!------------------------------------------------------------------------------
! 4. Read and distribute essential boundary conditions
!------------------------------------------------------------------------------

  numfix_pp = 0
  ! SEE xx7 for this section
  ! USE of FIXED_NODES 
!------------------------------------------------------------------------------
! 5. Set Initial Conditions
!------------------------------------------------------------------------------

! - scatter_noadd has no barriers in the subroutine
  CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
  CALL scatter_noadd(Dfield,x0_pp(1:))

  CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
  CALL scatter_noadd(Ufield,d1x0_pp(1:))

  CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
  CALL scatter_noadd(Afield,d2x0_pp(1:))

  CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

  !vu_pp      =  0._iwp
  !loads_pp   =  zero
  !pmul_pp    =  zero
  !temp_pp    =  zero

!-------------------------------------------------------------------------
! 7. Initialise the solution vector to 0.0
!-------------------------------------------------------------------------
  
  ! U_n = U
  xnew_pp = 0._iwp
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

 num_load_steps=1

  DO iload = 1,num_load_steps    
    converged = .FALSE.
   
   ! Commented to check validity
   fext_pp(1:) = FLOAT(iload)*fextpiece_pp(1:)
    
!------------------------------------------------------------------------------
!----------------------- Start Newton-Raphson iterations ----------------------
!------------------------------------------------------------------------------
    inewton = 0
    iterations: DO
      inewton = inewton + 1

      storefint_pp = 0._iwp

      xnewel_pp = zero
      CALL GATHER(xnew_pp(1:),xnewel_pp)


  !-------------------------------------------------------------------------
  ! x. Build Matricies (K , M , f_int)
  !-------------------------------------------------------------------------
  ! Clean [K] and [M] 
  storekm_pp  =  zero
  storemm_pp  =  zero

      DO iel = 1,nels_pp
        kmat_elem = 0._iwp
        kgeo_elem = 0._iwp
        DO i = 1,nod
          num(i) = g_num_pp(i,iel) !- nn_start + 1
        END DO
        !coord = TRANSPOSE(g_coord_pp(:,:,iel))
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

      END DO ! nels_pp
 

      ! F_int
      fint_pp(:) = .0_iwp
      CALL SCATTER(fint_pp(1:),storefint_pp)

!-------------------------------------------------------------------------
! X. Get residual
!-------------------------------------------------------------------------

      !r_pp(1:) = fext_pp(1:) - fint_pp(1:)       !Residual
      !r_pp(0) = .0_iwp
	  
      ! Compute maxdiff of residual 
      !maxdiff =  MAXABSVAL_P(r_pp(1:))

      ! Normalise residual vector and stiffness matrix for pcg
      !IF (maxdiff == 0.0) THEN
      !  EXIT
      !END IF

!-------------------------------------------------------------------------
! X. Set Newmark Parameters
!-------------------------------------------------------------------------

     alpha = 0.25
     beta = 0.0
     delta = 0.5

     a0 = 1.0/(alpha*(dtim**2.0))
     a1 = zero
     a2 = 1.0/(alpha*dtim)
     a3 = (1.0/( 2.0*alpha)) -1.0
     a4 = zero
     a5 = zero
     a6 = 1.0/(alpha*(dtim**2.0))
     a7 = -1.0/(alpha*dtim)
     a8 = -( (1.0/(2.0*alpha) ) -1.0)
     a9 = dtim*(1-delta)
     a10 = delta*dtim

     ! M_eff
     meff_pp = zero
     meff_pp(1:) = a0*(x0_pp(1:)-xnew_pp(1:)) + a2*d1x0_pp(1:) +a3*d2x0_pp(1:)

     temp_pp = storemm_pp
     pmul_pp = zero

     ! M*M_eff
     CALL GATHER(meff_pp(1:),pmul_pp) ; utemp_pp=zero

     DO iel=1,nels_pp
       CALL DGEMV('N',ntot,ntot,one,temp_pp(:,:,iel),ntot,                 &
                pmul_pp(:,iel),1,zero,utemp_pp(:,iel),1)
     END DO

     vu_pp = zero
     CALL scatter(vu_pp(1:),utemp_pp)

     r_pp(1:) = fext_pp(1:) - fint_pp(1:) + vu_pp(1:)

     storekm_pp = storekm_pp + a0*storemm_pp

!-------------------------------------------------------------------------
! X. Diagonal Preconditioner
!-------------------------------------------------------------------------
      diag_precon_tmp = .0_iwp
      DO iel = 1,nels_pp
        DO k = 1,ntot 
          diag_precon_tmp(k,iel)=diag_precon_tmp(k,iel) + storekm_pp(k,k,iel)
        END DO
      END DO
  
!     Input: diag_precon_tmp r(ntot,nels_pp): Diagonal preconditioner at 
!                                             element level
      diag_precon_pp(:) = .0_iwp
      CALL SCATTER(diag_precon_pp(1:),diag_precon_tmp)

!     Output: diag_precon_pp r(1:neq_pp) Diagonal preconditioner assembled

      diag_precon_pp(1:) = 1._iwp/diag_precon_pp(1:)
      diag_precon_pp(0)  = .0_iwp

!---------------------------------------------------------------------------
!------------------------------- Solve using PCG ---------------------------
!---------------------------------------------------------------------------
      deltax_pp = .0_iwp
      res_pp    = r_pp

      ! Include dynamic Effects

      CALL PCG_VER1(inewton,limit,tol,storekm_pp,r_pp(1:), &
                    diag_precon_pp(1:),rn0,deltax_pp(1:),iters)

      !IF(numpe .EQ. 1)WRITE(*,'(2(a,I3))'),"N-R: ",inewton," PCG iters: ",iters

      xnew_pp(1:) = xnew_pp(1:) + deltax_pp(1:)
      xnew_pp(0) = .0_iwp

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

! Check convergence for Newton-Raphson iterations 
      energy = ABS(DOT_PRODUCT_P(res_pp(1:),deltax_pp(1:)))
      IF (inewton==1) THEN
       energy1 = energy
      END IF

      !if(numpe .EQ. 1)PRINT*,iload,inewton,energy,energy/energy1

      IF (inewton>1) THEN
        IF ((energy/energy1)<=tol2) THEN
          converged = .TRUE.
        END IF 
      END IF 

      IF(converged .OR. inewton==100) THEN
        EXIT
      END IF

    END DO iterations
   
   IF(numpe .EQ. 1)WRITE(*,'(a,I3)')," Newton-Raphson Iters: ",inewton
 
  END DO !iload

   x1_pp=zero; d2x1_ppstar= zero; d1x1_pp= zero; d2x1_pp=zero;
 
   x1_pp(1:)       = xnew_pp(1:)

   d2x1_ppstar(1:) = a6*(x1_pp(1:)-x0_pp(1:)) + a7*d1x0_pp(1:) + a8 * d2x0_pp(1:)

   d1x1_pp(1:)     = d1x0_pp(1:) + a9*d2x0_pp(1:) + a10*d2x1_ppstar(1:)
   d2x1_pp(1:)     = d2x1_ppstar(1:)
  
   !x1_pp	=  xnew_pp
   !d1x1_pp	=  (x1_pp-x0_pp)/(theta*dtim)-d1x0_pp*(1._iwp-theta)/theta
   !d2x1_pp	=  (d1x1_pp-d1x0_pp)/(theta*dtim)-d2x0_pp*(1._iwp-theta)/theta
   !x0_pp	=  x1_pp;
   !d1x0_pp	=  d1x1_pp
   !d2x0_pp	=  d2x1_pp;
     
!------------------------------------------------------------------------------
! 11. Gather Data from ntot,nels_pp to ndim,nodes_pp
!------------------------------------------------------------------------------

  IF(.NOT.ALLOCATED(eld_pp))THEN
   ALLOCATE(eld_pp(ntot,nels_pp))
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

  END SUBROUTINE

  !------------------------------------------------------------------------
  !------------------------------------------------------------------------
  !------------------------------------------------------------------------


