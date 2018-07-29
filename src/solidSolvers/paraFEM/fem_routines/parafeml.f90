  !/****h* solvers/parafeml
  !*  NAME
  !*    parafeml - Routines to solve 3D linear elastic 
  !*               deformation using the small strain assumption.
  !*
  !*  SYNOPSIS
  !*    The group of subroutines makes up the routines required to
  !*    solve the forced vibration of a 3D linear elastic solid
  !*    assuming the material undergoes finite strain. 
  !*
  !*    The governing equations solved:
  !*    
  !*    [M]{a}+[C]{u}+[K]{d} = {f_t}
  !*
  !*    The equations are stepped through using the linear interpolation
  !*    in time using theta.
  !*
  !*  FUNCTION
  !*    These routines are based on the decomposition of program 12.9 
  !*    found in "Programming the Finite Element Method".
  !*    P12.9 is the parallel analysis of the forced vibration of a
  !*    linear elastic solid.
  !*    
  !*    Subroutine           Purpose    
  !*
  !*    initl                Generates initial matricies and arrays
  !*    runl                 Solves the governing equations
  !*
  !*  AUTHOR
  !*    S.Hewitt
  !*  COPYRIGHT
  !*    (c) University of Manchester 1996-2017
  !******
  !*/

  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
  !--------------------------------------------------------------------

  SUBROUTINE initl(g_coord,rest,nn,nr,g_num_pp,g_g_pp,g_coord_pp)
  !/****f* parafeml/initl
  !*  NAME
  !*    SUBROUTINE: initl
  !*
  !*  SYNOPSIS
  !*    Usage:  initl_(mPoints_,rest_,&gPoints_,&numRestrNodes_,
  !*                   g_num_pp_OF_,g_g_pp_OF_,g_coord_pp_OF_);
  !*            
  !*  FUNCTION
  !*    Initialises ParaFEM, this inculdes initialising the MPI
  !*    Passing mesh information from OpenFOAM into ParaFEM and  
  !*    creating the steering matrix (g_g_pp).
  !*
  !*  INPUTS
  !*    g_coord   (ndim,nn)             - Coordinates of the mesh
  !*    rest      (nr,nodof+1)          - Restrained Nodes, e.g.(# x y z)
  !*
  !*    nn                              - Number of Nodes
  !*    nr                              - Number of restrained Nodes 
  !*
  !*    g_num_pp  (nod,nels_pp)         - Distributed element steering array
  !*
  !*  OUTPUT
  !*    g_g_pp      (ntot,nels_pp)      - Distributed equation steering array
  !*    g_coord_pp  (nod,ndim,nels_pp)  - Distributed nodal cooridnates 
  !*
  !*  AUTHOR
  !*    S. Hewitt
  !******
  !*  COMMENTS
  !*  initl and initnl are the same and so one will be removed 
  !*  prior to the next release  
  !*/

  USE mpi_wrapper;    USE precision;  USE global_variables; 
  USE mp_interface;   USE input;      USE output; 
  USE loading;        USE timing;     USE maths; 
  USE gather_scatter; USE steering;   USE new_library;
  USE large_strain;

  IMPLICIT NONE

!----------------------------------------------------------------------
! 1. Declare variables
!----------------------------------------------------------------------

  INTEGER,PARAMETER        :: nodof=3,ndim=3,nst=6,nod=8
  REAL(iwp),PARAMETER      :: zero=0.0_iwp

  INTEGER,INTENT(IN)       :: nn,nr

  INTEGER,INTENT(INOUT)    :: g_g_pp(ndim*nod,nels_pp),rest(nr,nodof+1)
  INTEGER,INTENT(INOUT)    :: g_num_pp(nod,nels_pp)

  INTEGER                   :: iel,nip,npes_pp,partitioner,ndof
  INTEGER                   :: nlen,i

  REAL(iwp),INTENT(INOUT)  :: g_coord(ndim,nn)
  REAL(iwp),INTENT(INOUT)  :: g_coord_pp(nod,ndim,nels_pp)

  CHARACTER(LEN=50)        :: argv
  CHARACTER(LEN=15)        :: element
 
  LOGICAL                  :: initialised

!----------------------------------------------------------------------
! 2. Input and Initialisation
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
! 3. Populate  g_num_pp
!----------------------------------------------------------------------
  
  ! Convert from Foam-Extend to Smith Gritths format
  DO iel=1,nels_pp
    CALL of2sg(element,g_num_pp(:,iel),nod)
  ENDDO
  
!----------------------------------------------------------------------
! 4. Calculate g_g_pp
!----------------------------------------------------------------------

  ! Rearrange the rest Array
  CALL rearrange(rest)

  ! Clean arrays
  g_g_pp  =  zero
  neq     =  zero
  
  ! Find the global Steering Matrix
  elements_0: DO iel=1,nels_pp
    !CALL find_g(g_num_pp(:,iel),g_g_pp(:,iel),rest) ! Stable but slow
    CALL find_g3(g_num_pp(:,iel),g_g_pp(:,iel),rest) ! Unstable but Fast
  END DO elements_0
  
  ! Build GGL Array
  neq  =  MAXVAL(g_g_pp)
  neq  =  max_p(neq)
  CALL calc_neq_pp

  ! Most failures occur in this routine
  ! CALL calc_npes_pp(npes,npes_pp)
  ! Set npes_pp
  SELECT CASE (npes)
      CASE (1:15)
        npes_pp = npes
      CASE (16:32)
        npes_pp = npes
      CASE (33:256)
        npes_pp = npes/2
      CASE (257:1024)
        npes_pp = npes/4
      CASE DEFAULT
        npes_pp = npes/8
  END SELECT
 
  CALL make_ggl(npes_pp,npes,g_g_pp)

!----------------------------------------------------------------------
! 5. Populate g_coord_pp
!----------------------------------------------------------------------

  CALL POPULATE_G_COORD_PP2(g_coord,g_coord_pp,g_num_pp,nn,nod,ndim) 
 
  ! output g_g_pp,g_num_pp,g_coord_pp

  END SUBROUTINE

  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
  !--------------------------------------------------------------------

  SUBROUTINE runl(node,val,num_var,mat_prop,nr,loaded_nodes,time_step, &
                  g_g_pp,g_num_pp,g_coord_pp,gravlo_pp,Dfield,Ufield,Afield,flag)
  !/****f* parafeml/runl
  !*  NAME
  !*    SUBROUTINE: runl
  !*
  !*  SYNOPSIS
  !*    Usage:     runl_(forceNodes_,fext_OF_,numSchemes_,solidProps_,
  !*                     &numRestrNodes_,&numFixedForceNodes_, &dtim,
  !*                     g_g_pp_OF_,g_num_pp_OF_,g_coord_pp_OF_,gravlo_,
  !*                     ptDtemp_,ptUtemp_,ptAtemp_);
  !* 
  !*  FUNCTION
  !*    Reads in the current timesteps displacement, velocity, 
  !*    acceleration and external force field. Loads the structure    
  !*    and solves the governing equations of the problem
  !*
  !*        {F} = [M]{a} + [C]{u} + [K]{d} 
  !*
  !*    The new displacement, velocity and acceleration fields are
  !*    output. Note this subroutine is used for problems with infinte strain.
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
    
  USE mpi_wrapper;    USE precision; USE global_variables; 
  USE mp_interface;   USE input;     USE output; 
  USE loading;        USE timing;    USE maths; 
  USE gather_scatter; USE steering;  USE new_library; 
  
  IMPLICIT NONE
  
!------------------------------------------------------------------------------ 
! 1. Declare variables used in the main program
!------------------------------------------------------------------------------

  INTEGER,PARAMETER        :: nodof=3,ndim=3,nst=6,nod=8
  REAL(iwp),PARAMETER      :: zero=0.0_iwp,one=1.0_iwp
  
  INTEGER,INTENT(INOUT)    :: loaded_nodes,node(loaded_nodes)
  INTEGER,INTENT(INOUT)    :: g_num_pp(nod,nels_pp)
  INTEGER,INTENT(INOUT)    :: g_g_pp(ntot,nels_pp),nr
  
  INTEGER                  :: iel,i,j,k,l,m,n,iters,printres
  INTEGER                  :: limit,nels,node_end,node_start
  INTEGER                  :: nlen,myCount,disps(npes)
  INTEGER                  :: nodesCount(npes),RSS,VM,RSSa,VMa  

  REAL(iwp),INTENT(IN)      :: mat_prop(3),flag
  REAL(iwp),INTENT(IN)      :: num_var(5),time_step,gravlo_pp(neq_pp)

  REAL(iwp),INTENT(INOUT)  :: g_coord_pp(nod,ndim,nels_pp)
  REAL(iwp),INTENT(INOUT)  :: val(ndim,loaded_nodes)
  REAL(iwp),INTENT(INOUT)  :: Dfield(ntot,nels_pp),Ufield(ntot,nels_pp)
  REAL(iwp),INTENT(INOUT)  :: Afield(ntot,nels_pp)
  
  REAL(iwp)                :: tol,up,alpha,beta,alpha1,beta1,theta,dtim
  REAL(iwp)                :: c1,c2,c3,c4
  REAL(iwp)                :: X,Y,Z

  REAL(iwp),SAVE           :: time_step_orig
  
  LOGICAL                  :: converged
  CHARACTER(LEN=50)        :: argv
  CHARACTER(LEN=15)        :: element;
  CHARACTER(LEN=80)        :: FMT
  CHARACTER(LEN=1024)      :: filename  
  
!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------
  REAL(iwp),SAVE,ALLOCATABLE:: loads_pp(:),fext_pp(:),x1_pp(:),d1x1_pp(:)
  REAL(iwp),SAVE,ALLOCATABLE:: d2x1_pp(:),x0_pp(:),d1x0_pp(:),d2x0_pp(:)
  REAL(iwp),SAVE,ALLOCATABLE:: vu_pp(:),u_pp(:),p_pp(:),d_pp(:),x_pp(:)
  REAL(iwp),SAVE,ALLOCATABLE:: xnew_pp(:),pmul_pp(:,:),utemp_pp(:,:)
  REAL(iwp),SAVE,ALLOCATABLE:: temp(:),d1temp(:),d2temp(:),temp_pp(:,:,:)
  REAL(iwp),SAVE,ALLOCATABLE:: disp_pp(:),vel_pp(:),acel_pp(:),eld_pp(:,:)
  REAL(iwp),SAVE,ALLOCATABLE:: gDisp(:),gVel(:),gAcel(:)
  REAL(iwp),SAVE,ALLOCATABLE:: fext_o_pp(:),timest(:),fext_tmp_pp(:)
  REAL(iwp),SAVE,ALLOCATABLE:: diag_precon_pp(:),store_km_pp(:,:,:)
  REAL(iwp),SAVE,ALLOCATABLE:: store_mm_pp(:,:,:)
  
  
!------------------------------------------------------------------------------
! 3. Start Program
!------------------------------------------------------------------------------
  IF(numpe .EQ. 1)PRINT*,"ParaFEM Small Strain Solver:"

  IF(.NOT. ALLOCATED(diag_precon_pp))THEN
    
    time_step_orig = time_step

    ALLOCATE(diag_precon_pp(neq_pp),store_km_pp(ntot,ntot,nels_pp),&
             store_mm_pp(ntot,ntot,nels_pp))

    store_km_pp = zero; store_mm_pp=zero;diag_precon_pp=zero;
    CALL finddiagprecon(store_km_pp,store_mm_pp,g_coord_pp,num_var,mat_prop, &
                   diag_precon_pp,time_step)
  ENDIF

  IF(time_step_orig .NE. time_step)THEN
    ! NOTE : Accuracy of this in FSI context not been checked
    IF(numpe .EQ. 1)PRINT*,"Time Step Updated"
 
    time_step_orig = time_step
    store_km_pp = zero; store_mm_pp=zero;diag_precon_pp=zero;
    CALL finddiagprecon(store_km_pp,store_mm_pp,g_coord_pp,num_var,mat_prop, &
                   diag_precon_pp,time_step)
  ENDIF
  
  ! Barrier (may not be needed but safe)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
  
  ! Set Base paramenter
  argv     =  "Case"        ! Name files write to
  nlen     =  4             ! Length of Name
  element  =  "hexahedron"  ! Element Name
  
  ! Set Numerical and Material Values 
  alpha1  =  num_var(1)    ! Rayleigh damping parameters
  beta1   =  num_var(2)    ! Rayleigh damping parameters
  theta   =  num_var(3)    ! Linear time interpolator (0.5 to 1)
  tol     =  num_var(4)    ! Tolerance of PCG loop
  limit   =  num_var(5)    ! Max number of Interation in PCG

  dtim   =  time_step

  c1      =  (1._iwp-theta)*dtim
  c2      =  beta1-c1
  c3      =  alpha1+1._iwp/(theta*dtim);
  c4      =  beta1+theta*dtim

  printres=0

  ! Allocate memory required for the time loop
  IF(.NOT.ALLOCATED(timest))THEN
   CALL system_mem_usage(RSSa,VMa)
   ALLOCATE(x0_pp(neq_pp),d1x0_pp(neq_pp),x1_pp(neq_pp),vu_pp(neq_pp))
   ALLOCATE(u_pp(neq_pp),d2x0_pp(neq_pp),loads_pp(neq_pp))
   ALLOCATE(d1x1_pp(neq_pp),d2x1_pp(neq_pp),d_pp(neq_pp),p_pp(neq_pp))
   ALLOCATE(x_pp(neq_pp),timest(20),xnew_pp(neq_pp),fext_pp(neq_pp))
   ALLOCATE(pmul_pp(ntot,nels_pp))
   ALLOCATE(utemp_pp(ntot,nels_pp),temp_pp(ntot,ntot,nels_pp))
   ALLOCATE(fext_o_pp(neq_pp),fext_tmp_pp(neq_pp))
   fext_o_pp  =  zero; fext_tmp_pp = zero
  ENDIF
    
  timest=zero 
  timest(1)=elap_time()
  
  ! Clean Arrays
  x0_pp    =  zero;  d1x0_pp =  zero;  x1_pp    =  zero; 
  vu_pp    =  zero;  u_pp    =  zero;  d2x0_pp  =  zero; 
  loads_pp =  zero;  d1x1_pp =  zero;  d2x1_pp  =  zero;
  d_pp     =  zero;  p_pp    =  zero;  x_pp     =  zero; 
  xnew_pp  =  zero;

!------------------------------------------------------------------------------
! 4. Set Loads
!------------------------------------------------------------------------------

  timest(2)  =  elap_time()

  fext_pp    =  zero

  ! Load fext_pp based on global load vector
  CALL load(g_g_pp,g_num_pp,node,val,fext_pp)
 
  CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
 
  fext_pp=fext_pp+gravlo_pp  ! Adding gravity if set

  timest(3)=elap_time()
  
!------------------------------------------------------------------------------
! 5. Set Initial Conditions
!------------------------------------------------------------------------------

! - scatter_noadd has no barriers in the subroutine
  CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
  CALL scatter_noadd(Dfield,x0_pp)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
  CALL scatter_noadd(Ufield,d1x0_pp)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
  CALL scatter_noadd(Afield,d2x0_pp)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

  u_pp      =  zero
  vu_pp     =  zero
  loads_pp  =  zero
  pmul_pp   =  zero
  temp_pp   =  zero

  timest(4)=elap_time()
  
!------------------------------------------------------------------------------
! 6. Displacement
!------------------------------------------------------------------------------
   temp_pp  =  store_km_pp*c2+store_mm_pp*c3
   
   CALL gather(x0_pp,pmul_pp)
   
   DO iel=1,nels_pp
     CALL DGEMV('N',ntot,ntot,one,temp_pp(:,:,iel),ntot,                 &
                pmul_pp(:,iel),1,zero,utemp_pp(:,iel),1)
   END DO
   
   CALL scatter(u_pp,utemp_pp)
   
!------------------------------------------------------------------------------
! 7. Velocity
!------------------------------------------------------------------------------
   temp_pp=store_mm_pp/theta
   
   CALL gather(d1x0_pp,pmul_pp);utemp_pp=zero
   
   DO iel=1,nels_pp
     CALL DGEMV('N',ntot,ntot,one,temp_pp(:,:,iel),ntot,                 &
                pmul_pp(:,iel),1,zero,utemp_pp(:,iel),1)
   END DO
   
   CALL scatter(vu_pp,utemp_pp)

   timest(5) =  elap_time()
   
!------------------------------------------------------------------------------
! 8. Bulid Right Hand Side
!------------------------------------------------------------------------------

  loads_pp     =  fext_pp*theta*dtim+(1-theta)*dtim*fext_o_pp

  IF(flag .EQ. 2)THEN
      fext_o_pp    =  fext_tmp_pp
  ENDIF

  fext_tmp_pp    =  fext_pp
  loads_pp     =  u_pp+vu_pp+loads_pp
  temp_pp      =  store_mm_pp*c3+store_km_pp*c4

  timest(6) =  elap_time()
  
!------------------------------------------------------------------------------
! 9. PCG
!------------------------------------------------------------------------------

  IF(ABS(SUM_P(loads_pp)) .GE. 1E-10)THEN
   
    d_pp  =  diag_precon_pp*loads_pp
    p_pp  =  d_pp
    x_pp  =  zero
    iters =  0
   
    iterations: DO
      iters =  iters+1
      u_pp  =  zero
      vu_pp =  zero
     
      CALL gather(p_pp,pmul_pp)
     
      elements_4: DO iel=1,nels_pp
        CALL DGEMV('N',ntot,ntot,one,temp_pp(:,:,iel),ntot,   &
                   pmul_pp(:,iel),1,zero,utemp_pp(:,iel),1)
      END DO elements_4; 
      
      CALL scatter(u_pp,utemp_pp)
      
      up        =  DOT_PRODUCT_P(loads_pp,d_pp);
      alpha	=  up/DOT_PRODUCT_P(p_pp,u_pp)
      xnew_pp	=  x_pp+p_pp*alpha; 
      loads_pp	=  loads_pp-u_pp*alpha
      d_pp	=  diag_precon_pp*loads_pp; 
      beta	=  DOT_PRODUCT_P(loads_pp,d_pp)/up  
      p_pp	=  d_pp+p_pp*beta; 
      u_pp	=  xnew_pp
      
      CALL checon_par(xnew_pp,tol,converged,x_pp)
      
      IF(converged.OR.iters==limit)EXIT
    END DO iterations
   
!------------------------------------------------------------------------------
! 10. PCG END
!------------------------------------------------------------------------------
   !CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
   
   IF(numpe .EQ. 1)PRINT*,"Number of PCG Iterations: ",iters
   ENDIF
   x1_pp      =  xnew_pp
   utemp_pp   =  zero
   d1x1_pp    =  (x1_pp-x0_pp)/(theta*dtim)-d1x0_pp*(1._iwp-theta)/theta
   d2x1_pp    =  (d1x1_pp-d1x0_pp)/(theta*dtim)-d2x0_pp*(1._iwp-theta)/theta
   x0_pp      =  x1_pp;
   d1x0_pp    =  d1x1_pp
   d2x0_pp    =  d2x1_pp;
   
  timest(7) = elap_time()
  
!------------------------------------------------------------------------------
! 11. Gather Data from ntot,nels_pp to ndim
!------------------------------------------------------------------------------

  IF(.NOT.ALLOCATED(eld_pp))THEN
   ALLOCATE(eld_pp(ntot,nels_pp))
   printres=1
  END IF

  eld_pp   =  zero
 
  ! Displacement
  CALL gather(x1_pp(1:),eld_pp)
  Dfield=eld_pp
  
  ! Velocity
  CALL gather(d1x1_pp(1:),eld_pp)
  Ufield=eld_pp

  ! Acceleration
  CALL gather(d2x1_pp(1:),eld_pp)
  Afield=eld_pp

  timest(8)=elap_time()
  
!------------------------------------------------------------------------------
! 12. Print Runtime Information to File
!------------------------------------------------------------------------------

  IF(numpe==1 .AND. printres==1)THEN
    !CALL system_mem_usage(RSS,VM)
    printres=0
    CALL WRITE_SMALLSTRAIN(argv,timest,iters)
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
  END SUBROUTINE

  !------------------------------------------------------------------------
  !------------------------------------------------------------------------
  !------------------------------------------------------------------------

