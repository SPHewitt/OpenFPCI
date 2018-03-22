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
  !*    finddiagl            Finds the diagonal preconditioner
  !*    runl                 Solves the governing equations
  !*
  !*  AUTHOR
  !* 	  S.Hewitt
  !*  COPYRIGHT
  !*    (c) University of Manchester 1996-2017
  !******
  !*/

  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
  !--------------------------------------------------------------------

  SUBROUTINE initl(g_coord,g_num_pp,rest,nn,nels,nr,sProp,g_g_pp,store_km_pp,store_mm_pp)

  !/****f* parafeml/initl
  !*  NAME
  !*    SUBROUTINE: initl
  !*
  !*  SYNOPSIS
  !*    Usage: initl_(mPoints_,g_num_pp_OF_,rest_,&gPoints,
  !*                   &gCells,&numRestrNodes_,solidProps_,
  !*                   g_g_pp_OF_,store_km_pp_OF_,store_mm_pp_OF_);
  !*            
  !*  FUNCTION
  !*    Initialises ParaFEM, this inculdes initialising the MPI
  !*    calculating the stiffness matrix [k], Mass Matrix [M] 
  !*	and steering matrix (g_g_pp).
  !*
  !*  INPUTS
  !*    g_coord   (ndim,nn)     - Coordinates of the mesh		
  !*    g_num_pp  (nod,nels_pp) - Steering matrix				
  !*    rest      (nr,nodof+1)  - Restrained Nodes, e.g.(# x y z)
  !*
  !*    nn                      - Number of Nodes
  !*    nels                    - Number of elements
  !*    nr                      - Number of restrained Nodes
  !*
  !*    solidProp (e v rho)     - Solid Properties 	 
  !*	  			
  !*  OUTPUT
  !*    g_g_pp      (ntot,nels_pp)      - Global Steering Matrix
  !*   	store_km_pp (ntot,ntot,nels_pp) - Stiffness Matrix [k]
  !*   	store_mm_pp (ntot,ntot,nels_pp) - Mass Matrix [M]

  !*			
  !*  AUTHOR
  !*    S. Hewitt
  !******
  !*  COMMENTS
  !*    nels is currently deprecated, however it is possibly used to
  !*    decompose the mesh using calc_nels_pp.
  !*/

  USE mpi_wrapper;    USE precision;  USE global_variables; 
  USE mp_interface;   USE input;      USE output; 
  USE loading;        USE timing;     USE maths; 
  USE gather_scatter; USE steering;   USE new_library;

  IMPLICIT NONE

!----------------------------------------------------------------------
! 1. Declare variables
!----------------------------------------------------------------------

  INTEGER,PARAMETER		    :: nodof=3,ndim=3,nst=6,nod=8
  REAL(iwp),PARAMETER		  :: zero=0.0_iwp

  INTEGER,INTENT(IN) 		  :: nels,nn

  INTEGER,INTENT(INOUT)		:: g_g_pp(ndim*nod,nels_pp),rest(nr,nodof+1)
  INTEGER,INTENT(INOUT)		:: g_num_pp(nod,nels_pp)

  INTEGER                 :: iel,i,j,k,l,m,n,nr,nip,ndof,npes_pp,nlen
  INTEGER			            :: partitioner,printres,RSS,VM,RSSa,VMa

  REAL(iwp),INTENT(INOUT)	:: g_coord(ndim,nn),sProp(3)
  REAL(iwp),INTENT(INOUT)	:: store_km_pp(ndim*nod,ndim*nod,nels_pp)
  REAL(iwp),INTENT(INOUT)	:: store_mm_pp(ndim*nod,ndim*nod,nels_pp)

  REAL(iwp)			          :: det,period,volume,tol,e,v,rho	

  LOGICAL			            :: converged=.false.
  LOGICAL			            :: consistent=.TRUE.
  LOGICAL			            :: initialised	

  CHARACTER(LEN=50)		    :: argv
  CHARACTER(LEN=15)		    :: element
  CHARACTER(LEN=6) 		    :: ch
  CHARACTER(LEN=1024) 		:: filename  

!----------------------------------------------------------------------
! 2. Declare dynamic arrays
!----------------------------------------------------------------------

  REAL(iwp),ALLOCATABLE		:: points(:,:),dee(:,:),weights(:),timest(:)
  REAL(iwp),ALLOCATABLE		:: jac(:,:),der(:,:),deriv(:,:),bee(:,:)
  REAL(iwp),ALLOCATABLE		:: fun(:),emm(:,:),ecm(:,:),g_coord_pp(:,:,:)

  INTEGER,ALLOCATABLE		  :: node(:),localcount(:),readcount(:)
 
  ! Initilal System Memory Usage
  CALL system_mem_usage(RSSa,VMa)
 
!----------------------------------------------------------------------
! 3. Input and Initialisation
!---------------------------------------------------------------------- 

  ! Variables required for writing and partioning
  argv="Case"; nlen=4; nip=8; element="hexahedron"; partitioner=1


  ! Youndgs Modulus, Poissons Ratio and Density
  e=sProp(1);v=sProp(2);rho=sProp(3);

  ALLOCATE(timest(20))
  timest	  =  zero
  timest(1)	=  elap_time()

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

  ! Degrees of Freedon per Element
  ndof  =  nod*nodof
  ntot  =  ndof
 
!----------------------------------------------------------------------
! 4. Populate g_coord_pp and g_num_pp
!----------------------------------------------------------------------
  timest(2)=elap_time();
  
  ALLOCATE(g_coord_pp(nod,ndim,nels_pp))
  g_coord_pp  =  zero
  
  ! Poulate the Steering matrix
  !CALL POPULATE_G_NUM_PP(g_num,g_num_pp,npes,nod,nels)

  ! Convert from Foam-Extend to Smith Gritths format
  DO iel=1,nels_pp
    CALL of2sg(element,g_num_pp(:,iel),nod)
  ENDDO

  ! Populate Coorinate Matrix
  !CALL POPULATE_G_COORD_PP(g_coord,g_coord_pp,g_num_pp,npes,nn,nod,ndim)
  CALL POPULATE_G_COORD_PP2(g_coord,g_coord_pp,g_num_pp,nn,nod,ndim) 
 
  timest(3)=elap_time();
  
!----------------------------------------------------------------------
! 5. Find the Steering arrays and equations per core
!----------------------------------------------------------------------

  ALLOCATE(points(nip,ndim),dee(nst,nst),jac(ndim,ndim))
  ALLOCATE(der(ndim,nod),deriv(ndim,nod),bee(nst,ntot))
  ALLOCATE(weights(nip),ecm(ntot,ntot),emm(ntot,ntot),fun(nod))

  ! Rearrange the rest Array
  CALL rearrange(rest)

  ! Clean arrays
  g_g_pp  =  zero
  neq	    =  zero
  
  ! Find the global Steering Matrix
  elements_0: DO iel=1,nels_pp
  ! CALL find_g(g_num_pp(:,iel),g_g_pp(:,iel),rest) ! Stable but slow
    CALL find_g3(g_num_pp(:,iel),g_g_pp(:,iel),rest)! Unstable but Fast
  END DO elements_0
  
  ! Build GGL Array
  neq  =  MAXVAL(g_g_pp)
  neq  =  max_p(neq)
  CALL calc_neq_pp
  CALL calc_npes_pp(npes,npes_pp)
  CALL make_ggl(npes_pp,npes,g_g_pp)
  
  timest(4)=elap_time()	

!----------------------------------------------------------------------
! 6. Element Stiffness and Mass Integration
!----------------------------------------------------------------------
  
  ! [D] : Stress-Strain Matrix
  dee  =  zero 
  CALL deemat(dee,e,v)
  CALL sample(element,points,weights)
  
  ! Clean [K] and [M] 
  store_km_pp  =  zero
  store_mm_pp  =  zero
  
  elements_2: DO iel=1,nels_pp
    volume  	=  zero
    emm  	=  zero
    ecm  	=  zero
    
    gauss_points_1: DO i=1,nip  
      ! [N] : Shape Functions   
      CALL shape_der(der,points,i)
      
      ! [J] : Jacobian, global to local derivative
      jac=MATMUL(der,g_coord_pp(:,:,iel))
      
      ! det|J|
      det=determinant(jac)
      CALL invert(jac)
      deriv=matmul(jac,der)
      
      ! [B] : Strain - Displacement Matrix 
      CALL beemat(bee,deriv)
      
      ! Integration for [K]
      store_km_pp(:,:,iel)=store_km_pp(:,:,iel) +                     &
          MATMUL(MATMUL(TRANSPOSE(bee),dee),bee) * det*weights(i)
     volume=volume+det*weights(i)
     CALL shape_fun(fun,points,i)
     
     ! Integration for [M]
     IF(consistent)THEN
       CALL ecmat(ecm,fun,ntot,nodof)
       ecm  =  ecm*det*weights(i)*rho
       emm  =  emm+ecm
     END IF
    END DO gauss_points_1   
    IF(.NOT.consistent)THEN
      DO i=1,ntot; emm(i,i)=volume*rho/13._iwp; END DO
      DO i=1,19,6; emm(i,i)=emm(4,4)*.125_iwp; END DO
      DO i=2,20,6; emm(i,i)=emm(4,4)*.125_iwp; END DO
      DO i=3,21,6; emm(i,i)=emm(4,4)*.125_iwp; END DO
      DO i=37,55,6; emm(i,i)=emm(4,4)*.125_iwp; END DO
      DO i=38,56,6; emm(i,i)=emm(4,4)*.125_iwp; END DO
      DO i=39,57,6; emm(i,i)=emm(4,4)*.125_iwp; END DO
    END IF
    store_mm_pp(:,:,iel)  =  emm
  END DO elements_2  
  
  timest(5)=elap_time()	
  
!----------------------------------------------------------------------
! 6. Print Information about Runtime to File
!---------------------------------------------------------------------- 

  ! Manual Switch to turn on/off .res file
  printres  =  1
  
  IF(numpe==1 .AND. printres==1)THEN
  CALL system_mem_usage(RSS,VM)
  OPEN(11,FILE=argv(1:nlen)//".res",STATUS='REPLACE',ACTION='WRITE')
    WRITE(11,'(A,I7,A)') "This job ran on ",npes," processes"
    WRITE(11,'(A,3(I12,A))') "There are ",nn," nodes", nr, &
                            " restrained and ",neq," equations"
    WRITE(11,'(A,F10.4)') "-------------------------------------------------"
    WRITE(11,'(A,F10.4)') "Time to Populate g_coord/num_pp is:",timest(3)-timest(2)
    WRITE(11,'(A,F10.4)') "Time to find_g & make_ggl:",timest(4)-timest(3)
    WRITE(11,'(A,F10.4)') "Time to Build K and M:",timest(5)-timest(4)
    WRITE(11,'(A,F10.4)') "Time in Routine(Total):",elap_time()-timest(1)
    WRITE(11,'(A,F10.4)') "-------------------------------------------------"
    WRITE(11,'(A,I10)') "Virtual Memory Change(Kb): ",VM-VMa
    WRITE(11,'(A,I10)') "RSS Memory Change(Kb): ",RSS-RSSa
    WRITE(11,'(A,I10)') "Total Virtual Memory(Kb): ",VM
    WRITE(11,'(A,I10)') "Total RSS Memory (Kb): ",RSS
  CLOSE(11)
  END IF
  
  DEALLOCATE(points,dee,jac,der,deriv,bee)
  DEALLOCATE(weights,ecm,emm,fun,g_coord_pp,timest)
 
  END SUBROUTINE

  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
  !--------------------------------------------------------------------

  SUBROUTINE finddiagl(store_km_pp,store_mm_pp,numVar,diag_precon_pp)
  !/****f* parafeml/finddiagl
  !*  NAME
  !*    SUBROUTINE: finddiagl
  !*
  !*  SYNOPSIS
  !*    Usage: finddiagl_(store_km_pp_OF_,store_mm_pp_OF_,      &
  !*                            numSchemes_,diag_precon_pp_OF_); 
  !*
  !*  FUNCTION
  !*    Calculates the diagonal preconditioner vector used 
  !*    in the PCG solution.
  !*
  !*  INPUTS
  !*    store_km_pp (ntot,ntot,nels_pp)  - Stiffness Matrix	
  !*    store_mm_pp (ntot,ntot,nels_pp)  - Mass Matrix
  !*
  !*    numVar      (a1 b1 theta dTim)   - Numerical Variables			
  !*	 			
  !*  OUTPUT
  !*   	diag_precon_pp (neq_pp)          - Diagonal Preconditioner
  !*			
  !*  AUTHOR
  !*    S. Hewitt
  !******
  !*/
  
  USE mpi_wrapper;	  USE precision;	USE global_variables; 
  USE mp_interface; 	USE input;	    USE output; 
  USE loading; 		    USE timing; 	  USE maths; 
  USE gather_scatter;	USE steering; 	USE new_library; 
  
  IMPLICIT NONE
  
!------------------------------------------------------------------------------ 
! 1. Declare variables used in the main program
!------------------------------------------------------------------------------

  INTEGER,PARAMETER		    :: nodof=3,ndim=3,nst=6, nod=8
  REAL(iwp),PARAMETER		  :: zero=0.0_iwp
  
  INTEGER			            :: i,j,k,iel,ndof

  REAL(iwp),INTENT(IN)		:: store_km_pp(ndim*nod,ndim*nod,nels_pp)
  REAL(iwp),INTENT(IN)		:: store_mm_pp(ndim*nod,ndim*nod,nels_pp)
  REAL(iwp),INTENT(IN)		:: numVar(4)
  
  REAL(iwp),INTENT(INOUT)	:: diag_precon_pp(neq_pp)
  
  REAL(iwp)			          :: alpha1,beta1,theta
  REAL(iwp)			          :: dtim,c3,c4
  
  REAL(iwp),ALLOCATABLE		:: timest(:),diag_precon_tmp(:,:)
  
  ! Set Parameters
  alpha1  =  numVar(1)
  beta1   =  numVar(2)
  theta   =  numVar(3)
  dtim	  =  numVar(4)
  c3	    =  alpha1+1._iwp/(theta*dtim)
  c4	    =  beta1+theta*dtim
  ndof    =  ntot
  
!---------------------------------------------------------------------- 
! 2. Calculate Diagonal preconditioner
!----------------------------------------------------------------------
  ALLOCATE(diag_precon_tmp(ntot,nels_pp))
  
  diag_precon_pp  =  zero
  diag_precon_tmp =  zero
  
  elements_2: DO iel=1,nels_pp
    DO k=1,ndof
      diag_precon_tmp(k,iel)=diag_precon_tmp(k,iel)    +             &
                  store_mm_pp(k,k,iel)*c3+store_km_pp(k,k,iel)*c4
    END DO
  END DO elements_2
  
  CALL scatter(diag_precon_pp,diag_precon_tmp); 
  
  DEALLOCATE(diag_precon_tmp)
  
  diag_precon_pp=1.0_iwp/diag_precon_pp 
  
  END SUBROUTINE

  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
  !--------------------------------------------------------------------

  SUBROUTINE runl(numVar,val,node,loaded_nodes,time,nodes_pp, 	&
                      g_g_pp,g_num_pp,store_km_pp,store_mm_pp,diag_precon_pp,gravlo,Dfield,Ufield,Afield)
  
  !/****f* parafeml/runl
  !*  NAME
  !*    SUBROUTINE: runl
  !*
  !*  SYNOPSIS
  !*    Usage: runl_(numSchemes_,fext_OF_,forceNodes_,          &
  !*            &numFixedForceNodes_,&time,&lPoints,g_g_pp_OF_,       &
  !*            g_num_pp_OF_,store_km_pp_OF_,store_mm_pp_OF_,         &
  !*            diag_precon_pp_OF_,gravlo_,ptDtemp_,ptUtemp_,ptAtemp_);
  !* 	
  !*  FUNCTION
  !*    Reads in the current timesteps displacement, velocity, 
  !*    acceleration and external force field. Loads the structure    
  !*    and solves the governing equations of the problem
  !*
  !*        {F} = [M]{a} + [C]{u} + [K]{d} 
  !*
  !*    The new displacement, velocity and acceleration fields are
  !*    output.
  !*
  !*  INPUTS
  !*    numVar      (a1 b1 theta dTim)   - Numerical Variables
  !*    val         (ndim,loaded_nodes)	 - Force vector of loaded nodes
  !*    node        (loaded_nodes)       - Node # of loaded_nodes
  !*
  !*    loaded_nodes                     - # of loaded nodes		  
  !*    time                             - Current time
  !*    nodes_pp                         - # of nodes per cores
  !*
  !*    g_g_pp      (ntot,nels_pp)       - Equation steering matrix
  !*    g_num_pp    (nod,nels_pp)        - Element steering matrix
  !*    store_km_pp (ntot,ntot,nels_pp)  - Stiffness Matrix [k]
  !*    store_mm_pp (ntot,ntot,nels_pp)  - Mass Matrix [M]
  !*
  !*    diag_precon_pp (neq_pp)          - Diagonal Preconditioner	
  !*    gravlo         (neq_pp)          - Vector of gravity loads	
  !*				
  !*  OUTPUT
  !*    Dfield  (ntot,nels_pp)           - Nodal displacements
  !*    Ufield  (ntot,nels_pp)           - Nodal velocities
  !*    Afield  (ntot,nels_pp)           - Nodal accelerations
  !* 				
  !*  AUTHOR
  !*    S. Hewitt
  !******		
  !*/  
    
  USE mpi_wrapper;	  USE precision;	USE global_variables; 
  USE mp_interface; 	USE input;	    USE output; 
  USE loading; 		    USE timing; 	  USE maths; 
  USE gather_scatter;	USE steering; 	USE new_library; 
  
  IMPLICIT NONE
  
!------------------------------------------------------------------------------ 
! 1. Declare variables used in the main program
!------------------------------------------------------------------------------

  INTEGER,PARAMETER		      :: nodof=3,ndim=3,nst=6,nod=8
  REAL(iwp),PARAMETER		    :: zero=0.0_iwp,one=1.0_iwp
  
  INTEGER,INTENT(INOUT)		  :: loaded_nodes,node(loaded_nodes)
  INTEGER,INTENT(INOUT)		  :: g_num_pp(nod,nels_pp),nodes_pp
  INTEGER,INTENT(INOUT)		  :: g_g_pp(ntot,nels_pp)
  
  INTEGER			              :: iel,i,j,k,l,m,n,iters,printres
  INTEGER			              :: limit,nels,node_end,node_start
  INTEGER			              :: nlen,myCount,disps(npes),flag
  INTEGER			              :: nodesCount(npes),RSS,VM,RSSa,VMa  

  REAL(iwp),INTENT(INOUT) 	:: diag_precon_pp(neq_pp)
  REAL(iwp),INTENT(INOUT) 	:: store_km_pp(ntot,ntot,nels_pp)
  REAL(iwp),INTENT(INOUT) 	:: store_mm_pp(ntot,ntot,nels_pp)
  REAL(iwp),INTENT(INOUT) 	:: numVar(4),time,gravlo(neq_pp)
  REAL(iwp),INTENT(INOUT) 	:: val(ndim,loaded_nodes)
  REAL(iwp),INTENT(INOUT)   :: Dfield(ntot,nels_pp),Ufield(ntot,nels_pp)
  REAL(iwp),INTENT(INOUT)   :: Afield(ntot,nels_pp)
  
  REAL(iwp)			            :: tol,up,alpha,beta,alpha1,beta1,theta,dtim
  REAL(iwp)			            :: c1,c2,c3,c4
  REAL(iwp)			            :: X,Y,Z
  
  LOGICAL			              :: converged
  CHARACTER(LEN=50)		      :: argv
  CHARACTER(LEN=15)		      :: element;
  CHARACTER(LEN=80)		      :: FMT
  CHARACTER(LEN=1024) 		  :: filename  
  
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
  REAL(iwp),SAVE,ALLOCATABLE:: fext_o_pp(:),timest(:)
  
  
!------------------------------------------------------------------------------
! 3. Start Program
!------------------------------------------------------------------------------
  IF(numpe .EQ. 1)PRINT*,"ParaFEM: "
  
  ! Barrier (may not be needed but safe)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
  
  ! Set Base paramenter
  argv		=  "Case"			  ! Name files write to
  nlen		=  4				    ! Length of Name
  limit		=  5000				  ! Max number of Interation in PCG
  tol		  =  1e-6				  ! Tolerance of PCG loop
  element	=  "hexahedron"	! Element Name
  
  ! Set Numerical and Material Values 
  alpha1  =  numVar(1)
  beta1   =  numVar(2)
  theta   =  numVar(3)
  dtim	  =  numVar(4)  

  c1	    =  (1._iwp-theta)*dtim
  c2	    =  beta1-c1
  c3	    =  alpha1+1._iwp/(theta*dtim);
  c4	    =  beta1+theta*dtim
  
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
   ALLOCATE(fext_o_pp(neq_pp))
   fext_o_pp  =  zero
  ENDIF
    
  timest=zero 
  timest(1)=elap_time()
  
  ! Clean Arrays
  x0_pp    =  zero;	  d1x0_pp	=  zero; 	x1_pp    =  zero; 
  vu_pp    =  zero; 	u_pp    =  zero;  d2x0_pp  =  zero; 
  loads_pp =  zero; 	d1x1_pp =  zero;	d2x1_pp  =  zero;
  d_pp	   =  zero; 	p_pp    =  zero;	x_pp	   =  zero; 
  xnew_pp  =  zero; 	
  
!------------------------------------------------------------------------------
! 4. Set Loads
!------------------------------------------------------------------------------

  timest(2)	=  elap_time()

  fext_pp	=  zero

  ! Load fext_pp based on global load vector
  CALL load(g_g_pp,g_num_pp,node,val,fext_pp)
 
  CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
 
  fext_pp=fext_pp+gravlo  ! Adding gravity if set

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

  u_pp		  =  zero
  vu_pp		  =  zero
  loads_pp	=  zero
  pmul_pp	  =  zero
  temp_pp	  =  zero

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

  loads_pp 	  =  fext_pp*theta*dtim+(1-theta)*dtim*fext_o_pp
  fext_o_pp 	=  fext_pp
  loads_pp	  =  u_pp+vu_pp+loads_pp
  temp_pp	    =  store_mm_pp*c3+store_km_pp*c4

  timest(6) =  elap_time()
  
!------------------------------------------------------------------------------
! 9. PCG
!------------------------------------------------------------------------------

  IF(ABS(SUM_P(loads_pp)) .GE. 1E-10)THEN
    IF(numpe .EQ. 1)PRINT*,"Solving using PCG"
    !CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
   
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
      
      up		=  DOT_PRODUCT_P(loads_pp,d_pp);
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
   
   IF(numpe .EQ. 1)PRINT*,"Number of Iterations: ",iters
   ENDIF
   x1_pp	  =  xnew_pp
   utemp_pp	=  zero
   d1x1_pp	=  (x1_pp-x0_pp)/(theta*dtim)-d1x0_pp*(1._iwp-theta)/theta
   d2x1_pp	=  (d1x1_pp-d1x0_pp)/(theta*dtim)-d2x0_pp*(1._iwp-theta)/theta
   x0_pp	  =  x1_pp;
   d1x0_pp	=  d1x1_pp
   d2x0_pp	=  d2x1_pp;
   
  timest(7) = elap_time()
  
!------------------------------------------------------------------------------
! 11. Gather Data from ntot,nels_pp to ndim,nodes_pp
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

