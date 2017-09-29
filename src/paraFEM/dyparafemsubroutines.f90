  !/****h* /dyparafemsubroutines
  !*  NAME
  !*    dyparafemsubroutines
  !*
  !*  SYNOPSIS
  !*    Usage:	This group of subroutines needs to me added to the
  !*		Make file of the top level solver (fsiFoam)
  !*
  !* 		Make/options, The following libraries need to be added
  !*
  !* 		"path to folder"/dyparafemsubroutines.o \
  !*		-L/"path To paraFEM"/lib -lParaFEM_mpi.5.0.3 \
  !*		-L/"path To paraFEM"/lib -larpack_linuxdesktop \
  !*		-lgfortran \
  !*		-l/"Path To Mpi"/lib -lmpi_f90 \
  !*
  !*	
  !*  FUNCTION
  !*	This is a group of subroutines required for the OpenFPCI.
  !*	These routines are based on the decomposition of program 12.9 
  !* 	found in "Programming the Finite Element Method".
  !*	P12.9 is the parallel analysis of the forced vibration of a
  !*	linear elastic solid.
  !* 	
  !*    
  !*  SUBROUTINE            PURPOSE
  !*
  !*    initparafem:		      Generates Initial Matricies and arrays
  !*    finddiagprecon:	    Finds Diagonal Preconditioner
  !*	  runparafem:		      Solves the governing equations  	
  !*	  checkparafem 		    Write Mesh and geometry to file (ensi)
  !*	  forcecheck 		      Write the loads to file (ensi)
  !*	  of2sg 			          OpenFOAM to Smith Griffiths format
  !*	  gloads			          Gravity Loading
  !*	  writeToFile		      Writes float field to file(debugging)
  !*	  writeToFileInt		    Writes int field to file(debugging)
  !*	  system_mem_usage	    Track Memory usage
  !*    populate_g_coord_pp2 Populate g_coord_pp array
  !*		 
  !*	
  !*  FUNCTION         	PURPOSE
  !* 
  !*	findnelspp		    Returns nels_pp to c++
  !*	findneqpp		      Returns neq_pp to c++
  !*
  !*  AUTHORS
  !* 	S.Hewitt
  !*  COPYRIGHT
  !*    (c) University of Manchester 1996-2017
  !******
  !*/

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------

  SUBROUTINE initparafem(g_coord,g_num_pp,rest,nn,nels,nr,sProp,g_g_pp,store_km_pp,store_mm_pp)
  !/****f*dyparafemsubroutines/initparafem
  !*  NAME
  !*    SUBROUTINE: initparafem
  !*
  !*  SYNOPSIS
  !*    Usage:  initparafem_(MeshData,g_num_OF,rest,&numPoints,	&
  !*			          &numCells,&nr,&e,&pois,steerG,stiff);   
  !*            
  !*  FUNCTION
  !*    Initialises ParaFEM, calculating the stiffness matirx [k], 
  !*	  Mass Matrix [M] and steering matrix (g_g_pp).
  !*
  !*  INPUTS
  !*    g_coord(ndim,nn)		  : Mesh coordinates (stressMesh OpenFOAM)		
  !*	  g_num_pp(nod,nels_pp)	: Steering Matrix (OF Format) 				
  !*	  rest(nr,nodof+1)		  : Restrained Nodes example:(node# x y z)
  !* 	  nn						        : # of Nodes
  !*	  nels					        : # of elements
  !*	  nr 						        : # of restrained Nodes
  !*    sProp					        : Solid Properties (e v rho) 	 
  !*	  			
  !*  OUTPUTS:
  !*    g_g_pp(ntot,nels_pp)			      : Global Steering Matrix
  !*   	store_km_pp(ntot,ntot,nels_pp)	: Stiffness Matrix [k]
  !*   	store_mm_pp(ntot,ntot,nels_pp)	: Mass Matrix [M]

  !*			
  !*  AUTHORS
  !*    S. Hewitt
  !*/
  !*  COMMENTS
  !*    neq,ntot are now global variables - must not be declared
  !* --------------------------------------------------------------------

  USE mpi_wrapper;    USE precision;  USE global_variables; 
  USE mp_interface;	  USE input;      USE output; 
  USE loading;        USE timing;     USE maths; 
  USE gather_scatter; USE steering;   USE new_library;

  IMPLICIT NONE

!------------------------------------------------------------------------------ 
! 1. Declare variables
!------------------------------------------------------------------------------

  INTEGER,PARAMETER		    ::nodof=3,ndim=3,nst=6,nod=8
  REAL(iwp),PARAMETER		  ::zero=0.0_iwp

  INTEGER,INTENT(IN) 		  ::nels,nn

  INTEGER,INTENT(INOUT)		::g_g_pp(ndim*nod,nels_pp),rest(nr,nodof+1)
  INTEGER,INTENT(INOUT)		::g_num_pp(nod,nels_pp)

  INTEGER                 ::iel,i,j,k,l,m,n,nr,nip,ndof,npes_pp,nlen
  INTEGER			            ::partitioner,printres,RSS,VM,RSSa,VMa

  REAL(iwp),INTENT(INOUT)	::g_coord(ndim,nn),sProp(3)
  REAL(iwp),INTENT(INOUT)	::store_km_pp(ndim*nod,ndim*nod,nels_pp)
  REAL(iwp),INTENT(INOUT)	::store_mm_pp(ndim*nod,ndim*nod,nels_pp)

  REAL(iwp)			          ::det,period,volume,tol,e,v,rho	

  LOGICAL			            ::converged=.false.
  LOGICAL			            ::consistent=.TRUE.
  LOGICAL			            ::initialised	

  CHARACTER(LEN=50)		    ::argv
  CHARACTER(LEN=15)		    ::element
  CHARACTER(LEN=6) 		    ::ch
  CHARACTER(LEN=1024) 		::filename  

!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------

  REAL(iwp),ALLOCATABLE		::points(:,:),dee(:,:),weights(:),timest(:)
  REAL(iwp),ALLOCATABLE		::jac(:,:),der(:,:),deriv(:,:),bee(:,:)
  REAL(iwp),ALLOCATABLE		::fun(:),emm(:,:),ecm(:,:),g_coord_pp(:,:,:)

  INTEGER,ALLOCATABLE		  ::node(:),localcount(:),readcount(:)
 
  CALL system_mem_usage(RSSa,VMa)
 
!------------------------------------------------------------------------------
! 3. Input and Initialisation
!------------------------------------------------------------------------------ 

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
  ! nels_pp calculated in OpenFOAM
  ! Calculate iel_start
  CALL setielstart() 

  ! Degrees of Freedon per Element
  ndof  =  nod*nodof
  ntot  =  ndof;
 
 !------------------------------------------------------------------------------
! 4. Populate g_coord_pp and g_num_pp
!------------------------------------------------------------------------------ 
  timest(2)=elap_time();
  
  ALLOCATE(g_coord_pp(nod,ndim,nels_pp))
  g_coord_pp  =  zero
  
  ! Poulate the Steering matrix
  !CALL POPULATE_G_NUM_PP(g_num,g_num_pp,npes,nod,nels)

  ! Convert from OpenFOAM to Smith Gritths format
  DO iel=1,nels_pp
    CALL of2sg(element,g_num_pp(:,iel),nod)
  ENDDO

  ! Populate Coorinate Matrix
  !CALL POPULATE_G_COORD_PP(g_coord,g_coord_pp,g_num_pp,npes,nn,nod,ndim)
  CALL POPULATE_G_COORD_PP2(g_coord,g_coord_pp,g_num_pp,nn,nod,ndim) 
 
  timest(3)=elap_time();
  
!------------------------------------------------------------------------------
! 5. Find the Steering arrays and equations per core
!------------------------------------------------------------------------------ 

  ALLOCATE(points(nip,ndim),dee(nst,nst),jac(ndim,ndim))
  ALLOCATE(der(ndim,nod),deriv(ndim,nod),bee(nst,ntot))
  ALLOCATE(weights(nip),ecm(ntot,ntot),emm(ntot,ntot),fun(nod))

  ! Rearrange the rest Array
  CALL rearrange(rest)

  ! Clean arrays
  g_g_pp  =  zero
  neq	  =  zero
  
  ! Find the global Steering Matrix
  elements_0: DO iel=1,nels_pp
  ! CALL find_g(g_num_pp(:,iel),g_g_pp(:,iel),rest) // Stable but slow
    CALL find_g3(g_num_pp(:,iel),g_g_pp(:,iel),rest)
  END DO elements_0
  
  ! Build GGL Array
  neq  =  MAXVAL(g_g_pp)
  neq  =  max_p(neq)
  CALL calc_neq_pp
  CALL calc_npes_pp(npes,npes_pp)
  CALL make_ggl(npes_pp,npes,g_g_pp)
  
  timest(4)=elap_time()	

!------------------------------------------------------------------------------
! 5. Element Stiffness and Mass Integration
!------------------------------------------------------------------------------ 
  
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
      store_km_pp(:,:,iel)=store_km_pp(:,:,iel) +                        & 
                          MATMUL(MATMUL(TRANSPOSE(bee),dee),bee) *      &
                          det*weights(i)
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
  
!------------------------------------------------------------------------------
! 6. Print Information about Runtime to File
!------------------------------------------------------------------------------ 
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

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

  SUBROUTINE finddiagparafem(store_km_pp,store_mm_pp,numVar,diag_precon_pp)
  !/****f*dyparafemsubroutines/finddiagparafem
  !*  NAME
  !*    SUBROUTINE: finddiagparafem
  !*
  !*  SYNOPSIS
  !*    Usage:      finddiagparafem_(stiff,precon); 
  !*
  !*  FUNCTION
  !*    Calculates the diagonal preconditioner matrix
  !*
  !*  INPUTS
  !*     store_km_pp(ntot,ntot,nels_pp)	: Stiffness Matrix	
  !*     store_mm_pp(ntot,ntot,nels_pp)	: Mass Matrix
  !*     numVar							: Array(alpha1 beta1 theta dTim) 			
  !*	 			
  !*  OUTPUTS:
  !*   	diag_precon_pp(neq_pp)		: Diagonal Preconditioner
  !*			
  !*  AUTHORS
  !*    S. Hewitt
  !*/
  !* -------------------------------------------------------------------
  
  USE mpi_wrapper;	  USE precision;	USE global_variables; 
  USE mp_interface; 	USE input;	    USE output; 
  USE loading; 		    USE timing; 	  USE maths; 
  USE gather_scatter;	USE steering; 	USE new_library; 
  
  IMPLICIT NONE
  
!------------------------------------------------------------------------------ 
! 1. Declare variables used in the main program
!------------------------------------------------------------------------------

  INTEGER,PARAMETER		    ::nodof=3,ndim=3,nst=6
  REAL(iwp),PARAMETER		  ::zero=0.0_iwp
  
  INTEGER			            ::i,j,k,iel,ndof

  REAL(iwp),INTENT(IN)		::store_km_pp(ndim*8,ndim*8,nels_pp)
  REAL(iwp),INTENT(IN)		::store_mm_pp(ndim*8,ndim*8,nels_pp)
  REAL(iwp),INTENT(IN)		::numVar(4)
  
  REAL(iwp),INTENT(INOUT)	::diag_precon_pp(neq_pp)
  
  REAL(iwp)			          ::alpha1,beta1,theta
  REAL(iwp)			          ::dtim,c3,c4
  
  REAL(iwp),ALLOCATABLE		::timest(:),diag_precon_tmp(:,:)
  
  ! Set Parameters
  alpha1  =  numVar(1)
  beta1   =  numVar(2)
  theta   =  numVar(3)
  dtim	  =  numVar(4)
  c3	    =  alpha1+1._iwp/(theta*dtim)
  c4	    =  beta1+theta*dtim
  ndof    =  ntot
  
!------------------------------------------------------------------------------ 
! 2. Calculate Diagonal preconditioner
!------------------------------------------------------------------------------
  ALLOCATE(diag_precon_tmp(ntot,nels_pp))
  
  diag_precon_pp  =  zero
  diag_precon_tmp  =  zero
  
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

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

  SUBROUTINE runparafem(numVar,val,node,loaded_nodes,time,nodes_pp, 	&
                      g_g_pp,g_num_pp,store_km_pp,store_mm_pp,diag_precon_pp,gravlo,Dfield,Ufield,Afield)
  
  !/****f*dyparafemsubroutines/runparafem
  !*  NAME
  !*    SUBROUTINE: runparafem
  !*
  !*  SYNOPSIS
  !*    Usage:
  !* 	
  !*  FUNCTION
  !*    Solve the governing equations and output the results 
  !*
  !*  INPUTS
  !*	  numVar				                  : (alpha1 beta1 theta dtim)
  !*	  val(ndim,loaded_nodes)		      : Force vector of loaded nodes
  !*   	node(loaded_nodes)		          : Node # of loaded_nodes
  !*  	loaded_nodes			              : # of loaded nodes		  
  !*  	Time			                      : Current time
  !*    nodes_pp                        : # of nodes per cores
  !*    g_g_pp(ntot,els)		            : Global Steering Matrix
  !*	  g_num_pp(nod,nels)		          : Steering Matrix
  !*    storkm_pp(ntot,ntot,nels_pp)	  : Stiffness Matrix [k]
  !*    stormm_pp(ntot,ntot,nels_pp)	  : Mass Matrix [M]
  !*	  diag_precon_pp(neq_pp)		      : Diagonal Preconditioner	
  !*	  gravlo(neq_pp)			            : vector of gravity loads	
  !*				
  !*  INPUTS/OUTPUTS:
  !*  	Dfield(ntot,nels_pp)		        : Nodal displacements
  !*  	Ufield(ntot,nels_pp)		        : Nodal velocities
  !*  	Afield(ntot,nels_pp)		        : Nodal accelerations
  !* 				
  !*  AUTHOR
  !*    S. Hewitt
  !*		
  !*/
  !* -------------------------------------------------------------------!
  
    
  USE mpi_wrapper;	  USE precision;	USE global_variables; 
  USE mp_interface; 	USE input;	    USE output; 
  USE loading; 		    USE timing; 	  USE maths; 
  USE gather_scatter;	USE steering; 	USE new_library; 
  
  IMPLICIT NONE
  
!------------------------------------------------------------------------------ 
! 1. Declare variables used in the main program
!------------------------------------------------------------------------------

  INTEGER,PARAMETER		      ::nodof=3,ndim=3,nst=6,nod=8
  REAL(iwp),PARAMETER		    ::zero=0.0_iwp,one=1.0_iwp
  
  INTEGER,INTENT(INOUT)		  ::loaded_nodes,node(loaded_nodes)
  INTEGER,INTENT(INOUT)		  ::g_num_pp(nod,nels_pp),nodes_pp
  INTEGER,INTENT(INOUT)		  ::g_g_pp(ntot,nels_pp)
  
  INTEGER			              ::iel,i,j,k,l,m,n,iters,printres
  INTEGER			              ::limit,nels,node_end,node_start
  INTEGER			              ::nlen,myCount,disps(npes),flag
  INTEGER			              ::nodesCount(npes),RSS,VM,RSSa,VMa  

  REAL(iwp),INTENT(INOUT) 	::diag_precon_pp(neq_pp)
  REAL(iwp),INTENT(INOUT) 	::store_km_pp(ntot,ntot,nels_pp)
  REAL(iwp),INTENT(INOUT) 	::store_mm_pp(ntot,ntot,nels_pp)
  REAL(iwp),INTENT(INOUT) 	::numVar(4),time,gravlo(neq_pp)
  REAL(iwp),INTENT(INOUT) 	::val(ndim,loaded_nodes)
  REAL(iwp),INTENT(INOUT)   ::Dfield(ntot,nels_pp),Ufield(ntot,nels_pp)
  REAL(iwp),INTENT(INOUT)   ::Afield(ntot,nels_pp)
  
  REAL(iwp)			            ::tol,up,alpha,beta,alpha1,beta1,theta,dtim
  REAL(iwp)			            ::c1,c2,c3,c4
  REAL(iwp)			            ::X,Y,Z
  
  LOGICAL			              ::converged
  CHARACTER(LEN=50)		      ::argv
  CHARACTER(LEN=15)		      ::element;
  CHARACTER(LEN=80)		      ::FMT
  CHARACTER(LEN=1024) 		  ::filename  
  
!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------
  REAL(iwp),SAVE,ALLOCATABLE::loads_pp(:),fext_pp(:),x1_pp(:),d1x1_pp(:)
  REAL(iwp),SAVE,ALLOCATABLE::d2x1_pp(:),x0_pp(:),d1x0_pp(:),d2x0_pp(:)
  REAL(iwp),SAVE,ALLOCATABLE::vu_pp(:),u_pp(:),p_pp(:),d_pp(:),x_pp(:)
  REAL(iwp),SAVE,ALLOCATABLE::xnew_pp(:),pmul_pp(:,:),utemp_pp(:,:)
  REAL(iwp),SAVE,ALLOCATABLE::temp(:),d1temp(:),d2temp(:),temp_pp(:,:,:)
  REAL(iwp),SAVE,ALLOCATABLE::disp_pp(:),vel_pp(:),acel_pp(:),eld_pp(:,:)
  REAL(iwp),SAVE,ALLOCATABLE::gDisp(:),gVel(:),gAcel(:)
  REAL(iwp),SAVE,ALLOCATABLE::fext_o_pp(:),timest(:)
  
  
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
  CALL system_mem_usage(RSS,VM)
  printres=0
  OPEN(10,FILE=argv(1:nlen)//".run",STATUS='REPLACE',ACTION='WRITE')
    WRITE(10,'(A,F10.4)') "Time to load:",timest(3)-timest(2)
    WRITE(10,'(A,F10.4)') "Time to Set ICs:",timest(4)-timest(3)
    WRITE(10,'(A,F10.4)') "Displacement & Velocity:",timest(5)-timest(4)    
    WRITE(10,'(A,F10.4)') "Time to Build RHS:",timest(6)-timest(5)
    WRITE(10,'(A,F10.4)') "Time to Solve using PCG:",timest(7)-timest(6)
    WRITE(10,'(A,F10.4)') "Time to Gather Data:",timest(8)-timest(7)
    WRITE(10,'(A,F10.4)') "Time in Routine(Total):",elap_time()-timest(1)
    WRITE(10,'(A,F10.4)') "-------------------------------------------------"
    WRITE(10,'(A,I10)') "Virtual Memory Change(Kb): ",VM-VMa
    WRITE(10,'(A,I10)') "RSS Memory Change(Kb): ",RSS-RSSa
    WRITE(10,'(A,I10)') "Total Virtual Memory(Kb): ",VM
    WRITE(10,'(A,I10)') "Total RSS Memory (Kb): ",RSS
    CLOSE(10)
  END IF
  timest(11)=elap_time()
  CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
  END SUBROUTINE
  !------------------------------------------------------------------------
  !------------------------------------------------------------------------
  !------------------------------------------------------------------------


  SUBROUTINE checkforce(force,sense,node,solidPatchIDSize,nn)
  !/****f*dyparafemsubroutines/checkforce
  !*  NAME
  !*    SUBROUTINE: checkforce
  !*
  !*  SYNOPSIS
  !*    Usage:      checkforce_(force,sense,node_ensi,&solidPatchIDSize,&
  !*		    paraFemSolidDisp,&numPoints);                                 
  !*
  !*  FUNCTION
  !*    Prints the forces to file in ENSI GOLD format  
  !*
  !*  INPUTS
  !*    force(loaded_nodes*ndim)	: Mesh coordinates (stressMesh OpenFOAM)
  !* 	  sense(loaded_nodes*ndim)	: Vector to define x,y,z 
  !* 	  node(loaded_nodes*ndim)		: Node Numbers of loaded nodes
  !*   	solidPatchIDSize			    : Size of boundary Mesh	
  !*    nn							          : # of Nodes  			
  !*
  !*  OUTPUTS:
  !*	Writes to file "argv".ensi.NDLDS
  !*	Contains information about loads on nodes
  !*			
  !*  AUTHORS
  !*    S. Hewitt
  !*    L. Margetts
  !*
  !* -------------------------------------------------------------------
  !USE mpi_wrapper  !remove comment for serial compilation
  USE precision;  USE global_variables; USE mp_interface; 
  USE input;      USE output;           USE loading;
  USE timing;     USE maths;            USE gather_scatter
  USE steering;   USE new_library;
  
  IMPLICIT NONE
  
  !----------------------  Declarations --------------------------------!
  INTEGER,PARAMETER       ::nodof=3,ndim=3,nst=6,nod=8
  INTEGER                 :: solidPatchIDSize,nlen,nn,loaded_nodes,pos,i,j
  INTEGER,INTENT(IN)      :: sense(solidPatchIDSize*ndim),node(solidPatchIDSize*ndim)
  REAL(iwp),INTENT(IN)    :: force(solidPatchIDSize*ndim)
  REAL(iwp),PARAMETER     ::zero=0.0_iwp
  REAL(iwp)               :: loads(nn*ndim)
  CHARACTER(LEN=15)       :: argv

  !--------------------- Set Variables for Mesh_ensi -------------------!
  argv="Case"; nlen=4; loads=zero

  !--------------------- Write loaded nodes ----------------------------!
  OPEN(16,FILE=argv(1:nlen)//'.ensi.NDLDS', status = "replace")
  WRITE(16,'(A)')     "Alya Ensight Gold --- Vector per-node variable file"
  WRITE(16,'(A/A/A)') "part", "      1","coordinates"
  DO i=1,solidPatchIDSize*ndim
    pos = ((node(i))*ndim)+sense(i)
    loads(pos)=force(i)
  END DO 

  DO j=1,ndim; DO i=1,nn
    WRITE(16,'(E12.5)') loads(((i-1)*ndim)+j)
  ENDDO; ENDDO

  CLOSE(16)
  END SUBROUTINE
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------


  SUBROUTINE checkparafem(g_coord,g_num,rest_ensi,nn,els)
  !*  NAME
  !*    SUBROUTINE: checkparafem
  !*  SYNOPSIS
  !*    Usage:      checkparafem_(MeshData,g_num_OF,rest_ensi,	&
  !*			&numPoints,&numCells);        
  !*  FUNCTION
  !*    Prints the Solid Mesh to mesh_ensi format  
  !*
  !*  INPUTS
  !*    g_coord(ndim,nn)	    : Mesh coordinates (stressMesh OpenFOAM)			
  !*	  g_num(nod,nels)		    : Steering Matrix (x-y-z)				
  !*	  rest_ensi(nodof+1,nr)	: Restrained Nodes
  !*    nn			              : # of Nodes 
  !*	  els=nels		          : # of Elements/Cells
  !*			
  !*  OUTPUTS:
  !*	Writes to file "argv".ensi.case and oter. "argv".ensi.xxx files
  !*	Contains information about Mesh, Material, Geometry and Restraints
  !*			
  !*  AUTHORS
  !*    S. Hewitt
  !*
  !* --------------------------------------------------------------------
  !USE mpi_wrapper  !remove comment for serial compilation
  USE precision;  USE global_variables; USE mp_interface; 
  USE input;      USE output;           USE loading;
  USE timing;     USE maths;            USE gather_scatter
  USE steering;   USE new_library;
  
  IMPLICIT NONE

  !----------------------  Declarations --------------------------!
  INTEGER,INTENT(IN)    :: nn,els

  !---------------------- Mesh_Ensi Declarations -----------------!
  INTEGER               :: i,j,k,l,m,n,iel,nels,ndim
  INTEGER               :: nod,BC,nlen,prnwidth		
	INTEGER               :: nstep,npri,remainder
  INTEGER               :: etype(els),rest_ensi(4,nn)
	INTEGER               :: mesh_print,g_num_print, g_num(8,els)
    
  REAL(iwp)             :: dtim,g_coord(3,nn)
  LOGICAL               :: solid 
  CHARACTER(LEN=15)     :: argv,element  
  !----------------------------------------------------------------------
  ! 1. Initialisation
  !----------------------------------------------------------------------
    !nn   = UBOUND(g_coord,2) ; 
    ndim = UBOUND(g_coord,1)
    nels = UBOUND(g_num,2)   ; nod  = UBOUND(g_num,1)

!--------------------- Set Variables for Mesh_ensi --------------!
    etype(:)  = 0
    argv      = "Case"
    nlen      = 4
    element   = "hexahedron"
    nstep     = 1
    npri      = 1
    dtim      = 1
    solid     = .true.
  !------------------------------------------------------------------------------
  ! 2. Write case file
  !------------------------------------------------------------------------------
  
    OPEN(12,FILE=argv(1:nlen)//'.ensi.case')
  
    WRITE(12,'(A/A)')    "#", "# Post-processing file generated by subroutine &
                               &WRITE_ENSI in "
    WRITE(12,'(A,A,/A)') "#", " Smith, Griffiths and Margetts, 'Programming the &
                               &Finite Element Method',","# Wiley, 2013."        
    WRITE(12,'(A/A/A)')  "#","# Ensight Gold Format","#"
    WRITE(12,'(2A/A)')   "# Problem name: ",argv(1:nlen),"#"
    WRITE(12,'(A/A/A)')  "FORMAT","type:  ensight gold","GEOMETRY"
    WRITE(12,'(2A/A)')   "model: 1  ",argv(1:nlen)//'.ensi.geo',"VARIABLE"
    WRITE(12,'(2A)')     "scalar per element:  material      ",                &
                          argv(1:nlen)//'.ensi.MATID'
    IF(solid) THEN
      WRITE(12,'(2A)')   "scalar per node:     restraint     ",                &
                          argv(1:nlen)//'.ensi.NDBND'
      WRITE(12,'(2A)')   "vector per node:     displacement  ",                &
                          argv(1:nlen)//'.ensi.DISPL-******'
    ELSE
      WRITE(12,'(2A)')   "scalar per node:     pressure      ",                &
                          argv(1:nlen)//'.ensi.PRESSURE-******'
    END IF
    WRITE(12,'(2A)')     "vector per node:     load          ",                &
                          argv(1:nlen)//'.ensi.NDLDS'
    WRITE(12,'(A/A)')     "TIME","time set:     1"
    WRITE(12,'(A,I5)')    "number of steps:",nstep/npri
    WRITE(12,'(A,I5)')    "filename start number:",npri
    WRITE(12,'(A,I5)')    "filename increment:",npri
    WRITE(12,'(A)')       "time values:"
    prnwidth  = 5
    remainder = mod(nstep/npri,prnwidth)
    n         = ((nstep/npri) - remainder)/prnwidth
    IF(nstep/npri<=prnwidth) THEN
      DO i=1,nstep,npri
        IF(i==nstep) THEN
          WRITE(12,'(E12.5)') i*dtim
        ELSE
          WRITE(12,'(E12.5)',ADVANCE='no') i*dtim
        END IF
      END DO
    ELSE
      IF(remainder==0) THEN
        DO j=1,n
          m = ((j-1)*prnwidth)+1
          l = ((j-1)*prnwidth)+prnwidth
          WRITE(12,'(5E12.5)') (k*dtim,k=m,l)
        END DO
      ELSE
  !     DO j=1,n-1
        DO j=1,n
          m = ((j-1)*prnwidth)+1
          l = ((j-1)*prnwidth)+prnwidth
          WRITE(12,'(5E12.5)') (k*dtim,k=m,l)
        END DO
        m = (n*prnwidth)+1
        l = (n*prnwidth)+remainder
        DO i=m,l
          IF(i==l) THEN
            WRITE(12,'(E12.5)') dtim*i
          ELSE
            WRITE(12,'(E12.5)',ADVANCE='no') dtim*i
          END IF
        END DO
      END IF
    END IF
   
    CLOSE(12)
  
  !------------------------------------------------------------------------------
  ! 3. Write geometry file
  !------------------------------------------------------------------------------
  
    OPEN(13,FILE=argv(1:nlen)//'.ensi.geo')
    WRITE(13,'(/2A)')   "Problem name: ", argv(1:nlen)
    WRITE(13,'(A/A/A)') "Geometry files","node id given","element id given"
    WRITE(13,'(A/A)')   "part","      1"
    IF(ndim==2) WRITE(13,'(A)') "2d-mesh"
    IF(ndim==3) WRITE(13,'(A)') "Volume Mesh"
    WRITE(13,'(A)')     "coordinates"
    
    WRITE(13,'(I10)') nn
    DO j=1,ndim
      DO i=1,nn  
        WRITE(13,'(E12.5)') g_coord(j,i)
      END DO
    END DO
  
    IF(ndim==2) THEN ! ensight requires zeros for the z-ordinate
      DO i=1,nn
        WRITE(13,'(A)') " 0.00000E+00"
      END DO
    END IF
  
    SELECT CASE(element)
      CASE('triangle')
        SELECT CASE(nod)
          CASE(3)
            WRITE(13,'(A/I10)') "tria3", nels
            DO i = 1,nels
              WRITE(13,'(3I10)')g_num(3,i),g_num(2,i),g_num(1,i)
            END DO
          CASE DEFAULT
            WRITE(13,'(A)')   "# Element type not recognised"
        END SELECT
      CASE('quadrilateral')
        SELECT CASE(nod)
          CASE(4)
            WRITE(13,'(A/I10)') "quad4", nels
            DO i = 1,nels
              WRITE(13,'(4I10)')g_num(1,i),g_num(4,i),g_num(3,i),g_num(2,i)
            END DO
          CASE(8)
            WRITE(13,'(A/I10)') "quad8", nels
            DO i = 1,nels
              WRITE(13,'(8I10)')g_num(1,i),g_num(7,i),g_num(5,i),g_num(3,i),    &
                                g_num(8,i),g_num(6,i),g_num(4,i),g_num(2,i)
            END DO
          CASE DEFAULT
            WRITE(13,'(A)')   "# Element type not recognised"
        END SELECT
      CASE('hexahedron')
        SELECT CASE(nod)
          CASE(8)
            WRITE(13,'(A/I10)') "hexa8", nels
            DO i = 1,nels
	    WRITE(13,'(8I10)') g_num(1,i),g_num(2,i),g_num(6,i),g_num(5,i),   &
                               g_num(3,i),g_num(4,i),g_num(8,i),g_num(7,i)
            END DO
          CASE(20)
! COMMENTED OUT BECAUSE THE WARNINGS WERE BUGGING ME
!            WRITE(13,'(A/I10)') "hexa20", nels
!            DO i = 1,nels
!              WRITE(13,'(20I10)')                                               &
!                g_num(1,i), g_num(7,i), g_num(19,i),g_num(13,i),g_num(3,i),     &
!                g_num(5,i), g_num(17,i),g_num(15,i),g_num(8,i), g_num(12,i),    &
!                g_num(20,i),g_num(9,i), g_num(4,i), g_num(11,i),g_num(16,i),    &
!                g_num(10,i),g_num(2,i), g_num(6,i), g_num(18,i),g_num(14,i) 
!            END DO
          CASE DEFAULT
            WRITE(13,'(A)')   "# Element type not recognised"
        END SELECT
      CASE('tetrahedron')
        SELECT CASE(nod)
          CASE(4)
            WRITE(13,'(A/I10)') "tetra4", nels
            DO i = 1,nels
              WRITE(13,'(4I10)') g_num(1,i),g_num(3,i),g_num(2,i),g_num(4,i)
            END DO
          CASE DEFAULT
            WRITE(13,'(A)')   "# Element type not recognised"
        END SELECT
      CASE DEFAULT
        WRITE(13,'(A)')       "# Element type not recognised"
    END SELECT
  
    CLOSE(13)   
  
  !------------------------------------------------------------------------------
  ! 4. Write file containing material IDs
  !------------------------------------------------------------------------------
  
    OPEN(14,FILE=argv(1:nlen)//'.ensi.MATID')
    WRITE(14,'(A)') "Alya Ensight Gold --- Scalar per-element variable file"
    WRITE(14,'(A/A)') "part", "      1"
  
    SELECT CASE(element)
      CASE('triangle')
        SELECT CASE(nod) 
          CASE(3)
            WRITE(14,'(A)') "tria3"
          CASE DEFAULT
            WRITE(14,'(A)') "# Element type not recognised"
        END SELECT
      CASE('quadrilateral')
        SELECT CASE(nod) 
          CASE(4)
            WRITE(14,'(A)') "quad4"
          CASE(8)
            WRITE(14,'(A)') "quad8"
          CASE DEFAULT
            WRITE(14,'(A)') "# Element type not recognised"
        END SELECT
      CASE('hexahedron')
        SELECT CASE(nod) 
          CASE(8)
            WRITE(14,'(A)') "hexa8"
          CASE(20)
            WRITE(14,'(A)') "hexa20"
          CASE DEFAULT
            WRITE(14,'(A)') "# Element type not recognised"
        END SELECT
      CASE('tetrahedron')
        SELECT CASE(nod)
          CASE(4)
            WRITE(14,'(A)') "tetra4"
          CASE DEFAULT
          WRITE(14,'(A)') "# Element type not recognised"
        END SELECT
      CASE DEFAULT
        WRITE(14,'(A)')   "# Element type not recognised"
    END SELECT
   
    DO i=1,nels; WRITE(14,'(I10)') etype(i); END DO
  
    WRITE(14,'(A)')
  
    CLOSE(14)
  
  !------------------------------------------------------------------------------
  ! 5. Write boundary conditions. Encoded using formula: 4z + 2y + 1x
  !
  !    110 = 1   010 = 2   100 = 3   011 = 4   101 = 5   001 = 6   000 = 7
  !------------------------------------------------------------------------------
  
    IF(solid) THEN
      OPEN(15,FILE=argv(1:nlen)//'.ensi.NDBND')
      WRITE(15,'(A)')     "Alya Ensight Gold --- Scalar per-node variable file"
      WRITE(15,'(A/A/A)') "part", "      1","coordinates"
      IF(ndim==3) THEN
	DO i=1,UBOUND(g_coord,2)
	   BC = (rest_ensi(4,i)*4) + (rest_ensi(3,i)*2) + rest_ensi(2,i)
          WRITE(15,'(I2)') BC
	ENDDO 
      ELSE IF(ndim==2) THEN
        DO i=1,nn
	! WRITE 2D CASE
        END DO
      ELSE
        PRINT *, "Wrong number of dimensions in mesh_ensi"
      END IF   
    END IF
  
    CLOSE(15)
  
    RETURN
  
  END SUBROUTINE

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------


  SUBROUTINE POPULATE_G_COORD_PP(g_coord,g_coord_pp,g_num_pp,npes,nn,nod,ndim)
  !*  NAME
  !*    SUBROUTINE: populate_g_coord_pp
  !*  SYNOPSIS
  !*    Usage:      CALL populate_g_coord_pp(g_coord,g_coord_pp,g_num_pp,	&
  !* 					     npes,nn,nod,ndim)
  !*  FUNCTION
  !*	Populates g_coord_pp based on g_coord

  !*  INPUTS
  !*    g_coord(ndim,nn)	: Mesh coordinates (stressMesh OpenFOAM)		
  !*	  g_num(nod,nels)		: Steering Matrix (OF Format) 		
  !*	  nn			          : Number of solid Points 
  !*	  npes			        : Number of processes
  !*	  nod			          : Number of nodes per element
  !*	  ndim			        : Number of dimensions per node		  
  !*	  g_num_pp(nod,nels): Steering Matrix
  !*			
  !*  OUTPUTS:
  !*	  g_coord_pp(nod,ndim,nels_pp)	: Coordinate matrix
  !*    
  !*  AUTHOR
  !*    Sam Hewitt
  !*  CREATION DATE
  !*    28th September 2016
  !*  COPYRIGHT
  !*    (c) University of Manchester 2007-2016
  !******
  !* G_COORD_PP exists on all processors

  USE precision; USE global_variables;
  IMPLICIT NONE

  INTEGER, INTENT(IN)          :: nn,npes,nod,ndim
  INTEGER, INTENT(IN)          :: g_num_pp(nod,nels_pp)
  INTEGER                      :: nnStart         ! first node ID in g_coord
  INTEGER                      :: nnEnd           ! last node ID in g_coord
  INTEGER                      :: iel, i, j, k, l ! loop counters
  INTEGER                      :: readSteps
  INTEGER                      :: readCount
  INTEGER                      :: readRemainder    
  REAL(iwp), INTENT(INOUT)     :: g_coord_pp(nod,ndim,nels_pp) 
  REAL(iwp), INTENT(IN)        :: g_coord(ndim,nn)    
  REAL(iwp)                    :: zero = 0.0_iwp

  !------------------------------------------------------------------------------
  ! 1. Find READSTEPS, the number of steps in which the read will be carried
  !    out, READCOUNT, the size of each read and READREMAINDER, the number
  !    of entries to be read after the last READSTEP.
  !------------------------------------------------------------------------------

  IF(npes > nn) THEN
    readSteps     = 1
    readRemainder = 0
    readCount     = nn
  ELSE
    readSteps     = npes
    readRemainder = MOD(nn,readSteps)
    readCount     = (nn - readRemainder) / npes
  END IF

  !------------------------------------------------------------------------------
  ! 2. Go round READSTEPS loop, read data, broadcast and populate g_coord_pp
  !------------------------------------------------------------------------------
  DO i = 1, readSteps
    nnStart = (i-1) * readCount + 1
    nnEnd   =  i    * readCount
    DO iel = 1, nels_pp
      DO k = 1, nod
        IF(g_num_pp(k,iel) < nnStart) CYCLE
        IF(g_num_pp(k,iel) > nnEnd)   CYCLE
        l                  = g_num_pp(k,iel) !- nnStart + 1
        g_coord_pp(k,:,iel)= g_coord(:,l)
      END DO
    END DO
  END DO

  !------------------------------------------------------------------------------
  ! 3. If READREMAINDER > 0, collect remaining entries
  !------------------------------------------------------------------------------

  IF(readRemainder > 0) THEN
    bufsize  = ndim * readRemainder
    nnStart  = (readSteps * readCount) + 1
    nnEnd    = nnStart + readRemainder - 1

    IF(nnEnd > nn) THEN
      PRINT *, "Too many nodes"
      CALL closefem()
      STOP
    END IF

    DO iel = 1, nels_pp
      DO k = 1, nod
        IF(g_num_pp(k,iel) < nnStart) CYCLE
        IF(g_num_pp(k,iel) > nnEnd)   CYCLE
        l                  = g_num_pp(k,iel) !- nnStart + 1
        g_coord_pp(k,:,iel)= g_coord(:,l)
      END DO
    END DO
  END IF
  END SUBROUTINE


  SUBROUTINE POPULATE_G_NUM_PP(g_num,g_num_pp,npes,nod,nels)
  !*  NAME
  !*    SUBROUTINE: populate_g_num_pp
  !*  SYNOPSIS
  !*    Usage:      CALL populate_g_num_pp(g_num,g_num_pp,npes,nod,nels)
  !*                                     
  !*  FUNCTION
  !*    
  !*  INPUTS
  !*
  !*  AUTHOR
  !*    Sam Hewitt
  !*  CREATION DATE
  !*    28th September 2016
  !*  COPYRIGHT
  !*    (c) University of Manchester 2007-2016
  !******
  !* G_NUM exists on all processors may need to be changed.

  USE precision; USE global_variables; USE mp_interface; USE mpi_wrapper
  IMPLICIT NONE

  INTEGER, INTENT(IN)           :: npes, nels ,nod,g_num(nod,nels)
  INTEGER, INTENT(INOUT)        :: g_num_pp(nod,nels_pp)
  INTEGER                       :: i, k, ielstart ,ielend
  INTEGER, ALLOCATABLE          :: localCount(:),readCount(:)

  ALLOCATE(localcount(npes),readCount(npes))
  readCount = 0; localCount = 0; 
  localCount(numpe) = nels_pp
  CALL MPI_ALLREDUCE(localCount,readCount,npes,MPI_INTEGER,MPI_SUM,	&
                    MPI_COMM_WORLD,ier)
  DO i=2,npes
    readcount(i)=readcount(i)+readcount(i-1)
  END DO
  IF(numpe==1)THEN
   ielstart = 1;
  ELSE
   ielstart = readcount(numpe-1)+1
  END IF  
  ielend = readcount(numpe)
  g_num_pp(:,:)=g_num(:,ielstart:ielend);
  END SUBROUTINE


SUBROUTINE POPULATE_G_COORD_PP2(g_coord,g_coord_pp,g_num_pp,nn,nod,ndim)
  !*  NAME
  !*    SUBROUTINE: populate_g_coord_pp2
  !*  SYNOPSIS
  !*    Usage: CALL populate_g_coord_pp2(g_coord,g_coord_pp,g_num_pp,	&
  !* 					     nod,ndim)
  !*  FUNCTION
  !*	  Populates g_coord_pp based on g_coord

  !*  INPUTS
  !*    g_coord(ndim,nn)	: Mesh coordinate list		
  !*	  g_num(nod,nels)		: Steering Matrix (OF Format) 		
  !*	  nod			          : Number of nodes per element
  !*	  ndim			        : Number of dimensions per node		  
  !*	  g_num_pp(nod,nels): Steering Matrix
  !*			
  !*  OUTPUTS:
  !*	  g_coord_pp(nod,ndim,nels_pp)	: Coordinate matrix
  !*    
  !*  AUTHOR
  !*    Sam Hewitt
  !*  CREATION DATE
  !*    6th May 2017
  !*  COPYRIGHT
  !*    (c) University of Manchester 2007-2016
  !******

  USE precision; USE global_variables;
  IMPLICIT NONE

  INTEGER, INTENT(IN)          :: nn,nod,ndim
  INTEGER, INTENT(IN)          :: g_num_pp(nod,nels_pp)
  INTEGER                      :: iel, i, j, k, l ! loop counters 
  REAL(iwp), INTENT(INOUT)     :: g_coord_pp(nod,ndim,nels_pp) 
  REAL(iwp), INTENT(IN)        :: g_coord(ndim,nn)    
  REAL(iwp)                    :: zero = 0.0_iwp
  
  ! Simple Method used after OpenFOAM Decomposition

  DO iel=1,nels_pp
    DO j=1,nod 
	k = g_num_pp(j,iel)
	g_coord_pp(j,:,iel)=g_coord(:,k)
    ENDDO
  ENDDO
  
  END SUBROUTINE


  SUBROUTINE closefem()
   USE mp_interface
   IMPLICIT NONE
   CALL SHUTDOWN()	
  END SUBROUTINE closefem


  SUBROUTINE of2sg(element,vector,nod)
  !*  NAME
  !*    SUBROUTINE: of2sg
  !*  SYNOPSIS
  !*    Usage:      of2sg(element,g_num_pp(:,iel))       
  !*  FUNCTION
  !*    rearrange the OpenFOAM steering matrix to ParaFEM format  
  !*
  !*  INPUTS
  !*    element			: Element Type			
  !*	temp (INOUT)			: Vector				
  !*			
  !*  AUTHORS
  !*    S. Hewitt
  !* 
  !* ---------------------------------------------------------------------------

  USE precision; IMPLICIT NONE
  INTEGER               :: i,temp(nod)
  INTEGER,INTENT(IN)    :: nod
  INTEGER,INTENT(INOUT) :: vector(nod)
  CHARACTER(LEN=15)     :: element 

  temp=0.0_iwp

  SELECT CASE(element) 
   CASE('hexahedron')
     SELECT CASE(nod)
       CASE(8)
	! OpenFOAM to Abaqus Format
!            temp(:)=vector(:)
! 	    vector(1)=temp(5)
!  	    vector(2)=temp(6)
!  	    vector(3)=temp(2)
!  	    vector(4)=temp(1)
! 	    vector(5)=temp(7)
!	    vector(6)=temp(8)
!	    vector(7)=temp(4)
!	    vector(8)=temp(3)

	! Code from OpenFOAM2sg:
      	temp(:)      = vector(:)
        vector(1)  = temp(5)
        vector(2)  = temp(7)
        vector(3)  = temp(8)
        vector(4)  = temp(6)
        vector(5)  = temp(1)
        vector(6)  = temp(3)
        vector(7)  = temp(4)
        vector(8)  = temp(2)
    END SELECT
  END SELECT
  END SUBROUTINE of2sg

 !* ---------------------------------------------------------------------------
 !*  FUNCTIONS: 
 !*	PURPOSE: To Pass variables defined in global_varaiables.mod
 !*		 to OpenFOAM
 !*  AUTHORS
 !*    S. Hewitt
 !* ---------------------------------------------------------------------------
 INTEGER FUNCTION findnelspp()
     USE precision; USE global_variables;
     IMPLICIT NONE	
     findnelspp=nels_pp
     RETURN
     END FUNCTION


 INTEGER FUNCTION calcnelsppof(nels,npes)
     USE precision; USE global_variables;
     IMPLICIT NONE

    INTEGER, INTENT(IN)	:: nels,npes
    INTEGER             :: num_nels_pp1, nels_pp1, nels_pp2,iel_start
	
     IF (npes == 1) THEN
        nels_pp1 = 0
        nels_pp2 = nels
        nels_pp  = nels_pp2
        iel_start = 1
      ELSE
        nels_pp2     = nels/npes
        num_nels_pp1 = nels - nels_pp2*npes
        IF (num_nels_pp1 == 0) THEN
          nels_pp1   = nels_pp2
        ELSE
          nels_pp1   = nels_pp2 + 1
        ENDIF
        IF (numpe <= num_nels_pp1 .OR. num_nels_pp1 == 0) THEN
          nels_pp    = nels_pp1
          iel_start  = (numpe - 1)*nels_pp1 + 1
        ELSE
          nels_pp    = nels_pp2
          iel_start  = num_nels_pp1*nels_pp1+                                 &
                       (numpe-num_nels_pp1-1)*(nels_pp1-1)+1
        ENDIF
      ENDIF
     calcnelsppof=1
     RETURN
     END FUNCTION

  INTEGER FUNCTION setnelspp(nCells)
   USE precision; USE global_variables;
   IMPLICIT NONE
   INTEGER,INTENT(IN) :: nCells	
   nels_pp=nCells
   RETURN
  END FUNCTION

  SUBROUTINE setielstart()
  USE precision;	USE global_variables;
  USE mp_interface; USE mpi_wrapper
  USE gather_scatter 
  IMPLICIT NONE

  INTEGER, ALLOCATABLE             :: psize(:)
  INTEGER                          :: i,p

  ALLOCATE(psize(npes))
  psize = 0
    
 CALL MPI_ALLGATHER(nels_pp,1,MPI_INT,psize,1,MPI_INT,MPI_COMM_WORLD,ier)
 iel_start = 0

  IF(numpe==1) THEN
    iel_start = 1
  ELSE
    DO i = 2, numpe
      iel_start = iel_start + psize(i-1)
    END DO
    iel_start   = iel_start + 1
  END IF 

  WRITE(details,'(A,I8,A,I8)') 'PE no: ', numpe, ' iel_start: ', iel_start

  DEALLOCATE(psize)

  END SUBROUTINE setielstart


  INTEGER FUNCTION findneqpp()
   USE precision; USE global_variables;
   IMPLICIT NONE	
   findneqpp=neq_pp
   RETURN
  END FUNCTION

  SUBROUTINE gloads(gravlo,specWeight,nn,nodof,nod,ndim,nr,g_coord,g_num_pp,rest)
  !*  NAME
  !*    SUBROUTINE: GLOADS
  !*  SYNOPSIS
  !*    Usage:      GLOADS(element,g_num_pp(:,iel))       
  !*  FUNCTION
  !*    Supplies a gravitational force vector   
  !*
  !*  INPUTS/OUTPUTS
  !*    specWeight		: value of specWeight			
  !*	gravlo (OUT)		: Vector of gravity loads
  !*	nf			: Gloabal freedom Matrix	
  !*	g_num			: Element Steering matrix
  !*	g_coord_pp		: Element cooridnate matrix
  !*	element			: Element Type	
  !*			
  !*  AUTHORS
  !*    S. Hewitt
  !* 
  !*  WARNING
  !*	REST must be rearranged before function call
  !* ---------------------------------------------------------------------------
  USE mpi_wrapper  !remove comment for serial compilation#
  USE precision; USE global_variables; USE mp_interface; USE input
  USE output; USE loading; USE timing; USE maths; USE gather_scatter
  USE steering; USE new_library; IMPLICIT NONE

  INTEGER :: i,j,k,iel,ndof,nodof,nn,nip,nf(nodof,nn),node_end,node_start,nodes_pp
  INTEGER,INTENT(IN) :: nr,nod,ndim,g_num_pp(nod,nels_pp),rest(nr,nodof+1)
  REAL(iwp),PARAMETER::zero=0.0_iwp
  REAL(iwp)::lamba,det,g_coord_pp(nod,ndim,nels_pp),g_coord(ndim,nn),		&
            pmul_pp(ntot,nels_pp),load_pp(ndim*nn),gravlotmp(0:neq_pp)
  REAL(iwp),INTENT(IN)::specWeight
  REAL(iwp),INTENT(INOUT)::gravlo(neq_pp)
  CHARACTER(LEN=15) :: element

  INTEGER,ALLOCATABLE :: num(:),g(:)
  REAL(iwp),ALLOCATABLE :: g_g(:,:),points(:,:),weights(:),eld(:),fun(:),	&
	der(:,:),jac(:,:)

  !----------  find the steering array and equations per process ---------
  element="hexahedron";ndof=nodof*nod;nip=8;nf=zero;
  IF(numpe.NE.1)THEN
      PRINT*,"--------WARNING--------"
      PRINT*," GLOADS is Serial Only "
  ELSE
      PRINT*,"Calculating gravity loads:"
      PRINT*,"Element :",element
      PRINT*,"Number of Integration points:",nip	
  ENDIF
  CALL REST_TO_NF(rest,nf)
  CALL POPULATE_G_COORD_PP2(g_coord,g_coord_pp,g_num_pp,nn,nod,ndim)
  ALLOCATE(num(nod),g(ntot),g_g(ntot,nels_pp),points(nip,ndim),weights(nip),	&
	fun(nod),eld(ntot),der(ndim,nod),jac(ndim,ndim))
  DO iel=1,nels_pp
    num=g_num_pp(:,iel); CALL num_to_g(num,nf,g); g_g(:,iel)=g
  END DO
  CALL sample(element,points,weights); gravlotmp=zero; 
 
  elements_2: DO iel=1,nels_pp  
   g=g_g(:,iel); eld=zero;
   gauss_points_1: DO i=1,nip     
     CALL shape_der(der,points,i); jac=MATMUL(der,g_coord_pp(:,:,iel))
     det=determinant(jac);
     CALL shape_fun(fun,points,i)
     eld(2:ndof-1:3)=eld(2:ndof-1:3)+fun(:)*det*weights(i)
   END DO gauss_points_1   
   gravlotmp(g)=gravlotmp(g)-eld*specWeight
  END DO elements_2

  ! gravlo decalred from 0:neq_pp - loop to correct it
  DO i=1,neq_pp
     gravlo(i)=gravlotmp(i)
  ENDDO

  !pmul_pp = zero;
  !CALL calc_nodes_pp(nn,npes,numpe,node_end,node_start,nodes_pp) 
  !CALL gather(gravlo(:),pmul_pp)
  !CALL scatter_nodes(npes,nn,nels_pp,g_num_pp,nod,ndim,nodes_pp,       	&
   !                      node_start,node_end,pmul_pp,load_pp,1)

  !OPEN(16,FILE='Case.ensi.NDLDS', status = "replace")
  !WRITE(16,'(A)')     "Alya Ensight Gold --- Vector per-node variable file"
  !WRITE(16,'(A/A/A)') "part", "      1","coordinates"
  !DO j=1,ndim; DO i=1,nn
  ! WRITE(16,'(E12.5)') load_pp(((i-1)*ndim)+j)
  !ENDDO; ENDDO
  !CLOSE(16)
 
  END SUBROUTINE gloads

   SUBROUTINE writeToFile(fileNum,fileName,field,dim1)
    USE mpi_wrapper  !remove comment for serial compilation#
    USE precision; USE global_variables; USE mp_interface; USE input
    USE output; USE loading; USE timing; USE maths; USE gather_scatter
    USE steering; USE new_library; IMPLICIT NONE
    INTEGER::i
    INTEGER,INTENT(IN)::dim1,fileNum
    REAL(iwp),INTENT(INOUT)::field(dim1)
    CHARACTER(LEN=*),INTENT(IN)::fileName  
    OPEN(fileNum,FILE=fileName,status = "replace")
     DO i=1,dim1
      WRITE(fileNum,'(E12.5)') field(i)
     ENDDO
    CLOSE(fileNum)
  END SUBROUTINE writeToFile


  SUBROUTINE writeToFileInt(fileNum,fileName,field,dim1)
    USE mpi_wrapper  !remove comment for serial compilation#
    USE precision; USE global_variables; USE mp_interface; USE input
    USE output; USE loading; USE timing; USE maths; USE gather_scatter
    USE steering; USE new_library; IMPLICIT NONE
    INTEGER::i
    INTEGER,INTENT(IN)::dim1,fileNum
    INTEGER,INTENT(INOUT)::field(dim1)
    CHARACTER(LEN=*),INTENT(IN)::fileName
    OPEN(fileNum,FILE=fileName,status = "replace")
     DO i=1,dim1
      WRITE(fileNum,*) field(i)
     ENDDO
    CLOSE(fileNum)
  END SUBROUTINE writeToFileInt

  SUBROUTINE system_mem_usage(valueRSS,valueVM)
  !use ifport !if on intel compiler
  IMPLICIT NONE
  INTEGER, INTENT(out) :: valueRSS,valueVM
  CHARACTER(len=200):: filename=' '
  CHARACTER(len=80) :: line
  CHARACTER(len=8)  :: pid_char=' '
  INTEGER :: pid
  LOGICAL :: ifxst
  valueRSS=-1; valueVM =-1 ! return negative number if not found
  !--- get process ID
  pid=getpid()
  WRITE(pid_char,'(I8)') pid
  filename='/proc/'//trim(adjustl(pid_char))//'/status'
  !--- read system file
  INQUIRE (FILE=filename,EXIST=ifxst)
  IF (.NOT.ifxst) THEN
   WRITE (*,*) 'system file does not exist'
   RETURN
  ENDIF
  OPEN(UNIT=100, FILE=filename, ACTION='read')
  DO
  READ (100,'(a)',END=120) line
  ! VmSize printed first in proc/$pid/status file
  IF (line(1:7).EQ.'VmSize:') THEN
     READ (line(8:),*) valueVM
    ELSE IF (line(1:6).EQ.'VmRSS:')THEN
     READ (line(7:),*) valueRSS
     EXIT
   ENDIF
  ENDDO
  120 CONTINUE
  CLOSE(100)
  RETURN
  END SUBROUTINE system_mem_usage 

