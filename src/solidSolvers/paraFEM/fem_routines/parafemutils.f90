  !/****h* parafemutils
  !*  NAME
  !*    parafemutils - Collection of fortran subroutines                          
  !*
  !*  SYNOPSIS
  !*    This group of subroutines is a collection of routines
  !*    that is used by the parafem solvers.
  !*    
  !*    These routines are separated from the solver 
  !*    fortran files (parafeml.f90 and parafemnl.f90) 
  !*    to improve the structure of the software.
  !*
  !*  FUNCTION
  !*    This is a group of subroutines required for the OpenFPCI.
  !*    
  !*    The routines provide an interface between Foam-Extend and 
  !*    ParaFEM.
  !* 
  !*    Subroutine           Purpose    
  !*
  !*    finddiagprecon       Generate K,M and the diagonal preconditioner
  !*    checkparafem         Writes mesh and geometry to file (ensi)
  !*    checkforce           Writes the loads to file (ensi)
  !*    of2sg                Foam-Ex to Smith-Griffiths mesh formats
  !*    gloads               Applies gravity loads
  !*   
  !*    populate_g_coord_pp  Populate g_coord_pp array
  !*    populate_g_coord_pp2 Populate g_coord_pp array
  !*    populate_g_num_pp    Populate g_num_pp array
  !*
  !*    writeToFile          Writes float field to file
  !*    writeToFileInt       Writes int field to file
  !*    write_largestrain    Write out the data for femlargestrain.f90
  !*
  !*    system_mem_usage     Track Memory usage
  !* 
  !*    Functions            Purpose       
  !*
  !*    findnelspp           Returns nels_pp to c++
  !*    findneqpp            Returns neq_pp to c++
  !*    setielstart          Set element start number
  !*    setnelspp            Set number of elements/processor
  !*    calcnelsppof         Calculate number of elements/processor
  !*
  !*  AUTHOR
  !* 	  S.Hewitt
  !*  COPYRIGHT
  !*    (c) University of Manchester 1996-2017
  !******
  !*/

  SUBROUTINE finddiagprecon(store_km_pp,store_mm_pp,g_coord_pp,numVar,&
                            sProp,diag_precon_pp,time_step)
  !/****f* parafeml/finddiagprecon
  !*  NAME
  !*    SUBROUTINE: finddiagprecon
  !*
  !*  SYNOPSIS
  !*    Usage: finddiagprecon(store_km_pp_OF_,store_mm_pp_OF_,      &
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
  !*    time_step                        - time step
  !*    sProp       (e,v,rho)            - Rheology Properties
  !*
  !*  OUTPUT
  !*   	diag_precon_pp (neq_pp)          - Diagonal Preconditioner
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
! 1. Declare variables
!------------------------------------------------------------------------------

  INTEGER,PARAMETER         :: nodof=3,ndim=3,nst=6,nod=8
  REAL(iwp),PARAMETER       :: zero=0.0_iwp
  
  INTEGER                   :: i,j,k,iel,ndof,nip

  REAL(iwp),INTENT(INOUT)   :: store_km_pp(ndim*nod,ndim*nod,nels_pp)
  REAL(iwp),INTENT(INOUT)   :: store_mm_pp(ndim*nod,ndim*nod,nels_pp)
  REAL(iwp),INTENT(INOUT)   :: g_coord_pp(nod,ndim,nels_pp)
  REAL(iwp),INTENT(IN)      :: numVar(4),sProp(3),time_step
  
  REAL(iwp),INTENT(INOUT)   :: diag_precon_pp(neq_pp)
  
  REAL(iwp)                 :: alpha1,beta1,theta
  REAL(iwp)                 :: dtim,c3,c4
  REAL(iwp)                 :: det,e,rho,v,volume

  CHARACTER(LEN=15)         :: element

  LOGICAL                   :: consistent=.TRUE.

!----------------------------------------------------------------------
! 2. Declare dynamic arrays
!----------------------------------------------------------------------

  REAL(iwp),ALLOCATABLE :: diag_precon_tmp(:,:)
  REAL(iwp),ALLOCATABLE :: points(:,:),dee(:,:),weights(:)
  REAL(iwp),ALLOCATABLE :: jac(:,:),der(:,:),deriv(:,:),bee(:,:)
  REAL(iwp),ALLOCATABLE :: fun(:),emm(:,:),ecm(:,:)

  INTEGER,ALLOCATABLE   :: node(:),localcount(:),readcount(:)

  IF(numpe .EQ. 1)PRINT*,"Finding Diagonal Preconditioner: "
  
  ! Set Parameters
  alpha1  =  numVar(1)
  beta1   =  numVar(2)
  theta   =  numVar(3)
  dtim    =  time_step !numVar(4)
  c3      =  alpha1+1._iwp/(theta*dtim)
  c4      =  beta1+theta*dtim
  ndof    =  ntot

!----------------------------------------------------------------------
! 6. Element Stiffness and Mass Integration
!----------------------------------------------------------------------

  ! Variables required for writing and partioning
  nip=8; element="hexahedron"

  ALLOCATE(points(nip,ndim),dee(nst,nst),jac(ndim,ndim))
  ALLOCATE(der(ndim,nod),deriv(ndim,nod),bee(nst,ntot))
  ALLOCATE(weights(nip),ecm(ntot,ntot),emm(ntot,ntot),fun(nod))

  ! Youngs Modulus, Poissons Ratio and Density
  e=sProp(1);v=sProp(2);rho=sProp(3);

  ! [D] : Stress-Strain Matrix
  dee  =  zero 
  CALL deemat(dee,e,v)
  CALL sample(element,points,weights)
  
  ! Clean [K] and [M] 
  store_km_pp  =  zero
  store_mm_pp  =  zero
  
  elements_1: DO iel=1,nels_pp
    volume   =  zero
    emm      =  zero
    ecm      =  zero
    
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
  END DO elements_1  
    
  DEALLOCATE(points,dee,jac,der,deriv,bee)
  DEALLOCATE(weights,ecm,emm,fun)
  
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

  SUBROUTINE checkforce(force,sense,node,solidPatchIDSize,nn)

  !/****f* parafemutils/checkforce
  !*  NAME
  !*    SUBROUTINE: checkforce
  !*
  !*  SYNOPSIS
  !*    Usage: checkforce_(force,sense,node_ensi,&solidPatchIDSize,&
  !*		                   paraFemSolidDisp,&numPoints);                                 
  !*
  !*  FUNCTION
  !*    Prints the forces to file in ENSI GOLD format  
  !*
  !*  INPUTS
  !*    force  (loaded_nodes*ndim)  - Mesh coordinates 
  !*    sense  (loaded_nodes*ndim)  - Vector to define x,y,z 
  !*    node   (loaded_nodes*ndim)  - Node Numbers of loaded nodes
  !*    solidPatchIDSize            - Size of boundary Mesh	
  !*    nn                          - # of Nodes  			
  !*
  !*  OUTPUT
  !*	  Writes to file "argv".ensi.NDLDS
  !*	  Contains information about loads on nodes
  !*			
  !*  AUTHOR
  !*    S. Hewitt
  !*    L. Margetts
  !******
  !*
  !*/

  !USE mpi_wrapper  !remove comment for serial compilation
  USE precision;  USE global_variables; USE mp_interface; 
  USE input;      USE output;           USE loading;
  USE timing;     USE maths;            USE gather_scatter
  USE steering;   USE new_library;
  
  IMPLICIT NONE
  
  !----------------------  Declarations --------------------------------!
  INTEGER,PARAMETER       :: nodof=3,ndim=3,nst=6,nod=8
  INTEGER                 :: solidPatchIDSize,nlen,nn,loaded_nodes,pos,i,j
  INTEGER,INTENT(IN)      :: sense(solidPatchIDSize*ndim),node(solidPatchIDSize*ndim)
  REAL(iwp),INTENT(IN)    :: force(solidPatchIDSize*ndim)
  REAL(iwp),PARAMETER     :: zero=0.0_iwp
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

  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
  !--------------------------------------------------------------------


  SUBROUTINE checkparafem(g_coord,g_num,rest_ensi,nn,els)

  !/****f* parafemutils/checkparafem
  !*  NAME
  !*    SUBROUTINE: checkparafem
  !*
  !*  SYNOPSIS
  !*    Usage: checkparafem_(MeshData,g_num_OF,rest_ensi,	&
  !*			&numPoints,&numCells);
  !*        
  !*  FUNCTION
  !*    Writes the mesh to file in ensi gold format 
  !*
  !*  INPUTS
  !*    g_coord   (ndim,nn)    - Mesh coordinates 		
  !*    g_num     (nod,nels)   - Steering Matrix (x-y-z)				
  !*    rest_ensi (nodof+1,nr) - Restrained Nodes
  !*    nn                     - # of Nodes 
  !*    els=nels               - # of Elements/Cells
  !*			
  !*  OUTPUT
  !*    Writes to file "argv".ensi.case and oter. "argv".ensi.xxx files
  !*    Contains information about Mesh, Material, Geometry and Restraints
  !*
  !*    File    Description 
  !*
  !*    MATID   Material ID's
  !*    NDLDS   Loaded Nodes
  !*    NDBND   Restrained Nodes 
  !*    geo     Geometery File
  !*    case    Case file
  !*
  !*  WARNINGS
  !*    It currently only works in serial.
  !*    
  !*  TODO
  !*    * Clean and update
  !*    * Write to work in parallel
  !*			
  !*  AUTHOR
  !*    S. Hewitt
  !******
  !*
  !*/

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

  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
  !--------------------------------------------------------------------

  SUBROUTINE POPULATE_G_COORD_PP(g_coord,g_coord_pp,g_num_pp,npes,nn,nod,ndim)

  !/****f* parafemutils/populate_g_coord_pp
  !*  NAME
  !*    SUBROUTINE: populate_g_coord_pp
  !*
  !*  SYNOPSIS
  !*    Usage: CALL populate_g_coord_pp(g_coord,g_coord_pp,g_num_pp,	&
  !* 					     npes,nn,nod,ndim)
  !*  FUNCTION
  !*    Populates g_coord_pp based on g_coord. It takes the global
  !*    list of coordinates and reduces them to the local coordinate 
  !*    lists based on the steering matrix g_num_pp.
  !*    It is used when the decomposition is done using ParaFEMs
  !*    methods
  !*
  !*  INPUTS
  !*    g_coord  (ndim,nn)  - Mesh coordinates	
  !*    g_num_pp (nod,nels) - Steering Matrix	
  !*	
  !*    npes                - Number of processes
  !*    nn                  - Number of nodes 
  !*    nod                 - Number of nodes per element
  !*    ndim                - Number of dimensions per node		  
  !*			
  !*  OUTPUT
  !*	  g_coord_pp (nod,ndim,nels_pp) - Coordinate matrix
  !*
  !*  WARNINGS
  !*    This version has been outdated by populate_g_coord_pp2.
  !*    
  !*  AUTHOR
  !*    Sam Hewitt
  !*
  !*  CREATION DATE
  !*    28th September 2016
  !*
  !*  COPYRIGHT
  !*    (c) University of Manchester 2007-2016
  !******
  !* G_COORD_PP exists on all processors
  !*/
  

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

  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
  !--------------------------------------------------------------------

  SUBROUTINE POPULATE_G_COORD_PP2(g_coord,g_coord_pp,g_num_pp,nn,nod,ndim)

  !/****f* parafemutils/populate_g_coord_pp2
  !*  NAME
  !*    SUBROUTINE: populate_g_coord_pp2
  !*
  !*  SYNOPSIS
  !*    Usage: CALL populate_g_coord_pp2(g_coord,g_coord_pp,g_num_pp,	&
  !* 					     nod,ndim)
  !*  FUNCTION
  !*    Populates g_coord_pp based on g_coord. It takes the global
  !*    list of coordinates and reduces them to the local coordinate 
  !*    lists based on the steering matrix g_num_pp.
  !*    It is used when decomposition has been done by Foam-Extend
  !*
  !*  INPUTS
  !*    g_coord  (ndim,nn)  - Mesh coordinates		
  !*    g_num_pp (nod,nels) - Steering Matrix		
  !*    nn                  - Number of nodes 
  !*    nod                 - Number of nodes per element
  !*    ndim                - Number of dimensions per node		
  !*
  !*  OUTPUT
  !*	  g_coord_pp (nod,ndim,nels_pp)	- Coordinate matrix
  !*    
  !*  AUTHOR
  !*    Sam Hewitt
  !*
  !*  CREATION DATE
  !*    6th May 2017
  !*
  !*  COPYRIGHT
  !*    (c) University of Manchester 2007-2016
  !******
  !* This version uses a simple method to populate g_coord_pp
  !*
  !*/

  USE precision; USE global_variables;

  IMPLICIT NONE

  INTEGER, INTENT(IN)          :: nn,nod,ndim
  INTEGER, INTENT(IN)          :: g_num_pp(nod,nels_pp)
  INTEGER                      :: iel, i, j, k, l ! loop counters 
  REAL(iwp), INTENT(INOUT)     :: g_coord_pp(nod,ndim,nels_pp) 
  REAL(iwp), INTENT(IN)        :: g_coord(ndim,nn)    
  REAL(iwp)                    :: zero = 0.0_iwp
  
  ! Simple Method used after Foam-Extend Decomposition
  DO iel=1,nels_pp
    DO j=1,nod 
	    k = g_num_pp(j,iel)
	    g_coord_pp(j,:,iel)=g_coord(:,k)
    ENDDO
  ENDDO
  
  END SUBROUTINE

  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
  !--------------------------------------------------------------------

  SUBROUTINE POPULATE_G_NUM_PP(g_num,g_num_pp,npes,nod,nels)

  !/****f* parafemutils/populate_g_num_pp
  !*
  !*  NAME
  !*    SUBROUTINE: populate_g_num_pp
  !*
  !*  SYNOPSIS
  !*    Usage: CALL populate_g_num_pp(g_num,g_num_pp,npes,nod,nels)
  !*                                     
  !*  FUNCTION
  !*    This subroutine was used in early versions of the code to
  !*    distribute a global steering matrix to a local one
  !*        
  !*  INPUTS
  !*    g_num  (nod,nels)      - Global steering matrix		
  !*    g_num_pp (nod,nels_pp) - Local Steering Matrix		
  !*    npes                   - Number of Processors 
  !*    nod                    - Number of nodes per element
  !*    nels                   - Total Number elements	
  !*
  !*  WARNINGS
  !*    This subroutine is currently not used it was relevant in
  !*    a previous version.
  !*   
  !*  AUTHOR
  !*    Sam Hewitt
  !*
  !*  CREATION DATE
  !*    28th September 2016
  !*
  !*  COPYRIGHT
  !*    (c) University of Manchester 2007-2016
  !******
  !*
  !*/

  USE precision;    USE global_variables; 
  USE mp_interface; USE mpi_wrapper

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

  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
  !--------------------------------------------------------------------

  SUBROUTINE closefem()
   USE mp_interface
   IMPLICIT NONE
   CALL SHUTDOWN()	
  END SUBROUTINE closefem

  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
  !--------------------------------------------------------------------

  SUBROUTINE of2sg(element,vector,nod)
   
  !/****f* parafemutils/of2sg
  !*  NAME
  !*    SUBROUTINE: of2sg
  !*
  !*  SYNOPSIS
  !*    Usage:      of2sg(element,g_num_pp(:,iel)) 
  !*      
  !*  FUNCTION
  !*    It rearranges the steering matrix from Foam-Extend format
  !*    to Smith-Griffiths format
  !*
  !*  INPUTS
  !*    element - Element Type	
  !*    nod     - Number of nodes per element
  !*
  !*  OUTPUT
  !*    vector  - New element node orders
  !*
  !*  TODO
  !*    * Include tetrahedral elements
  !* 			
  !*  AUTHOR
  !*    S. Hewitt
  !******
  !*
  !*/

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
      ! Foam-Extend to Smith-Griffiths:
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
    
    CASE('tetrahedron')
      SELECT CASE(nod)
        CASE(4)
        temp(:)      = vector(:)
        vector(1)  = temp(1)
        vector(2)  = temp(3)
        vector(3)  = temp(2)
        vector(4)  = temp(4)
      END SELECT

  END SELECT

  END SUBROUTINE of2sg

  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
  !--------------------------------------------------------------------

  SUBROUTINE gloads(gravlo_pp,specWeight,nn,nodof,nod,ndim,nr,g_coord,g_num_pp,rest)

  !/****f* parafemutils/gloads
  !*  NAME
  !*    SUBROUTINE: gloads
  !*
  !*  SYNOPSIS
  !*    Usage: gloads_(gravlo_,&specWeight,&gPoints_,&nodof,&nod,&ndim,
  !*                   &numRestrNodes_, mPoints_,g_num_pp_OF_,rest_);     
  !*
  !*  FUNCTION
  !*    Supplies a gravitational force vector   
  !*
  !*  INPUTS
  !*    specWeight   - Specific Weight (density  * gravitaional loading)
  !*    nn           - Total number of nodes in mesh
  !*    nodof        - Number of degrees of freedom per node
  !*    nod          - Number of Nodes per Element
  !*    ndim         - Number of dimensions		
  !*    nr           - Number of restrained nodes	
  !*
  !*    g_num_pp (nod,nels_pp)  - Element Steering matrix
  !*    g_coord  (ndim,nn)      - Global Coordinate matrix
  !*    rest     (nr,nodof+1)   - Matrix of restrained nodes
  !*
  !*  OUTPUT
  !*    gravlo_pp   (neq_pp)     - Vector of loads
  !*		
  !*  WARNINGS
  !*    * Default direction is negative Y
  !*    * REST must be rearranged before function call	
  !*
  !*  TODO
  !*    * Remove need for global coordinate matrix
  !*
  !*  AUTHOR
  !*    S. Hewitt
  !* 
  !******
  !*
  !*/
  
  USE mpi_wrapper  !remove comment for serial compilation
  USE precision;   USE global_variables; USE mp_interface; 
  USE input;       USE output;           USE loading;
  USE timing;      USE maths;            USE gather_scatter;
  USE steering;    USE new_library; 
  
  IMPLICIT NONE

  INTEGER                 :: i,j,k,iel,ndof,nodof,nn,nip
  INTEGER                 :: node_end,node_start,nodes_pp

  INTEGER,INTENT(IN)      :: nr,nod,ndim,g_num_pp(nod,nels_pp)
  INTEGER,INTENT(IN)      :: rest(nr,nodof+1)

  REAL(iwp),PARAMETER     :: zero=0.0_iwp

  REAL(iwp)               :: lamba,det,g_coord_pp(nod,ndim,nels_pp)
  REAL(iwp)               :: g_coord(ndim,nn),pmul_pp(ntot,nels_pp)
  REAL(iwp)               :: load_pp(ndim*nn)

  REAL(iwp),INTENT(IN)    :: specWeight

  REAL(iwp),INTENT(INOUT) :: gravlo_pp(neq_pp)

  CHARACTER(LEN=15)       :: element

  REAL(iwp),ALLOCATABLE   :: points(:,:),weights(:)
  REAL(iwp),ALLOCATABLE   :: fun(:),der(:,:),jac(:,:)

  !----------  Apply loading to integration Points ---------
  element="hexahedron" ; ndof = nodof*nod; nip = 8

  IF(numpe==1)PRINT*,"Calculating gravity loads:"
  
  CALL POPULATE_G_COORD_PP2(g_coord,g_coord_pp,g_num_pp,nn,nod,ndim)

  ALLOCATE(points(nip,ndim))
  ALLOCATE(weights(nip),fun(nod))
  ALLOCATE(der(ndim,nod),jac(ndim,ndim))
  
  pmul_pp=zero

  CALL sample(element,points,weights)
  
  elements_1: DO iel=1,nels_pp  
    gauss_points_1: DO i=1,nip     
      CALL shape_der(der,points,i)
      jac=MATMUL(der,g_coord_pp(:,:,iel))
      det=determinant(jac)
      CALL shape_fun(fun,points,i)
      pmul_pp(2:ndof-1:3,iel)=pmul_pp(2:ndof-1:3,iel)+fun(:)*det*weights(i)
    END DO gauss_points_1

    ! Sign included here to donate negative Y
    pmul_pp(:,iel)=-pmul_pp(:,iel)*specWeight
  END DO elements_1

  gravlo_pp = zero 
  CALL scatter(gravlo_pp,pmul_pp)

  ! Write the loads out to file
  !IF(.false.)
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
  !ENDIF
 
  END SUBROUTINE gloads

  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
  !--------------------------------------------------------------------

  SUBROUTINE writeToFile(fileNum,fileName,field,dim1)
  
  !/****f* parafemutils/writeToFile
  !*  NAME
  !*    SUBROUTINE: writeToFile
  !*
  !*  SYNOPSIS
  !*    Usage: CALL writeToFile(10,"g_num_pp.txt",g_num_pp,nod*nels_pp)    
  !*
  !*  FUNCTION
  !*    This function writes a float field to file. Its main purpose
  !*     is debugging and is not currrently called in the code.  
  !*
  !*  INPUTS
  !*    filenum   - File Number, any positive Integer
  !*    filename  - Name of output file
  !*    field     - Field to loop through
  !*    dim1      - Dimensions of the field
  !*
  !*  OUTPUT
  !*    File of the name "filename" with the field written out		
  !*
  !*  AUTHOR
  !*    S. Hewitt
  !* 
  !******
  !*
  !*/

  USE mpi_wrapper  !remove comment for serial compilation#
  USE precision; USE global_variables; USE mp_interface;
  USE input;     USE output;           USE loading; 
  USE timing;    USE maths;            USE gather_scatter;
  USE steering;  USE new_library; 
  
  IMPLICIT NONE
  
  INTEGER                     :: i
  INTEGER,INTENT(IN)          :: dim1,fileNum
  REAL(iwp),INTENT(INOUT)     :: field(dim1)
  CHARACTER(LEN=*),INTENT(IN) :: fileName  

  OPEN(fileNum,FILE=fileName,status = "replace")
    DO i=1,dim1
      WRITE(fileNum,'(E12.5)') field(i)
    ENDDO
  CLOSE(fileNum)

  END SUBROUTINE writeToFile

  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
  !--------------------------------------------------------------------

  SUBROUTINE writeToFileInt(fileNum,fileName,field,dim1)

  !/****f* parafemutils/writeToFileInt
  !*  NAME
  !*    SUBROUTINE: writeToFileInt
  !*
  !*  SYNOPSIS
  !*    Usage: CALL writeToFileInt(10,"g_num_pp.txt",g_num_pp,nod*nels_pp)    
  !*
  !*  FUNCTION
  !*    This function writes an Integer field to file. Its main purpose
  !*     is debugging and is not currrently called in the code.  
  !*
  !*  INPUTS
  !*    filenum   - File Number, any positive Integer
  !*    filename  - Name of output file
  !*    field     - Field to loop through
  !*    dim1      - Dimensions of the field
  !*
  !*  OUTPUT
  !*    File of the name "filename" with the field written out		
  !*
  !*  AUTHOR
  !*    S. Hewitt
  !* 
  !******
  !*
  !*/
  USE mpi_wrapper  !remove comment for serial compilation#
  USE precision; USE global_variables; USE mp_interface; 
  USE input;     USE output;           USE loading; 
  USE timing;    USE maths;            USE gather_scatter;
  USE steering;  USE new_library; 
    
  IMPLICIT NONE
    
  INTEGER                     :: i
  INTEGER,INTENT(IN)          :: dim1,fileNum
  INTEGER,INTENT(INOUT)       :: field(dim1)
  CHARACTER(LEN=*),INTENT(IN) :: fileName

  OPEN(fileNum,FILE=fileName,status = "replace")
    DO i=1,dim1
      WRITE(fileNum,*) field(i)
    ENDDO
  CLOSE(fileNum)

  END SUBROUTINE writeToFileInt

  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
  !--------------------------------------------------------------------

  SUBROUTINE WRITE_LARGESTRAIN(job_name,nn,nr,loaded_nodes,timest,nr_timest,inewton,nr_iters)
  !/****f* parafemnl/write_largestrain
  !*  NAME
  !*    SUBROUTINE: write_largestrain
  !*
  !*  SYNOPSIS
  !*    Usage:      CALL write_largestrain(timest)
  !*
  !*  FUNCTION
  !*    Master processor writes out brief details about the problem and 
  !*    some performance data
  !*
  !*  INPUTS
  !*
  !*    The following dynamic real array has the INTENT(IN) attribute:
  !*
  !*    timest(:)              : Holds timing information
  !*
  !*  AUTHOR
  !*    Sam Hewitt
  !*
  !*  CREATION DATE
  !*    21.03.2018
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/  

  USE mpi_wrapper;    USE precision;  USE global_variables; 
  USE mp_interface;   USE input;      USE output; 
  USE loading;        USE timing;     USE maths; 
  USE gather_scatter; USE steering;   USE new_library;
  USE large_strain;
  
  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)  :: job_name
  INTEGER, INTENT(IN)            :: nn,nr,inewton
  INTEGER, INTENT(IN)            :: loaded_nodes,nr_iters(inewton,1)
  REAL(iwp), INTENT(IN)          :: timest(17)
  REAL(iwp), INTENT(IN)          :: nr_timest(10,20)

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------
  
  CHARACTER(LEN=50)              :: fname
  INTEGER                        :: i          ! loop counter
 
  IF(numpe==1) THEN

    fname       = job_name(1:INDEX(job_name, " ")-1) // ".run"
    OPEN(11,FILE=fname,STATUS='REPLACE',ACTION='WRITE')     

!------------------------------------------------------------------------------
! 2. Write basic details about the problem
!------------------------------------------------------------------------------

    WRITE(11,'(/A)')   "BASIC JOB DATA                                  "     
 
    WRITE(11,'(A,I12)')    "Number of processors used                   ",npes 
    !WRITE(11,'(A,I12)')    "Number of nodes in the mesh                 ",nn
    WRITE(11,'(A,I12)')    "Number of nodes that were restrained        ",nr
    WRITE(11,'(A,I12)')    "Number of equations solved                  ",neq
    IF(loaded_nodes > 0) THEN
      WRITE(11,'(A,I12)')    "Number of loaded nodes                      ",   &
                              loaded_nodes 
    END IF
!------------------------------------------------------------------------------
! 3. Output timing data
!------------------------------------------------------------------------------

    WRITE(11,'(/3A)')   "PROGRAM SECTION EXECUTION TIMES                    ",&
                        "SECONDS  ", "%TOTAL    "
    WRITE(11,'(A,F12.6,F8.2)') "Load the Structure                          ",&
                           timest(3)-timest(2),                               &
                           ((timest(3)-timest(2))/(timest(17)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Set Starting conditions                     ",&
                           timest(4)-timest(3),                               &
                           ((timest(4)-timest(3))/(timest(17)-timest(1)))*100  

    WRITE(11,'(/A,I1,A)')   "NEWTON RAPHSON ITERATIONS (",inewton,")"
    WRITE(11,'(/A)')   "PCG ITERATIONS "
    DO i=1,inewton
        WRITE(11,'(I4)') nr_iters(i,1)
    ENDDO

    WRITE(11,'(A,F12.6,F8.2)') "Gather Displacement Increment               ",&
                           SUM(nr_timest(:,1)),                               &
                           ((SUM(nr_timest(:,1)))/(timest(17)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Build Stiffness and Mass Matricies          ",&
                           SUM(nr_timest(:,2)),                               &
                           ((SUM(nr_timest(:,2)))/(timest(17)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Set internal Forces                         ",&
                           SUM(nr_timest(:,3)),                               &
                           ((SUM(nr_timest(:,3)))/(timest(17)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Set Newmark and Residual                    ",&
                           SUM(nr_timest(:,4)),                               &
                           ((SUM(nr_timest(:,4)))/(timest(17)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Build the preconditioner                    ",&
                           SUM(nr_timest(:,5)),                               &
                           ((SUM(nr_timest(:,5)))/(timest(17)-timest(1)))*100                             
    WRITE(11,'(A,F12.6,F8.2)') "Solve equations                             ",&
                           SUM(nr_timest(:,6)),                               &
                          ((SUM(nr_timest(:,6)))/(timest(17)-timest(1)))*100 
    WRITE(11,'(A,F12.6,F8.2)') "Check convergence                           ",&
                           SUM(nr_timest(:,7)),                             &
                          ((SUM(nr_timest(:,7)))/(timest(17)-timest(1)))*100
    WRITE(11,'(/A,F12.6,F8.2)') "Total Time in N-R loop                      ",&
                           timest(15)-timest(5),                               &
                          ((timest(15)-timest(5))/(timest(17)-timest(1)))*100
    WRITE(11,'(A,F12.6,F8.2)') "Update Velocity and Acceleration            ",&
                            timest(16)-timest(15),                            &
                          ((timest(16)-timest(15))/(timest(17)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Gather Data to pass out                     ",&
                            timest(17)-timest(16),                            &
                          ((timest(17)-timest(16))/(timest(17)-timest(1)))*100  
    WRITE(11,'(A,F12.6,A/)')   "Total execution time                        ",&
                          timest(17)-timest(1),"  100.00"
    CLOSE(11)
    
  END IF
  
  RETURN

  END SUBROUTINE WRITE_LARGESTRAIN

  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
  !--------------------------------------------------------------------

  SUBROUTINE WRITE_SMALLSTRAIN(job_name,timest,iters)
  !/****f* parafemnl/write_smallstrain
  !*  NAME
  !*    SUBROUTINE: write_smallstrain
  !*
  !*  SYNOPSIS
  !*    Usage:      CALL write_smallstrain(timest)
  !*
  !*  FUNCTION
  !*    Master processor writes out brief details about the problem and 
  !*    some performance data
  !*
  !*  INPUTS
  !*
  !*    The following dynamic real array has the INTENT(IN) attribute:
  !*
  !*    timest(:)              : Holds timing information
  !*
  !*  AUTHOR
  !*    Sam Hewitt
  !*
  !*  CREATION DATE
  !*    21.03.2018
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/  

  USE mpi_wrapper;    USE precision;  USE global_variables; 
  USE mp_interface;   USE input;      USE output; 
  USE loading;        USE timing;     USE maths; 
  USE gather_scatter; USE steering;   USE new_library;
  USE large_strain;
  
  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)  :: job_name
  REAL(iwp), INTENT(IN)          :: timest(8)
  INTEGER,INTENT(IN)             :: iters

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------

  CHARACTER(LEN=50)              :: fname
  INTEGER                        :: i          ! loop counter


  IF(numpe==1) THEN

    fname       = job_name(1:INDEX(job_name, " ")-1) // ".run"
    OPEN(11,FILE=fname,STATUS='REPLACE',ACTION='WRITE')     

!------------------------------------------------------------------------------
! 2. Output timing data
!------------------------------------------------------------------------------

    WRITE(11,'(/3A)')   "PROGRAM SECTION EXECUTION TIMES                    ",&
                        "SECONDS  ", "%TOTAL    "
    WRITE(11,'(A,F12.6,F8.2)') "Load the Structure                          ",&
                           timest(3)-timest(2),                               &
                           ((timest(3)-timest(2))/(timest(8)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Set Starting conditions                     ",&
                           timest(4)-timest(3),                               &
                           ((timest(4)-timest(3))/(timest(8)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Set Displacement & Velocity                 ",&
                            timest(5)-timest(4),                            &
                          ((timest(5)-timest(4))/(timest(8)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Build RHS                                   ",&
                            timest(6)-timest(5),                            &
                          ((timest(6)-timest(5))/(timest(8)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)')"Solve Equations                             ",&
                            timest(7)-timest(6),                            &
                          ((timest(7)-timest(6))/(timest(8)-timest(1)))*100
    WRITE(11,'(A,F12.6,F8.2)') "Gather data to pass out                     ",&
                            timest(8)-timest(7),                            &
                          ((timest(8)-timest(7))/(timest(8)-timest(1)))*100  
    WRITE(11,'(A,F12.6,A/)')   "Total execution time                        ",&
                          timest(8)-timest(1),"  100.00"

    WRITE(11,'(A,12I3)')       "Number of Iterations                        ",&
                            iters
    CLOSE(11)
    
  END IF
  
  RETURN

  END SUBROUTINE WRITE_SMALLSTRAIN

  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
  !--------------------------------------------------------------------

  SUBROUTINE system_mem_usage(valueRSS,valueVM)

  !/****f* parafemutils/system_mem_usage
  !*  NAME
  !*    SUBROUTINE: system_mem_usage
  !*
  !*  SYNOPSIS
  !*    Usage: CALL system_mem_usage(ValueA,ValueB)    
  !*
  !*  FUNCTION
  !*    This function reads the /proc/procnum file of the
  !*    current process and outputs the RSS and VM memory
  !*    values, copying them into the inputs
  !*  
  !*
  !*  INPUTS
  !*    valueRSS   - RSS Memory (Resident Set Size)
  !*    valueVM    - VM Memory  (Virtual Memory)	
  !*
  !*  WARNINGS
  !*    If you are using a intel compiler, you must include 
  !*    the ifport module
  !*
  !*  NOTES
  !*    https://stackoverflow.com/questions/22028571/track-memory-usage-in-fortran-90
  !*
  !*  AUTHOR
  !*    S. Hewitt
  !* 
  !******
  !*
  !*/

  !use ifport !if on intel compiler

  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: valueRSS,valueVM
  CHARACTER(len=200)   :: filename=' '
  CHARACTER(len=80)    :: line
  CHARACTER(len=8)     :: pid_char=' '
  INTEGER              :: pid
  LOGICAL              :: ifxst

  valueRSS=-1; valueVM =-1 ! return negative number if not found

  !--- get process ID
  pid = getpid()
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

  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
  !--------------------------------------------------------------------


  INTEGER FUNCTION findnelspp()

  !/****f* parafemutils/findnelspp
  !*  NAME
  !*    findnelspp
  !*
  !*  SYNOPSIS
  !*    Unused in current version
  !*
  !*  FUNCTION 
  !*    This function acts an access fuction for the
  !*    Foam-Extend code, it returns the value of 
  !*    nels_pp, the number of elements per processor 
  !*    
  !*  AUTHOR
  !*    S. Hewitt
  !*
  !******
  !*
  !*/
    USE precision; USE global_variables;
    IMPLICIT NONE	
      findnelspp=nels_pp
      RETURN
  END FUNCTION

  INTEGER FUNCTION calcnelsppof(nels,npes)

  !/****f* parafemutils/calcnelsppof
  !*  NAME
  !*    calcnelsppof
  !*
  !*  SYNOPSIS
  !*    Unused in current version
  !*
  !*  FUNCTION 
  !*    This function is used to calculate the number 
  !*    of elements per processor, it is a very primitive 
  !*    form of mesh decomposition assigning elements to
  !*    processors in order
  !*
  !*  WARNINGS
  !*    This is a deprecated routine
  !*   
  !*  AUTHOR
  !*    S. Hewitt
  !*
  !******
  !*
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

  !/****f* parafemutils/setnelspp
  !*  NAME
  !*    setnelspp
  !*
  !*  SYNOPSIS
  !*    Usage: setnelspp_(&nels_pp_OF);
  !*
  !*  FUNCTION
  !*    This function is used to set an input variable
  !*    to the number of elements per processor, nels_pp
  !*
  !*  INPUTS
  !*    nCells  - Number of Cells
  !*   
  !*  AUTHOR
  !*    S. Hewitt
  !*
  !******
  !*
   USE precision; USE global_variables;
   IMPLICIT NONE
   INTEGER,INTENT(IN) :: nCells
   nels_pp=nCells
   setnelspp=1
   RETURN
  END FUNCTION

  SUBROUTINE setielstart()

  !/****f* parafemutils/setielstart
  !*  NAME
  !*    setielstart
  !*
  !*  SYNOPSIS
  !*    Usage: CALL setielstart()
  !*
  !*  FUNCTION
  !*    This function is used by initparafem to set
  !*    the starting element number on each processor 
  !*   
  !*  AUTHOR
  !*    S. Hewitt
  !*
  !******
  !*
  
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
  !/****f* parafemutils/findneqpp
  !*  NAME
  !*    setnelspp
  !*
  !*  SYNOPSIS
  !*    Usage: const int neq_pp_OF = findneqpp_();
  !*
  !*  FUNCTION
  !*    This function is used to return the neq_pp, the
  !*    number of equations per processor to Foam-Extend
  !*    
  !*  AUTHOR
  !*    S. Hewitt
  !*
  !******
  !*
   USE precision; USE global_variables;
   IMPLICIT NONE	
   findneqpp=neq_pp
   RETURN
  END FUNCTION
