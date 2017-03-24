/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "DyParaFEMSolid.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "tractionDisplacementFvPatchVectorField.H"
#include "velocityTractionDisplacementFvPatchVectorField.H"
#include "multiMaterial.H"
#include "twoDPointCorrector.H"
#include "fixedValuePointPatchFields.H"
#include "RBFInterpolation.H"
#include "TPSFunction.H"
#include "fortranRest.H"
#include <ctime>
#include <iostream>
#include <fstream>

using namespace rbf;


// * * * * * * * * * * * * ParaFEM Fortran Subroutines* * * * * * * * * * * * //
// - Will need editing when more complex elements introduced
const int nod = 8;       // Element type
const int ndim = 3;	 // Number of Dimensions
const int ntot=ndim*nod; // ntot

using namespace std;

extern"C"
{
    void initparafem_
    (
        double* g_coord,
        int* g_num,
        int* rest,
        int* nn, 
        const int* nels,
        const int* nr,
	double* SolidProperties,
        int* g_g_pp,
	double* stiff,
	double* mass,
	int* g_num_pp
    );

    int findneqpp_();
    int findnelspp_();

    int calcnelsppof_
    (
	int* numElements,
	int* numProcessors
    );

    void finddiagparafem_
    (
        double* stiff,
	double* mass,
	double* NumericalVariables,
        double* precon
    );

    void checkforce_
    (
        double* force,
        int* sense, 
        int* node,
        int* solidPatchIDSize,
        int* solidMeshSize
    );

    void runparafem_
    (
	double* NumericalVariables,
        double* f_ext, 
        int* f_node,
        int* solidPatchIDSize,
	double* paraFemSolidDisp,
	double* paraFemSolidVel,
	double* paraFemSolidAcel,
	double* time,
        int* nn,
        int* g_g_pp,
        int* g_num_pp,
	double* stiff,
	double* mass,
        double* precon,
	double* gravlo
    );

    void checkparafem_
    (
        double* mesh,
        int* g_num,
    	int* restraint,
        int* nn,
        int* nels
    );

    void gloads_
    (
	double* gravlo,
	double* gravity,
	int*	nn,
	int*	nodof,
	const int*	nod,
	const int*	ndim,
	int*	nr,
	double*	g_coord,
	int*	g_num_pp,
	int*	rest
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidSolvers
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(DyParaFEMSolid, 0);
addToRunTimeSelectionTable(solidSolver, DyParaFEMSolid, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DyParaFEMSolid::DyParaFEMSolid(const fvMesh& mesh)
:
    solidSolver(typeName, mesh),
    mPoints_(NULL),
    solidProps_(NULL),
    numSchemes_(NULL),
    g_num_OF_(NULL), 
    g_num_pp_OF_(NULL),
    g_g_pp_OF_(NULL),
    store_km_pp_OF_(NULL),
    store_mm_pp_OF_(NULL),
    diag_precon_pp_OF_(NULL),
    gravlo_(NULL),
    d_OF_(NULL),
    u_OF_(NULL),
    a_OF_(NULL),
    numRestrNodes_(0),
    rest_(NULL),
    rest_ensi_(NULL),
    numFixedForceNodes_(0),
    forceNodes_(NULL),
    fext_OF_(NULL),
    nodeensi_(NULL),
    sense_(NULL),
    E_(0),
    nu_(0),
    rhotmp_(0),
    D_
    (
        IOobject
        (
            "D",
            runTime().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    U_
    (
        IOobject
        (
//             "U",
            "ddt(" + D_.name() + ")",
            runTime().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
//         fvc::ddt(D_)
        mesh,
        dimensionedVector("0", dimVelocity, vector::zero)
    ),
    pMesh_(mesh),
    pointD_
    (
        IOobject
        (
            "pointD",
            runTime().timeName(),
            mesh,
            // IOobject::READ_IF_PRESENT,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh_
	// dimensionedVector("0", dimLength, vector::zero)
    ),
    pointU_
    (
        IOobject
        (
            "pointU",
            runTime().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            //IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh_,
	dimensionedVector("0", dimVelocity, vector::zero)
    ),
    pointA_
    (
        IOobject
        (
            "pointA",
            runTime().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            //IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh_,
	dimensionedVector("0", dimAcceleration, vector::zero)
    ),
    sigma_
    (
        IOobject
        (
            "sigma",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    ),
    epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    rheology_(sigma_, D_),
    volToPoint_(mesh),
    pointToVol_(pMesh_,mesh),
    rho_("rho", rheology_.rho()),
    mu_(rheology_.mu()),
    lambda_(rheology_.lambda()),
    interface_(NULL)
{
    pointD_.oldTime();

    if (rheology_.law().type() == multiMaterial::typeName)
    {
        interface_.set(new TLMaterialInterface(D_, pointD_));
    }

    // - ParaFEM construction 
    E_ = rheology_.law().E()()[0];
    nu_ = rheology_.law().nu()()[0];
    rhotmp_ = rheology_.law().rho()()[0];

    // - Turn on/off the solid mesh Motion
    Switch moveMesh(solidProperties().lookup("moveMesh"));

    label numPoints = mesh.points().size();
    label numCells = mesh.nCells();

    mPoints_ = new double[numPoints*ndim]; 
    
    label globalIndex = 0;
    forAll(mesh.points(), pointI)
    {
	mPoints_[globalIndex++] = mesh.points()[pointI].x();
	mPoints_[globalIndex++] = mesh.points()[pointI].y();		
	mPoints_[globalIndex++] = mesh.points()[pointI].z();
    }

    g_num_OF_ = new int [numCells*nod];
   
    globalIndex = 0;
    const labelListList& cellPoints = mesh.cellPoints();
    forAll(cellPoints, cellI)
    {
        const labelList& curCellPoints = cellPoints[cellI];
	
	if (curCellPoints.size() != nod)
	{
            Info << "Not a hex cell!" << endl;
	}

	for(label i=0;i<nod;i++)
	{
	    g_num_OF_[globalIndex++] = curCellPoints[i]+1; // +1 fortran
	}	 
    }

// ------------- RESTRAINED ARRAYS ------------- //
// Create RESTRAINED array (0=Restrained 1=Unrestrained)
// (ENSI GOLD format is reversed)


    label nRows = mesh.points().size();
    rest_ensi_ = new int [nRows*4];
    fortranRest rest_ensi_obj(nRows);

    for (label i=0; i<numPoints; i++)
    {
	rest_ensi_obj.addNode(i,0,0,0);
    }


    forAll(pointD_.boundaryField(), patchI)
    {
 	//- Restrain the fixedValue pointD fields
        if 
        (
            isA<fixedValuePointPatchVectorField>
            (
                pointD_.boundaryField()[patchI]
            ) 
        )
	{
	    const labelList& mp = mesh.boundaryMesh()[patchI].meshPoints();

	    forAll(mp, pI)
	    {
                label i = mp[pI];
		rest_ensi_obj.editNode(i,1,1,1);		
	    }
	}
	
	if 
	(
            (mesh.boundary()[patchI].type() == "empty") || 
	    (mesh.boundary()[patchI].type() == "symmetryPlane")
	)

	{
            //numRestrNodes_ += mesh.boundaryMesh()[patchI].meshPoints().size();
	    const labelList& mp = mesh.boundaryMesh()[patchI].meshPoints();

	    // empty basically implies symmetry in Z
	    int x=0;
	    int y=0;
	    int z=0;
	    
            // symmetry boundaries must have these names
	    if (mesh.boundaryMesh()[patchI].name() == "symmetry-x")
	    {
	        x=1;
		y=2;
   		z=2;
	    }
	    else if (mesh.boundaryMesh()[patchI].name() == "symmetry-y")
	    {
	        x=2;
		y=1;
   		z=2;
	    }
	    else if (mesh.boundaryMesh()[patchI].type() == "empty")
	    {
	        x=2;
		y=2;
   		z=1;
	    }

	    forAll(mp, pI)
	    {
		// WARNING CURRENTLY SET UP FOR FIXED Z
		label i = mp[pI];
		rest_ensi_obj.editNode(i,x,y,z);
	    }
	}
    }

    rest_ensi_obj.printENSI(rest_ensi_);
    numRestrNodes_ = rest_ensi_obj.getNumRestNodes();
    rest_ = new int [numRestrNodes_*4];
    rest_ensi_obj.printFromENSI(rest_, numRestrNodes_);

// ------------- DECLARE ARRAYS ------------- //
    int size_ = Pstream::nProcs();
    calcnelsppof_(&numCells,&size_);
    const int nels_pp_OF = findnelspp_();
    
    g_num_pp_OF_ = new int [nod*nels_pp_OF];
    g_g_pp_OF_ = new int [ntot*nels_pp_OF];
    store_km_pp_OF_ = new double [ntot*ntot*nels_pp_OF];
    store_mm_pp_OF_ = new double [ntot*ntot*nels_pp_OF];

// ------------- SOLID SOLUTION PROPERTIES ------------- //
    double alpha1 (readScalar(solidProperties().lookup("alpha1")));
    double beta1 (readScalar(solidProperties().lookup("beta1")));
    double timestep (readScalar(solidProperties().lookup("timeStep")));
    double theta (readScalar(solidProperties().lookup("theta")));
    
    numSchemes_ = new double[4];
    numSchemes_[0] = alpha1;
    numSchemes_[1] = beta1;
    numSchemes_[2] = theta;
    numSchemes_[3] = timestep;

// ------------- SOLID RHEOLOGY PROPERTIES ------------- //
    solidProps_ = new double[3];
    solidProps_[0] = E_;
    solidProps_[1] = nu_;
    solidProps_[2] = rhotmp_;

//    if((Pstream::master()==true) or (Pstream::parRun()==false))
//    {
//        checkparafem_
//        (
//            mPoints_, 
//            g_num_OF_, 
//            rest_ensi_, 
//            &numPoints, 
//            &numCells
//        );
//    }
//    std::clock_t start;
//    double initFEM = 0.0;
//    start = std::clock();

    label tmp = Pstream::myProcNo();
    reduce(tmp,sumOp<label>());
    initparafem_
    (
        mPoints_,
        g_num_OF_,
        rest_,
        &numPoints,
        &numCells,
        &numRestrNodes_,
        solidProps_,
        g_g_pp_OF_,
        store_km_pp_OF_,
	store_mm_pp_OF_,
	g_num_pp_OF_
    );
    reduce(tmp,sumOp<label>());

//    initFEM =  (std::clock()-start)/(double) CLOCKS_PER_SEC;
//    Info << "initFEM: " << initFEM << endl;

    const int neq_pp_OF = findneqpp_();
    diag_precon_pp_OF_ = new double [neq_pp_OF];
    
    double gravity(readScalar(solidProperties().lookup("gravity")));
    gravlo_ = new double [neq_pp_OF];


    if(gravity > 1e-6)
    { 
	    Info << "Gravity Loading, gravity: " << gravity << " m/s^2" << endl; 

	    // Specific weight lambda = rho * g (gloads: direction of negative y implied)
	    double specWeight=gravity*rhotmp_;;
	    int nodof=3;	    

	    gloads_ 
	    (
	    gravlo_,
	    &specWeight,
	    &numPoints,
	    &nodof,
	    &nod,
	    &ndim,
	    &numRestrNodes_,
	    mPoints_,
	    g_num_pp_OF_,
	    rest_
	    );
    }
    else
    {
    	for(int i=0;i<neq_pp_OF;i++)
    	{
	    gravlo_[i]=0.0;
 	}
    }
    
    finddiagparafem_(store_km_pp_OF_,store_mm_pp_OF_,numSchemes_,diag_precon_pp_OF_);

// ------------- CALCULATING FORCE ARRAYS ------------- //
    // Identify traction specified points
    forAll(D_.boundaryField(), patchI)
    {
        if
        (
            isA<tractionDisplacementFvPatchVectorField>
            (
                D_.boundaryField()[patchI]
            )
        )
        {
            numFixedForceNodes_ += 
                mesh.boundaryMesh()[patchI].meshPoints().size();
        }
    }

    forceNodes_ = new int [numFixedForceNodes_*ndim];
    fext_OF_ = new double [numFixedForceNodes_*ndim];
    nodeensi_ = new int [numFixedForceNodes_*ndim];
    sense_ = new int [numFixedForceNodes_*ndim];

    label gi = 0;
    forAll(D_.boundaryField(), patchI)
    {
        if
        (
            isA<tractionDisplacementFvPatchVectorField>
            (
                D_.boundaryField()[patchI]
            )
        )
        {
	    const labelList& mp = 
                mesh.boundaryMesh()[patchI].meshPoints();

	    jj_=0;
	    forAll(mp, pI)
            {	      
                forceNodes_[gi++] = mp[pI]+1;
                nodeensi_[jj_] = mp[pI];
                nodeensi_[jj_+1] = mp[pI];
                nodeensi_[jj_+2] = mp[pI];
        	sense_[jj_]=1;
	        sense_[jj_+1]=2;
	        sense_[jj_+2]=3;
 		jj_=jj_+3;
	    }
        }
    }

    // Set force to zero
    for(int i=0; i<numFixedForceNodes_*ndim; i++)
    {
 	fext_OF_[i]=0;
    }	

    d_OF_ = new double [numPoints*ndim];
    u_OF_ = new double [numPoints*ndim];
    a_OF_ = new double [numPoints*ndim];

    // Set dispField to zero
    for(int i=0; i<numPoints*ndim; i++)
    {
        d_OF_[i]=0;
        u_OF_[i]=0;
        a_OF_[i]=0;
    }

}

DyParaFEMSolid::~DyParaFEMSolid()
{
    delete[] store_km_pp_OF_;
    delete[] store_mm_pp_OF_;
    delete[] diag_precon_pp_OF_;
    delete[] d_OF_;
    delete[] u_OF_;
    delete[] a_OF_;
    delete[] rest_;
    delete[] solidProps_;
    delete[] mPoints_;
    delete[] g_num_OF_;
    delete[] g_g_pp_OF_;
    delete[] rest_ensi_;
    delete[] numSchemes_;
    delete[] g_num_pp_OF_;
    delete[] forceNodes_;
    delete[] fext_OF_;
    delete[] nodeensi_;
    delete[] sense_;
    delete[] gravlo_;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vector DyParaFEMSolid::pointU(label pointID) const
{
    pointVectorField pointU
    (
        IOobject
        (
            "pointU",
            runTime().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pMesh_,
        dimensionedVector("0", dimVelocity, vector::zero)
    );

    volToPoint_.interpolate(U_, pointU);

    return pointU.internalField()[pointID];
}

//- Patch point displacement
tmp<vectorField> DyParaFEMSolid::patchPointDisplacementIncrement
(
    const label patchID
) const
{
    tmp<vectorField> tPointDisplacement
    (
        new vectorField
        (
            mesh().boundaryMesh()[patchID].localPoints().size(), 
            vector::zero
        )
    );

    tPointDisplacement() = 
        vectorField
        (
            pointD_.internalField() - pointD_.oldTime().internalField(), 
            mesh().boundaryMesh()[patchID].meshPoints()
        );

    return tPointDisplacement;
}

//- Face zone point displacement
tmp<vectorField> DyParaFEMSolid::faceZonePointDisplacementIncrement
(
    const label zoneID
) const
{

    tmp<vectorField> tPointDisplacement
    (
        new vectorField
        (
            mesh().faceZones()[zoneID]().localPoints().size(), 
            vector::zero
        )
    );

    vectorField& pointDisplacement = tPointDisplacement();

    const vectorField& pointDI = pointD_.internalField();
    const vectorField& oldPointDI = pointD_.oldTime().internalField();
    
    label globalZoneIndex = findIndex(globalFaceZones(), zoneID);
 
   
    // global zone index is set in the decomposeParDict
    if (globalZoneIndex != -1)
//    if (false)
    {
	Info << "WARNING FIELDS ARE STORED ON EVERY PROCESSOR" << endl;
        // global face zone

        const labelList& curPointMap =
            globalToLocalFaceZonePointMap()[globalZoneIndex];
           
        const labelList& zoneMeshPoints =
            mesh().faceZones()[zoneID]().meshPoints();

        vectorField zonePointsDisplGlobal
        (
            zoneMeshPoints.size(),
            vector::zero
        );

        //- Inter-proc points are shared by multiple procs
        //  pointNumProc is the number of procs which a point lies on
        scalarField pointNumProcs(zoneMeshPoints.size(), 0);
        
        forAll(zonePointsDisplGlobal, globalPointI)
        {
            label localPoint = curPointMap[globalPointI];
	    
            if(zoneMeshPoints[localPoint] < mesh().nPoints())
            {
                label procPoint = zoneMeshPoints[localPoint];
                
                zonePointsDisplGlobal[globalPointI] = 
                    pointDI[procPoint] - oldPointDI[procPoint];

                pointNumProcs[globalPointI] = 1;
            }
        }

// The following block was commented out and im not sure why 

        if (Pstream::parRun())
        {
            reduce(zonePointsDisplGlobal, sumOp<vectorField>());
            reduce(pointNumProcs, sumOp<scalarField>());

            //- now average the displacement between all procs
            zonePointsDisplGlobal /= pointNumProcs;
        }

//---------------------
        forAll(pointDisplacement, globalPointI)
        {
            label localPoint = curPointMap[globalPointI];
	    
            pointDisplacement[localPoint] = 
                zonePointsDisplGlobal[globalPointI];
        }
    }
    else
    {

        tPointDisplacement() = 
            vectorField
            (
                pointDI - oldPointDI, 
                mesh().faceZones()[zoneID]().meshPoints()
            );

    }
    return tPointDisplacement;
}

//- Patch point displacement
tmp<vectorField> DyParaFEMSolid::patchPointDisplacement
(
    const label patchID
) const
{
    Info << "patchPointDisplacement" << endl;
    tmp<vectorField> tPointDisplacement
    (
        new vectorField
        (
            mesh().boundaryMesh()[patchID].localPoints().size(), 
            vector::zero
        )
    );

    tPointDisplacement() = 
        vectorField
        (
            pointD_.oldTime().internalField(), 
            mesh().boundaryMesh()[patchID].meshPoints()
        );

    return tPointDisplacement;
}


//- Face zone point displacement
tmp<vectorField> DyParaFEMSolid::faceZonePointDisplacement
(
    const label zoneID
) const
{
    tmp<vectorField> tPointDisplacement
    (
        new vectorField
        (
            mesh().faceZones()[zoneID]().localPoints().size(), 
            vector::zero
        )
    );
    vectorField& pointDisplacement = tPointDisplacement();

    const vectorField& oldPointDI = pointD_.oldTime().internalField();
//    const vectorField& oldPointDI = pointD_.prevIter();

    label globalZoneIndex = findIndex(globalFaceZones(), zoneID);


    if (globalZoneIndex != -1)
    {
        // global face zone

        const labelList& curPointMap =
            globalToLocalFaceZonePointMap()[globalZoneIndex];

        const labelList& zoneMeshPoints =
            mesh().faceZones()[zoneID]().meshPoints();

        vectorField zonePointsDisplGlobal
        (
            zoneMeshPoints.size(),
            vector::zero
        );
        
        //- Inter-proc points are shared by multiple procs
        //  pointNumProc is the number of procs which a point lies on
        scalarField pointNumProcs(zoneMeshPoints.size(), 0);

        forAll(zonePointsDisplGlobal, globalPointI)
        {
            label localPoint = curPointMap[globalPointI];
	    
            if(zoneMeshPoints[localPoint] < mesh().nPoints())
            {
                label procPoint = zoneMeshPoints[localPoint];
                
                zonePointsDisplGlobal[globalPointI] = 
                    oldPointDI[procPoint];

                pointNumProcs[globalPointI] = 1;
            }
        }



//        if (Pstream::parRun())
//        {
//            reduce(zonePointsDisplGlobal, sumOp<vectorField>());
//            reduce(pointNumProcs, sumOp<scalarField>());

            //- now average the displacement between all procs
//            zonePointsDisplGlobal /= pointNumProcs;
//        }

        forAll(pointDisplacement, globalPointI)
        {
            label localPoint = curPointMap[globalPointI];

            pointDisplacement[localPoint] = 
                zonePointsDisplGlobal[globalPointI];
        }
    }
    else
    {
        tPointDisplacement() = 
            vectorField
            (
                oldPointDI, 
                mesh().faceZones()[zoneID]().meshPoints()
            );
    }

    return tPointDisplacement;
}

//- Patch face acceleration
tmp<Foam::vectorField> DyParaFEMSolid::patchFaceAcceleration
(
    const label patchID
) const
{
    tmp<vectorField> tAcceleration
    (
        new vectorField
        (
            mesh().boundary()[patchID].size(),
            vector::zero
        )
    );

    volVectorField a = fvc::ddt(U_);

    tAcceleration() = a.boundaryField()[patchID];

    return tAcceleration;
}

//- Patch face acceleration
tmp<vectorField> DyParaFEMSolid::faceZoneAcceleration
(
    const label zoneID,
    const label patchID
) const
{
    tmp<vectorField> tAcceleration
    (
        new vectorField
        (
            mesh().faceZones()[zoneID]().size(),
            vector::zero
        )
    );
    vectorField& acceleration = tAcceleration();

    volVectorField a = fvc::ddt(U_);

    vectorField patchAcceleration = a.boundaryField()[patchID];

    const label patchStart = 
        mesh().boundaryMesh()[patchID].start();

    forAll(patchAcceleration, i)
    {
        acceleration[mesh().faceZones()[zoneID].whichFace(patchStart + i)] = 
            patchAcceleration[i];
    }

    // Parallel data exchange: collect pressure field on all processors
    reduce(acceleration, sumOp<vectorField>());

    return tAcceleration;
}


//- Patch face acceleration
tmp<vectorField> DyParaFEMSolid::faceZoneVelocity
(
    const label zoneID,
    const label patchID
) const
{
    tmp<vectorField> tVelocity
    (
        new vectorField
        (
            mesh().faceZones()[zoneID]().size(),
            vector::zero
        )
    );
    vectorField& velocity = tVelocity();

    vectorField patchVelocity = U_.boundaryField()[patchID];

    const label patchStart = 
        mesh().boundaryMesh()[patchID].start();

    forAll(patchVelocity, i)
    {
        velocity[mesh().faceZones()[zoneID].whichFace(patchStart + i)] = 
            patchVelocity[i];
    }

    // Parallel data exchange: collect pressure field on all processors
    reduce(velocity, sumOp<vectorField>());

    return tVelocity;
}


//- Face zone point displacement
tmp<tensorField> DyParaFEMSolid::faceZoneSurfaceGradientOfVelocity
(
    const label zoneID,
    const label patchID
) const
{
    tmp<tensorField> tVelocityGradient
    (
        new tensorField
        (
            mesh().faceZones()[zoneID]().size(),
            tensor::zero
        )
    );

    return tVelocityGradient;
}


tmp<vectorField> 
DyParaFEMSolid::currentFaceZonePoints(const label zoneID) const
{
    vectorField pointDisplacement
    (
        mesh().faceZones()[zoneID]().localPoints().size(), 
        vector::zero
    );

    const vectorField& pointDI = pointD_.internalField();

    label globalZoneIndex = findIndex(globalFaceZones(), zoneID);

    if (globalZoneIndex != -1)
    {
        // global face zone
        const labelList& curPointMap =
            globalToLocalFaceZonePointMap()[globalZoneIndex];

        const labelList& zoneMeshPoints =
            mesh().faceZones()[zoneID]().meshPoints();

        vectorField zonePointsDisplGlobal
        (
            zoneMeshPoints.size(),
            vector::zero
        );
        
        //- Inter-proc points are shared by multiple procs
        //  pointNumProc is the number of procs which a point lies on
        scalarField pointNumProcs(zoneMeshPoints.size(), 0);

        forAll(zonePointsDisplGlobal, globalPointI)
        {
            label localPoint = curPointMap[globalPointI];
	    
            if(zoneMeshPoints[localPoint] < mesh().nPoints())
            {
                label procPoint = zoneMeshPoints[localPoint];
                
                zonePointsDisplGlobal[globalPointI] = 
                    pointDI[procPoint];

                pointNumProcs[globalPointI] = 1;
            }
        }

        if (Pstream::parRun())
        {
            reduce(zonePointsDisplGlobal, sumOp<vectorField>());
            reduce(pointNumProcs, sumOp<scalarField>());

            //- now average the displacement between all procs
            zonePointsDisplGlobal /= pointNumProcs;
        }

        forAll(pointDisplacement, globalPointI)
        {
            label localPoint = curPointMap[globalPointI];
	    
            pointDisplacement[localPoint] = 
                zonePointsDisplGlobal[globalPointI];
        }
    }
    else
    {
        pointDisplacement = 
            vectorField
            (
                pointDI, 
                mesh().faceZones()[zoneID]().meshPoints()
            );
    }

    tmp<vectorField> tCurrentPoints
    (
        new vectorField
        (
            mesh().faceZones()[zoneID]().localPoints() 
          + pointDisplacement
        )
    );

    return tCurrentPoints;
}


//- Face zone point displacement
tmp<vectorField> DyParaFEMSolid::faceZoneNormal
(
    const label zoneID,
    const label patchID
) const
{
    tmp<vectorField> tNormals
    (
        new vectorField
        (
            mesh().faceZones()[zoneID]().size(),
            vector::zero
        )
    );
    vectorField& normals = tNormals();

    const faceList& localFaces =
        mesh().boundaryMesh()[patchID].localFaces();

    vectorField localPoints =
        mesh().boundaryMesh()[patchID].localPoints();
    localPoints += pointD_.boundaryField()[patchID].patchInternalField();

    PrimitivePatch<face, List, const pointField&> patch
    (
        localFaces,
        localPoints
    );

    vectorField patchNormals(patch.size(), vector::zero);

    forAll(patchNormals, faceI)
    {
        patchNormals[faceI] =
            localFaces[faceI].normal(localPoints);
    }

    label globalZoneIndex = findIndex(globalFaceZones(), zoneID);

    if (globalZoneIndex != -1)
    {
        // global face zone

        const label patchStart = 
            mesh().boundaryMesh()[patchID].start();

        forAll(patchNormals, i)
        {
            normals
            [
                mesh().faceZones()[zoneID].whichFace(patchStart + i)
            ] = 
                patchNormals[i];
        }

        // Parallel data exchange: collect field on all processors
        reduce(normals, sumOp<vectorField>());
    }
    else
    {
        normals = patchNormals;
    }

    return tNormals;
}

void DyParaFEMSolid::setTraction
(
    const label patchID,
    const vectorField& traction
)
{
    if
    (
        D_.boundaryField()[patchID].type()
     != tractionDisplacementFvPatchVectorField::typeName
    )
    {
        FatalErrorIn("void DyParaFEMSolid::setTraction(...)")
            << "Bounary condition on " << D_.name() 
                <<  " is " 
                << D_.boundaryField()[patchID].type() 
                << "for patch" << mesh().boundary()[patchID].name()
                << ", instead " 
                << tractionDisplacementFvPatchVectorField::typeName
                << abort(FatalError);
    }

    tractionDisplacementFvPatchVectorField& patchU =
        refCast<tractionDisplacementFvPatchVectorField>
        (
            D_.boundaryField()[patchID]
        );

    patchU.traction() = traction;
}

void DyParaFEMSolid::setPressure
(
    const label patchID,
    const scalarField& pressure
)
{
    if
    (
        D_.boundaryField()[patchID].type()
     != tractionDisplacementFvPatchVectorField::typeName
    )
    {
        FatalErrorIn("void DyParaFEMSolid::setTraction(...)")
            << "Bounary condition on " << D_.name() 
                <<  " is " 
                << D_.boundaryField()[patchID].type() 
                << "for patch" << mesh().boundary()[patchID].name()
                << ", instead " 
                << tractionDisplacementFvPatchVectorField::typeName
                << abort(FatalError);
    }

    tractionDisplacementFvPatchVectorField& patchU =
        refCast<tractionDisplacementFvPatchVectorField>
        (
            D_.boundaryField()[patchID]
        );

    patchU.pressure() = pressure;
}

void DyParaFEMSolid::setTraction
(
    const label patchID,
    const label zoneID,
    const vectorField& faceZoneTraction
)
{
    vectorField patchTraction(mesh().boundary()[patchID].size(), vector::zero);

    const label patchStart = 
        mesh().boundaryMesh()[patchID].start();

    forAll(patchTraction, i)
    {
        patchTraction[i] =
            faceZoneTraction
            [
                mesh().faceZones()[zoneID].whichFace(patchStart + i)
            ];
    }

    setTraction(patchID, patchTraction);
}



void DyParaFEMSolid::setPressure
(
    const label patchID,
    const label zoneID,
    const scalarField& faceZonePressure
)
{
    scalarField patchPressure(mesh().boundary()[patchID].size(), 0.0);

    const label patchStart = 
        mesh().boundaryMesh()[patchID].start();

    forAll(patchPressure, i)
    {
        patchPressure[i] =
            faceZonePressure
            [
                mesh().faceZones()[zoneID].whichFace(patchStart + i)
            ];
    }

    setPressure(patchID, patchPressure);
}


//- Set traction at specified patch
void DyParaFEMSolid::setVelocityAndTraction
(
    const label patchID,
    const vectorField& traction,
    const vectorField& velocity,
    const vectorField& normal
)
{
    if
    (
        D_.boundaryField()[patchID].type()
     != velocityTractionDisplacementFvPatchVectorField::typeName
    )
    {
        FatalErrorIn
        (
            "void DyParaFEMSolid::"
            "setVelocityAndTraction(...)"
        )
            << "Bounary condition on " << D_.name() 
                <<  " is " 
                << D_.boundaryField()[patchID].type() 
                << "for patch" << mesh().boundary()[patchID].name()
                << ", instead " 
                << velocityTractionDisplacementFvPatchVectorField::typeName
                << abort(FatalError);
    }

    velocityTractionDisplacementFvPatchVectorField& patchU =
        refCast<velocityTractionDisplacementFvPatchVectorField>
        (
            D_.boundaryField()[patchID]
        );

    patchU.traction() = traction;
    patchU.velocity() = velocity;
    patchU.normal() = normal;
}
   
//- Set traction at specified patch
void DyParaFEMSolid::setVelocityAndTraction
(
    const label patchID,
    const label zoneID,
    const vectorField& faceZoneTraction,
    const vectorField& faceZoneVelocity,
    const vectorField& faceZoneNormal
)
{
    vectorField patchTraction(mesh().boundary()[patchID].size(), vector::zero);
    vectorField patchVelocity(mesh().boundary()[patchID].size(), vector::zero);
    vectorField patchNormal(mesh().boundary()[patchID].size(), vector::zero);

    const label patchStart = 
        mesh().boundaryMesh()[patchID].start();

    forAll(patchTraction, i)
    {
        patchTraction[i] =
            faceZoneTraction
            [
                mesh().faceZones()[zoneID].whichFace(patchStart + i)
            ];
        patchVelocity[i] =
            faceZoneVelocity
            [
                mesh().faceZones()[zoneID].whichFace(patchStart + i)
            ];
        patchNormal[i] =
            faceZoneNormal
            [
                mesh().faceZones()[zoneID].whichFace(patchStart + i)
            ];
    }

    setVelocityAndTraction
    (
        patchID, 
        patchTraction, 
        patchVelocity, 
        patchNormal
    );
}

tmp<vectorField> DyParaFEMSolid::predictTraction
(
    const label patchID,
    const label zoneID
)
{
    // Predict traction on patch
    if
    (
        D_.boundaryField()[patchID].type()
     != tractionDisplacementFvPatchVectorField::typeName
    )
    {
        FatalErrorIn("void DyParaFEMSolid::setTraction(...)")
            << "Bounary condition on " << D_.name() 
                <<  " is " 
                << D_.boundaryField()[patchID].type() 
                << "for patch" << mesh().boundary()[patchID].name()
                << ", instead " 
                << tractionDisplacementFvPatchVectorField::typeName
                << abort(FatalError);
    }

    tractionDisplacementFvPatchVectorField& patchUo =
        refCast<tractionDisplacementFvPatchVectorField>
        (
            D_.oldTime().boundaryField()[patchID]
        );

    tractionDisplacementFvPatchVectorField& patchUoo =
        refCast<tractionDisplacementFvPatchVectorField>
        (
            D_.oldTime().oldTime().boundaryField()[patchID]
        );


    vectorField ptF = 2*patchUo.traction() - patchUoo.traction();

    tmp<vectorField> ttF
    (
        new vectorField(mesh().faceZones()[zoneID].size(), vector::zero)
    );
    vectorField& tF = ttF();

    const label patchStart =
        mesh().boundaryMesh()[patchID].start();

    forAll(ptF, i)
    {
        tF[mesh().faceZones()[zoneID].whichFace(patchStart + i)] = ptF[i];
    }

    // Parallel data exchange: collect pressure field on all processors
    reduce(tF, sumOp<vectorField>());

    return ttF;
}


tmp<scalarField> DyParaFEMSolid::predictPressure
(
    const label patchID, 
    const label zoneID
)
{
    // Predict pressure field on patch
    if
    (
        D_.boundaryField()[patchID].type()
     != tractionDisplacementFvPatchVectorField::typeName
    )
    {
        FatalErrorIn("void DyParaFEMSolid::setTraction(...)")
            << "Bounary condition on " << D_.name() 
                <<  " is " 
                << D_.boundaryField()[patchID].type() 
                << "for patch" << mesh().boundary()[patchID].name()
                << ", instead " 
                << tractionDisplacementFvPatchVectorField::typeName
                << abort(FatalError);
    }

    tractionDisplacementFvPatchVectorField& patchUo =
        refCast<tractionDisplacementFvPatchVectorField>
        (
            D_.oldTime().boundaryField()[patchID]
        );

    tractionDisplacementFvPatchVectorField& patchUoo =
        refCast<tractionDisplacementFvPatchVectorField>
        (
            D_.oldTime().oldTime().boundaryField()[patchID]
        );

    scalarField pPF = 2*patchUo.pressure() - patchUoo.pressure();

    tmp<scalarField> tpF
    (
        new scalarField(mesh().faceZones()[zoneID].size(), 0)
    );
    scalarField& pF = tpF();

    const label patchStart = 
        mesh().boundaryMesh()[patchID].start();

    forAll(pPF, i)
    {
        pF[mesh().faceZones()[zoneID].whichFace(patchStart + i)] = pPF[i];
    }

    // Parallel data exchange: collect pressure field on all processors
    reduce(pF, sumOp<scalarField>());

    return tpF;
}

bool DyParaFEMSolid::evolve()
{
    Info << "Evolving solid solver: " 
        << DyParaFEMSolid::typeName << endl;

    label numPoints = mesh().points().size();
    double time = mesh().time().value();

    label tmp = Pstream::myProcNo();
    reduce(tmp,sumOp<label>());
 
    #include "updateForce.H"

    vectorField& oldPointDI = pointD_.oldTime().internalField();
    vectorField& oldPointUI = pointU_.oldTime().internalField();
    vectorField& oldPointAI = pointA_.oldTime().internalField();
     
    forAll(oldPointDI, pointI)
    {
	 d_OF_[pointI*ndim + 0]=oldPointDI[pointI].x();
	 d_OF_[pointI*ndim + 1]=oldPointDI[pointI].y();
	 d_OF_[pointI*ndim + 2]=oldPointDI[pointI].z();

	 u_OF_[pointI*ndim + 0]=oldPointUI[pointI].x();
	 u_OF_[pointI*ndim + 1]=oldPointUI[pointI].y();
	 u_OF_[pointI*ndim + 2]=oldPointUI[pointI].z();

	 a_OF_[pointI*ndim + 0]=oldPointAI[pointI].x();
	 a_OF_[pointI*ndim + 1]=oldPointAI[pointI].y();
	 a_OF_[pointI*ndim + 2]=oldPointAI[pointI].z();
     } 

    // Print Force to Parafem 
//    checkforce_
//    (
//        fext_OF_,
//        sense_, 
//        nodeensi_,
//        &numFixedForceNodes_,
//        &numPoints
//    );

    reduce(tmp,sumOp<label>());

    // Inputs to Parafem 
    runparafem_
    (
	numSchemes_,
        fext_OF_,
        forceNodes_,
        &numFixedForceNodes_,
        d_OF_,
	u_OF_,
	a_OF_,
	&time,	
        &numPoints,
        g_g_pp_OF_,
        g_num_pp_OF_,
        store_km_pp_OF_,
	store_mm_pp_OF_,
        diag_precon_pp_OF_,
	gravlo_
    );

    reduce(tmp,sumOp<label>());
   // ParaFEM Outputs from both processors the full
   // ptD, ptU and ptA arrays
    vectorField& pointDI = pointD_.internalField();
    vectorField& pointUI = pointU_.internalField();
    vectorField& pointAI = pointA_.internalField();

    // EDIT
    // Change to memcpy or std::copy()
    forAll(pointDI, pointI)
    {
	 pointDI[pointI].x() = d_OF_[pointI*ndim + 0];
	 pointDI[pointI].y() = d_OF_[pointI*ndim + 1];
	 pointDI[pointI].z() = d_OF_[pointI*ndim + 2];

	 pointUI[pointI].x() = u_OF_[pointI*ndim + 0];
	 pointUI[pointI].y() = u_OF_[pointI*ndim + 1];
	 pointUI[pointI].z() = u_OF_[pointI*ndim + 2];

	 pointAI[pointI].x() = a_OF_[pointI*ndim + 0];
	 pointAI[pointI].y() = a_OF_[pointI*ndim + 1];
	 pointAI[pointI].z() = a_OF_[pointI*ndim + 2];
    }

    D_ = pointToVol_.interpolate(pointD_);
//    U_ = pointToVol_.interpolate(pointU_); 
//    A_ = pointToVol_.interpolate(pointA_);
    return true;
}

tmp<volScalarField> DyParaFEMSolid::hydPressure() const
{
    tmp<volScalarField> tHydPressure
    (
        new volScalarField
        (
            tr(sigma_)/3
        )
    );

    return tHydPressure;
}

void DyParaFEMSolid::predict()
{
    Info << "Predicting solid:" << endl;
    D_ = D_ + U_*runTime().deltaT();
}

void DyParaFEMSolid::updateFields()
{

    Info << "updateFields:" << endl;
    volToPoint_.interpolate(D_, pointD_);

    rho_ = rheology_.rho();
//    mu_ = rheology_.mu();
//    lambda_ = rheology_.lambda();
}

tmp<surfaceVectorField> DyParaFEMSolid::traction() const
{
    tmp<surfaceVectorField> tTraction
    (
        new surfaceVectorField
        (
            (mesh().Sf() & fvc::interpolate(sigma_))
	   /mesh().magSf()
        )
    );

    return tTraction;
}

bool DyParaFEMSolid::writeObject
(
    IOstream::streamFormat,
    IOstream::versionNumber,
    IOstream::compressionType
) const
{

Switch moveMesh(solidProperties().lookup("moveMesh"));

//    nonLinearGeometry::nonLinearType nonLinear =
//        nonLinearGeometry::nonLinearNames_.read
//        (
//            solidProperties().lookup("nonLinear")
//        );
//     Switch nonLinear(solidProperties().lookup("nonLinear"));

    if (moveMesh)
    {

        pointIOField curPoints
        (
            IOobject
            (
                "points",
                runTime().timeName(),
                polyMesh::meshSubDir,
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh().allPoints()
        );

        const vectorField& pointDI = pointD_.internalField();

        forAll (pointDI, pointI)
        {
            curPoints[pointI] += pointDI[pointI];
        }

        // Unused points (procedure developed by Philip Cardiff, UCD) 
        forAll(globalFaceZones(), zoneI)
        {
            const label curZoneID = globalFaceZones()[zoneI];

            const labelList& curMap = 
                globalToLocalFaceZonePointMap()[zoneI];

            const labelList& curZoneMeshPoints =
                mesh().faceZones()[curZoneID]().meshPoints();

            vectorField curGlobalZonePointDispl
            (
                curZoneMeshPoints.size(), 
                vector::zero
            );

            //-Inter-proc points are shared by multiple procs
            // pointNumProc is the number of procs which a point lies on
            scalarField pointNumProcs(curZoneMeshPoints.size(), 0);

            forAll(curGlobalZonePointDispl, globalPointI)
            {
                label localPoint = curMap[globalPointI];
	    
                if(curZoneMeshPoints[localPoint] < mesh().nPoints())
                {
                    label procPoint = curZoneMeshPoints[localPoint];
                
                    curGlobalZonePointDispl[globalPointI] = pointDI[procPoint];
                
                    pointNumProcs[globalPointI] = 1;
                }
            }

            if (Pstream::parRun())
            {
                reduce(curGlobalZonePointDispl, sumOp<vectorField>());
                reduce(pointNumProcs, sumOp<scalarField>());

                //- now average the displacement between all procs
                curGlobalZonePointDispl /= pointNumProcs;
            }

            //- The curZonePointsDisplGlobal now contains the correct 
            //  face zone displacement in a global master processor order, 
            //  now convert them back into the local proc order

            vectorField curZonePointDispl
            (
                curZoneMeshPoints.size(), 
                vector::zero
            );

            forAll(curGlobalZonePointDispl, globalPointI)
            {
                label localPoint = curMap[globalPointI];

                curZonePointDispl[localPoint] = 
                    curGlobalZonePointDispl[globalPointI];
            }

            forAll(curZonePointDispl, pointI)
            {
                // unused points
                if (curZoneMeshPoints[pointI] >= mesh().nPoints())
                {
                    curPoints[curZoneMeshPoints[pointI]] +=
                        curZonePointDispl[pointI];
               }
            }
        }
        
        twoDPointCorrector twoDCorrector(mesh());
        twoDCorrector.correctPoints(curPoints);

        curPoints.write();

    }
    volScalarField sigmaEq
    (
        IOobject
        (
            "sigmaEq",
            runTime().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt((3.0/2.0)*magSqr(dev(sigma_)))
    );
    sigmaEq.write();

    Info<< "Max sigmaEq = " << max(sigmaEq).value()
        << endl;

    Info<< "SigmaEq, max: " << gMax(sigmaEq.internalField()) 
        << ", avg: " << gAverage(sigmaEq.internalField()) 
        << ", min: " << gMin(sigmaEq.internalField()) << endl;


    // Write point sigma field
    if (false)
    {
        pointSymmTensorField pointSigma
        (
            IOobject
            (
                "pointSigam",
                runTime().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pMesh_,
            dimensioned<symmTensor>("0", sigma_.dimensions(), symmTensor::zero)
        );

        for (direction cmpt = 0; cmpt < symmTensor::nComponents; cmpt++)
        {
            volScalarField cmptSigma = sigma_.component(cmpt);

            pointScalarField cmptPointSigma
            (
                IOobject
                (
                    "cmptPointSigma",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ
                ),
                pMesh_,
                dimensioned<scalar>
                (
                    "0",
                    sigma_.dimensions(),
                    0
                )
            );

            volToPoint_.interpolate(cmptSigma, cmptPointSigma);

            pointSigma.internalField().replace
            (
                cmpt, 
                cmptPointSigma.internalField()
            );            
        }

        pointSigma.write();
    }
   
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidSolvers
} // End namespace Foam

// ************************************************************************* //
