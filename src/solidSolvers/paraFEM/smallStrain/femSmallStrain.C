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

#include "femSmallStrain.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "tractionDisplacementFvPatchVectorField.H"
#include "velocityTractionDisplacementFvPatchVectorField.H"
#include "multiMaterial.H"
#include "twoDPointCorrector.H"
#include "fixedValuePointPatchFields.H"
#include "emptyPointPatchFields.H"
#include "RBFInterpolation.H"
#include "TPSFunction.H"
#include "ZoneIDs.H"
#include "primitivePatchInterpolation.H"

using namespace rbf;


// * * * * * * * * * * * * ParaFEM Fortran Subroutines* * * * * * * * * * * * //
// - Will need editing when more complex elements introduced
const int nod   =  8;        // Element type
const int ndim  =  3;        // Number of Dimensions
const int ntot  =  ndim*nod; // ntot

using namespace std;


// Declaration of Fortran Subroutines 
extern"C"
{

    // Initialises ParaFEM
    // Called at construction
    void initl_
    (
        double* g_coord,
        int* rest,
        int* nn, 
        const int* nr,
        int* g_num_pp,
        int* g_g_pp,
        double* g_coord_pp
    );

    // return number of equations/proc 
    int findneqpp_();

    // return number of cells/proc
    int setnelspp_
    (
        const int* numCells
    );

    // Solve the structural equation
    void runl_
    (
        int* node,
        double* val, 
        double* numVar,
        double* matProp,
        int* nr,
        int* loadedNodes,
        double* time,
        int* g_g_pp,
        int* g_num_pp,
        double* g_coord_pp,
        double* gravlo,
        double* ptDtemp_,
        double* ptUtemp_,
        double* ptAtmep_,
        double* flag
    );

    // Calculate Gravitational Loads 
    void gloads_
    (
        double* gravlo,
        double* gravity,
        int* nn,
        int* nodof,
        const int* nod,
        const int* ndim,
        int* nr,
        double* g_coord,
        int* g_num_pp,
        int* rest
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidSolvers
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(femSmallStrain, 0);
addToRunTimeSelectionTable(solidSolver, femSmallStrain, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

femSmallStrain::femSmallStrain(const fvMesh& mesh)
:
    solidSolver(typeName, mesh),
    mPoints_(NULL),
    solidProps_(NULL),
    numSchemes_(NULL),
    g_num_pp_OF_(NULL),
    g_g_pp_OF_(NULL),
    g_coord_pp_OF_(NULL),
    gravlo_(NULL),
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
    timeflag_(0),
    gCells_(0),
    gPoints_(0),
    ptDtemp_(0),
    ptUtemp_(0),
    ptAtemp_(0),
    rbfUpdate_(false),
    twoDimensional_(false),
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
            //IOobject::READ_IF_PRESENT,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh_
	//dimensionedVector("0", dimVelocity, vector::zero)
    ),
    pointA_
    (
        IOobject
        (
            "pointA",
            runTime().timeName(),
            mesh,
            //IOobject::READ_IF_PRESENT,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh_
	//dimensionedVector("0", dimAcceleration, vector::zero)
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
    pointProcAddressing_
    (
        IOobject
        (
            "pointProcAddressing",
             mesh.facesInstance(),
             mesh.meshSubDir,
             mesh,
             IOobject::READ_IF_PRESENT,
             IOobject::NO_WRITE
        )
    ),
    cellProcAddressing_
    (
        IOobject
        (
            "cellProcAddressing",
            mesh.facesInstance(),
            mesh.meshSubDir,
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    of2pfmap_(mesh.points().size()),
    rheology_(sigma_, D_),
    volToPoint_(mesh),
    pointToVol_(pMesh_,mesh),
    rho_("rho", rheology_.rho()),
    mu_(rheology_.mu()),
    lambda_(rheology_.lambda()),
    interface_(NULL)
{
    pointD_.oldTime();
    pointU_.oldTime();
    pointA_.oldTime();

    if (rheology_.law().type() == multiMaterial::typeName)
    {
        interface_.set(new TLMaterialInterface(D_, pointD_));
    }

    if (solidProperties().found("rbfUpdate"))
    {
        rbfUpdate_ = Switch(solidProperties().lookup("rbfUpdate"));
    }

    Info << "rbfUpdate: " << rbfUpdate_ << endl;

    // - ParaFEM construction 
    E_ 		=  rheology_.law().E()()[0];
    nu_ 	=  rheology_.law().nu()()[0];
    rhotmp_ 	=  rheology_.law().rho()()[0];

    // - Turn on/off the solid mesh Motion
    Switch moveMesh(solidProperties().lookup("moveMesh"));

    
//------------------------------------------------------------------------------
//  ParaFEM: Create Steering and Coordinate Matricies
//------------------------------------------------------------------------------

    // Find globalPoints and local Points(note can be done with global MESH)
    if(Pstream::parRun()==true)
    {
	// Global Mesh Points
	label index 	=  findMax(pointProcAddressing_);
	gPoints_ 	=  pointProcAddressing_[index];

	reduce(gPoints_,maxOp<label>());

	gPoints_ 	=  gPoints_ + 1; // Fields Start at 0
	index 		=  findMax(cellProcAddressing_);
	gCells_ 	=  cellProcAddressing_[index];

	reduce(gCells_,maxOp<label>());

	gCells_ 	=  gCells_ + 1; // Fields Start at 0
	mPoints_ 	=  new double[gPoints_*ndim];
    }
    else
    {
        gPoints_ 	=  mesh.points().size();
	gCells_ 	=  mesh.nCells();
	mPoints_ 	=  new double [gPoints_*ndim];
    }

    // Populate mPoints
    pointIOField startPoints
    (
        IOobject
        (
    	    "points",
             runTime().caseConstant(),
             polyMesh::meshSubDir,
             mesh,
             IOobject::MUST_READ,
             IOobject::NO_WRITE
        )
    );
    	
    label globalIndex = 0;

    forAll(startPoints, pointI)
    {
	mPoints_[globalIndex++]  =  startPoints[pointI].x();
	mPoints_[globalIndex++]  =  startPoints[pointI].y();		
	mPoints_[globalIndex++]  =  startPoints[pointI].z();
    }
    	
    // nels_pp
    const int nels_pp_OF = mesh.nCells();
    
    // Set nels_pp
    setnelspp_(&nels_pp_OF);

    g_num_pp_OF_ = new int [nod*nels_pp_OF]; 
    ptDtemp_     = new double [ntot*nels_pp_OF];
    ptUtemp_     = new double [ntot*nels_pp_OF];
    ptAtemp_     = new double [ntot*nels_pp_OF];

    // cellPoints() returns steering array
    const labelListList& cellPoints = mesh.cellPoints();

    label localIndex = 0;
 
    if(Pstream::parRun()==true)
    {
	forAll(cellPoints, cellI)
	{
	    const labelList& curCellPoints = cellPoints[cellI]; 

	    if (curCellPoints.size() != nod)
	    {
		Info << "Not a hex cell!" << endl;
	    }

	    for(label i=0;i<nod;i++)
	    {
		g_num_pp_OF_[localIndex++] = pointProcAddressing_[curCellPoints[i]]+1; // +1 fortran
	    }	 
	}
    }
    else
    {
        forAll(cellPoints, cellI)
        {
            const labelList& curCellPoints = cellPoints[cellI];
            
 	    if (curCellPoints.size() != nod)
            {
                Info << "Not a hex cell!" << endl;
            }

            for(label i=0;i<nod;i++)
            {
                g_num_pp_OF_[localIndex++] = curCellPoints[i]+1; // +1 fortran
            }
    	}
     }


//------------------------------------------------------------------------------
//  ParaFEM: Create Restrained Arrays
//------------------------------------------------------------------------------
// paraFem	    :  0=Restrained 1=Unrestrained
// ensi_gold	:  1=Restrained 0=Unrestrained

    numRestrNodes_  =  0;
    // First Time Get size
    forAll(pointD_.boundaryField(), patchI)
    {
        if 
        (
            isA<fixedValuePointPatchVectorField>
            (
                pointD_.boundaryField()[patchI]
            )
        )
        {
	        numRestrNodes_ += mesh.boundaryMesh()[patchI].meshPoints().size();
	    }
        
        if 
        (
            isA<emptyPointPatchVectorField>
            (
                pointD_.boundaryField()[patchI]
            )
        )
	    {
	        numRestrNodes_ += mesh.boundaryMesh()[patchI].meshPoints().size();
	        twoDimensional_ = true;
	    }	
    }
 
    if(twoDimensional_)
    {
	    Info << "Simulation is 2D" << endl;
    }
    // Declare local Restrained List    
    labelListList localRest(numRestrNodes_);

//------------------------------------------------------------------------------
//  ParaFEM: Boundary Conditions 
//------------------------------------------------------------------------------
// Special care Needs to be taken here, the current method isn't robust 

    int counter = 0;    

    forAll(pointD_.boundaryField(), patchI)
    {

        // ------ Fixed Z ------
        if 
        (
                isA<emptyPointPatchVectorField>
                (
                    pointD_.boundaryField()[patchI]
                )
        )
        {
            const labelList& mp = mesh.boundaryMesh()[patchI].meshPoints();

            forAll(mp, pI)
            {
                labelList myList(4);
                if(Pstream::parRun()==true)
                {
                    myList[0]  =  pointProcAddressing_[mp[pI]]+1;
                }
                else
                {
                    myList[0]  =  mp[pI]+1;
                }

                myList[1]  =  1;
                myList[2]  =  1;
                myList[3]  =  0;
                localRest[counter] =  myList;
                counter++;
            }
        }

        // ------ Fixed X Y Z ------
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
                labelList myList(4);
                if(Pstream::parRun()==true)
                {
                    myList[0]  =  pointProcAddressing_[mp[pI]]+1;
                }
                else
                {
                    myList[0]  =  mp[pI]+1;
                }

                if(mesh.boundary()[patchI].name() == "frontAndBack")
                {
                    myList[1]  =  1;
                    myList[2]  =  1;
                    myList[3]  =  0;
                }
                else
                {
                    myList[1]  =  0;
                    myList[2]  =  0;
                    myList[3]  =  0;
                }
                localRest[counter] =  myList;
                counter++;
            }
        }
    }

//------------------------------------------------------------------------------
//  ParaFEM: Generating Restrained Array 
//------------------------------------------------------------------------------

    // Sort lists
    sort(localRest);
    numRestrNodes_ = localRest.size();

    // Gather All Lists
    List<List<List<label>>> globalRest(Pstream::nProcs(), localRest);   
    //globalRest[Pstream::myProcNo()] = localRest;

    Pstream::gatherList(globalRest);
    Pstream::scatterList(globalRest); 

    // Create master Rest;
    reduce(numRestrNodes_,sumOp<double>());
    labelListList masterRest (numRestrNodes_);

    label restIndex = 0;
    forAll(globalRest,procI)
    {
        forAll(globalRest[procI],listI)
        {
            masterRest[restIndex]=globalRest[procI][listI];
            restIndex++;
        }
    }

    //Sort and find duplicate nodes
    sort(masterRest);
    labelList duplicateNodes;
    duplicateOrder(masterRest,duplicateNodes);

   // Recount number of restrained nodes    
    numRestrNodes_ = 0;
    forAll(masterRest,listI)
    {
        if (listI == 0)
        {
            numRestrNodes_++;
            continue;
        }

        if (masterRest[listI][0] == masterRest[listI-1][0])
        {
            // If exists continue
            continue;
        }
        else
        {
            // If doesnt exist insert into rest
            numRestrNodes_++;
        }
    }

    rest_ = new int [numRestrNodes_*4];
    restIndex = 0;

    // Set rest to ZERO
    for(int i=0; i<numRestrNodes_*4;i++)
    {
        rest_[i]=0.0;
    }

    // Populate rest array 
    forAll(masterRest,listI)
    {
        if (listI != 0 && masterRest[listI][0] == masterRest[listI-1][0])
        {
            continue;
        }
        else
        {
            rest_[ numRestrNodes_* 0 + restIndex ] =  masterRest[listI][0];
            rest_[ numRestrNodes_* 1 + restIndex ] =  masterRest[listI][1];
            rest_[ numRestrNodes_* 2 + restIndex ] =  masterRest[listI][2];
            rest_[ numRestrNodes_* 3 + restIndex ] =  masterRest[listI][3];
            restIndex++;	
        }
    }

    // Debugging
    // Most Errors occur from boundary conditions 
    if(false)
    {
        fileName outputFile("rest.txt");
        OFstream os(db().time().system()/outputFile);
        os << "Restrained Array.\n" << endl;
        os << numRestrNodes_ << "\n{" << endl;
        for(int i = 0; i < numRestrNodes_; i++)
        {
            os << rest_[numRestrNodes_* 0 + i] << " ";
            os << rest_[numRestrNodes_* 1 + i] << " ";
            os << rest_[numRestrNodes_* 2 + i] << " ";
            os << rest_[numRestrNodes_* 3 + i]  << endl;
        }
        os << "}" << endl;
    }   


//------------------------------------------------------------------------------
//  ParaFEM: Create Arrays for intialisation
//------------------------------------------------------------------------------

    g_g_pp_OF_ 	    =  new int [ntot*nels_pp_OF];
    g_coord_pp_OF_  =  new double [nod*ndim*nels_pp_OF];

    double alpha1 (readScalar(solidProperties().lookup("alpha1")));
    double beta1 (readScalar(solidProperties().lookup("beta1")));
    double theta (readScalar(solidProperties().lookup("theta")));
    
    // Tolerances    
    // Default Values
    double tol   = 1e-6;
    double limit = 1000;
    
    if (solidProperties().found("PCGTolerance"))
    {
        tol = readScalar(solidProperties().lookup("PCGTolerance"));
    }
    
    if (solidProperties().found("PCGLimit"))
    {
        limit = readScalar(solidProperties().lookup("PCGLimit"));
    }
    
    numSchemes_ 	=  new double[5];
    numSchemes_[0] 	=  alpha1;
    numSchemes_[1] 	=  beta1;
    numSchemes_[2] 	=  theta;
    numSchemes_[3] 	=  tol;
    numSchemes_[4] 	=  limit;

    solidProps_ 	=  new double[3];
    solidProps_[0] 	=  E_;
    solidProps_[1] 	=  nu_;
    solidProps_[2] 	=  rhotmp_;

    label tmp = Pstream::myProcNo();
    reduce(tmp,sumOp<label>());

//------------------------------------------------------------------------------
//  ParaFEM: Initialise ParaFEM
//------------------------------------------------------------------------------

    initl_
    (
        mPoints_,
        rest_,
        &gPoints_,
        &numRestrNodes_,
        g_num_pp_OF_,
        g_g_pp_OF_,
        g_coord_pp_OF_
    );

    reduce(tmp,sumOp<label>());

//------------------------------------------------------------------------------
//  ParaFEM: Find Gravitry loading if set in dictionary
//------------------------------------------------------------------------------

    const int neq_pp_OF = findneqpp_();
    double gravity(readScalar(solidProperties().lookup("gravity")));
    gravlo_ = new double [neq_pp_OF];

    if(gravity > 1e-6)
    { 
	Info << "Gravity Loading, gravity: " << gravity << " m/s^2" << endl; 

	// Specific weight lambda = rho * g (gloads: negative y implied)
	double specWeight=gravity*rhotmp_;;
	int nodof=3;	    

	gloads_ 
	(
	    gravlo_,
	    &specWeight,
	    &gPoints_,
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
    

//------------------------------------------------------------------------------
//  ParaFEM: Create Force Arrays
//------------------------------------------------------------------------------

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


    forceNodes_   =  new int [numFixedForceNodes_*ndim];
    fext_OF_ 	  =  new double [numFixedForceNodes_*ndim];

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
  	    // Local Point
	    const labelList& mp = 
                mesh.boundaryMesh()[patchI].meshPoints();

            if (Pstream::parRun())
            {
	    	jj_=0;
	    	forAll(mp, pI)
            	{			  
		    // Global Point
		    label myP = pointProcAddressing_[mp[pI]];
	 	    forceNodes_[gi++] = myP+1;
		    jj_=jj_+3;
	    	}
	    }
	    else
	    {
	    	jj_=0;
	    	forAll(mp, pI)
                {	
		    forceNodes_[gi++] = mp[pI]+1;
 		    jj_=jj_+3;
	        }		
	    }
    	} // if
    } // forAll


    // Set force to zero
    for(int i=0; i<numFixedForceNodes_*ndim; i++)
    {
 	fext_OF_[i]=0;
    }

//------------------------------------------------------------------------------
//  ParaFEM: Create Map from OpenFOAM local to ParaFEM 
//------------------------------------------------------------------------------

    vectorField& pointDI = pointD_.internalField();

    forAll(pointDI, pointI)
	{
	    label value=0;

	    if(Pstream::parRun()==true)
    	    {
	        value  =  pointProcAddressing_[pointI]+1;
	    }
	    else
 	    {
		value  =  pointI+1;
	    }

	    label index  =  0;
	    label iel    =  0;  

	    label resizeval=0;
	    label counter=0;

	    // nod represnets the max number of elements a node may hold
	    labelList myLabel (nod,0);

	    // for each value find position in steering matrix
	    for(int i=0; i<nels_pp_OF*nod;i++)
	    {
		if(i % nod == 0 && i!=0)
		{
		    iel++;
		    index=0;
		}

		if(value==g_num_pp_OF_[i])		
		{ 
		    resizeval++;
		    myLabel[counter]=(iel*ntot)+(index*ndim);
		    counter++;
		}

			index++;
	    }	// end for loop	
	    myLabel.resize(resizeval);
	    of2pfmap_[pointI] = myLabel;

	} // end for all

    delete[] mPoints_;
    delete[] rest_ensi_;
    delete[] rest_;
    delete[] nodeensi_;
    delete[] sense_;

} // End Constructor


femSmallStrain::~femSmallStrain()
{
    delete[] solidProps_;
    delete[] g_num_pp_OF_;
    delete[] g_g_pp_OF_;
    delete[] numSchemes_;
    delete[] forceNodes_;
    delete[] fext_OF_;
    delete[] gravlo_;
    delete[] ptDtemp_;
    delete[] ptUtemp_;
    delete[] ptAtemp_;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vector femSmallStrain::pointU(label pointID) const
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
tmp<vectorField> femSmallStrain::patchPointDisplacementIncrement
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
tmp<vectorField> femSmallStrain::faceZonePointDisplacementIncrement
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
	//Info << "WARNING FIELDS ARE STORED ON EVERY PROCESSOR" << endl;
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
tmp<vectorField> femSmallStrain::patchPointDisplacement
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
tmp<vectorField> femSmallStrain::faceZonePointDisplacement
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
tmp<Foam::vectorField> femSmallStrain::patchFaceAcceleration
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
tmp<vectorField> femSmallStrain::faceZoneAcceleration
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
tmp<vectorField> femSmallStrain::faceZoneVelocity
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
tmp<tensorField> femSmallStrain::faceZoneSurfaceGradientOfVelocity
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
femSmallStrain::currentFaceZonePoints(const label zoneID) const
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
tmp<vectorField> femSmallStrain::faceZoneNormal
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

void femSmallStrain::setTraction
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
        FatalErrorIn("void femSmallStrain::setTraction(...)")
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

void femSmallStrain::setPressure
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
        FatalErrorIn("void femSmallStrain::setTraction(...)")
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

void femSmallStrain::setTraction
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



void femSmallStrain::setPressure
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
void femSmallStrain::setVelocityAndTraction
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
            "void femSmallStrain::"
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
void femSmallStrain::setVelocityAndTraction
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

tmp<vectorField> femSmallStrain::predictTraction
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
        FatalErrorIn("void femSmallStrain::setTraction(...)")
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


tmp<scalarField> femSmallStrain::predictPressure
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
        FatalErrorIn("void femSmallStrain::setTraction(...)")
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

bool femSmallStrain::evolve()
{
    Info << "Evolving solid solver: " 
        << femSmallStrain::typeName << endl;
        
    double femFlag = 1;
        
    if(timeflag_ == runTime().value())
    {
       femFlag = 1;
    }
    else //Next time Step
    {
       timeflag_= runTime().value();
       femFlag = 2;
    }
        

    double dtim = runTime().deltaT().value();

    label tmp = Pstream::myProcNo();
    reduce(tmp,sumOp<label>());

    // Interpolating Face to Point Forces
    #include "updateForce.H"
 
    vectorField& oldPointDI = pointD_.oldTime().internalField();
    vectorField& oldPointUI = pointU_.oldTime().internalField();
    vectorField& oldPointAI = pointA_.oldTime().internalField();

    // Copy data OF-PF
    forAll(oldPointDI, pointI)
    {
	
	forAll(of2pfmap_[pointI],value)
	{
	    ptDtemp_[of2pfmap_[pointI][value] + 0]  =  oldPointDI[pointI].x();
	    ptDtemp_[of2pfmap_[pointI][value] + 1]  =  oldPointDI[pointI].y();
	    ptDtemp_[of2pfmap_[pointI][value] + 2]  =  oldPointDI[pointI].z();

	    ptUtemp_[of2pfmap_[pointI][value] + 0]  =  oldPointUI[pointI].x();
	    ptUtemp_[of2pfmap_[pointI][value] + 1]  =  oldPointUI[pointI].y();
	    ptUtemp_[of2pfmap_[pointI][value] + 2]  =  oldPointUI[pointI].z();

	    ptAtemp_[of2pfmap_[pointI][value] + 0]  =  oldPointAI[pointI].x();
	    ptAtemp_[of2pfmap_[pointI][value] + 1]  =  oldPointAI[pointI].y();
	    ptAtemp_[of2pfmap_[pointI][value] + 2]  =  oldPointAI[pointI].z(); 
	}
	
    }

    
//------------------------------------------------------------------------------
//  ParaFEM: Run ParaFEM Code
//------------------------------------------------------------------------------
    runl_
    (
        forceNodes_,
        fext_OF_,
        numSchemes_,
        solidProps_,
        &numRestrNodes_,
        &numFixedForceNodes_,
        &dtim,
        g_g_pp_OF_,
        g_num_pp_OF_,
        g_coord_pp_OF_,
        gravlo_,
        ptDtemp_,
        ptUtemp_,
        ptAtemp_,
        &femFlag
    );


    vectorField& pointDI = pointD_.internalField();
    vectorField& pointUI = pointU_.internalField();
    vectorField& pointAI = pointA_.internalField();

    // Copy data PF-OF
    forAll(pointDI, pointI)
    {
        pointDI[pointI].x() = ptDtemp_[of2pfmap_[pointI][0] + 0];
        pointDI[pointI].y() = ptDtemp_[of2pfmap_[pointI][0] + 1];
        pointDI[pointI].z() = ptDtemp_[of2pfmap_[pointI][0] + 2];

        pointUI[pointI].x() = ptUtemp_[of2pfmap_[pointI][0] + 0];
        pointUI[pointI].y() = ptUtemp_[of2pfmap_[pointI][0] + 1];
        pointUI[pointI].z() = ptUtemp_[of2pfmap_[pointI][0] + 2];

        pointAI[pointI].x() = ptAtemp_[of2pfmap_[pointI][0] + 0];
        pointAI[pointI].y() = ptAtemp_[of2pfmap_[pointI][0] + 1];
        pointAI[pointI].z() = ptAtemp_[of2pfmap_[pointI][0] + 2];
    }
    

    // Calculate Cauchy Green Stress Tensor
    {
	D_ = pointToVol_.interpolate(pointD_);
        epsilon_ = symm(fvc::grad(D_));

        sigma_ = 2*mu_*epsilon_ + I*(lambda_*tr(epsilon_));
    }

    //    U_ = pointToVol_.interpolate(pointU_); 
    //    A_ = pointToVol_.interpolate(pointA_);
    return true;
}

tmp<volScalarField> femSmallStrain::hydPressure() const
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

void femSmallStrain::predict()
{
    Info << "Predicting solid:" << endl;
    D_ = D_ + U_*runTime().deltaT();
}

void femSmallStrain::updateFields()
{
    Info << "updateFields:" << endl;

    rho_ = rheology_.rho();
    mu_ = rheology_.mu();
    lambda_ = rheology_.lambda();
}

tmp<surfaceVectorField> femSmallStrain::traction() const
{
    Info << "traction:" << endl;
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

bool femSmallStrain::writeObject
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

	// startPoints is the reference coordinates of 
	// initial configuration

	pointIOField startPoints
        (
            IOobject
            (
                "points",
                runTime().caseConstant(),
                polyMesh::meshSubDir,
                mesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

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
	   // curPoints[pointI]=startPoints[pointI]+pointDI[pointI];
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
