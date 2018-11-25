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

Class
    femLargeStrain

Description
   The femLargeStrain solver is a wrapper to integrate the ParaFEM fortran
   routines into Foam-Extend. The routines solve the deformation for a
   linear elastic solid.

Authors
    Sam Hewitt, University of Manchester.
    Zeljko Tukovic, FSB Zagreb.

SourceFiles
    femLargeStrain.C

\*---------------------------------------------------------------------------*/

#ifndef femLargeStrain_H
#define femLargeStrain_H

#include "solidSolver.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "constitutiveModel.H"
#include "leastSquaresVolPointInterpolation.H"
#include "pointVolInterpolation.H"
#include "TLMaterialInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidSolvers
{

/*---------------------------------------------------------------------------*\
                        Class femLargeStrain Declaration
\*---------------------------------------------------------------------------*/

class femLargeStrain
:
    public solidSolver
{
    // Private data
  
        //- An array for the StressMesh Points Data
        double* mPoints_;

        //- An array for the Solid Properties
        double* solidProps_;

        //- An array for the Numerical Schemes
        double* numSchemes_;

        //- Steering matrix
        int* g_num_pp_OF_;

        //- Global neq matrix
        int*  g_g_pp_OF_; 

        //- Coordinate matrix per processor
         double* g_coord_pp_OF_;

        //-
        int jj_;


	    //- bodyForce
	    double* gravlo_;

        //- Number of restrained nodes
        int numRestrNodes_;

        //- restrained array
        int* rest_;

        //- restrained array ensi gold
        int* rest_ensi_;

        //- Number of force specified nodes
        int numFixedForceNodes_;

        //- Node numbers for interface
        int* forceNodes_;

        //- Global list of force specified nodes
        int* globalForceNodes_;

        //- Mapping between Force local and global nodes
        int* forceLocalGlobalMap_;
    
        //- Total number of global forced nodes
        int numGlobalForceNodes_;

        //- processor counter
        int* processorCount_;

        //- Force val at interface
        double* fext_OF_;

        //-
        int* nodeensi_;

        //-
        int* sense_;

        //- Youngs Modulus
        double E_;

        //- Poissons ratio
        double nu_;

        //- temporary Density
        double rhotmp_;

        //- Total Cells
        int gCells_;

        //- Total Points
        int gPoints_;

	    //- Displacement Field PF format
	    double* ptDtemp_;

	    //- Velocity Field PF format
        double* ptUtemp_;

        //- Acceleration Field PF format
        double* ptAtemp_;

        //- Use rbf interpolation
        Switch rbfUpdate_;

        //- Use rbf interpolation
        Switch primitivePatchUpdate_;
        
        //- Check if Simulation 2D
        Switch twoDimensional_;

        //- Displacement field
        volVectorField D_;

        //- Velocity field
        volVectorField U_;

        //- Point mesh
        pointMesh pMesh_;

        //- Point displacement field
        pointVectorField pointD_;

        //- Point velocity field
        pointVectorField pointU_;

        //- Point acceleration field
        pointVectorField pointA_;

        //- Stress field
        volSymmTensorField sigma_;

        //- Strain field
        volSymmTensorField epsilon_;

        //- Point Processor Addressing 
        labelIOList pointProcAddressing_;

        //- Cell Processor Addressing
        labelIOList cellProcAddressing_;

        //- An array for global to Local point Map (evolve Solid)
        labelListList of2pfmap_;
	
        //- Rheology model
        constitutiveModel rheology_;

        //- Cell-to-point interpolation
        leastSquaresVolPointInterpolation volToPoint_;

	// Point-to-cell interpolation
	pointVolInterpolation pointToVol_;

        //- Density
        volScalarField rho_;

        //- First Lame's coefficient
        volScalarField mu_;

        //- Second Lame's coefficient
        volScalarField lambda_;

        //- Solid interface
        autoPtr<materialInterface> interface_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        femLargeStrain(const femLargeStrain&);

        //- Disallow default bitwise assignment
        void operator=(const femLargeStrain&);


public:

    //- Runtime type information
    TypeName("femLargeStrain");

    // Constructors

        //- Construct from components
        femLargeStrain(const fvMesh& mesh);


    // Destructor

        virtual ~femLargeStrain();



    // Member Functions

        // Access

            //- Return cell displacement field
            const volVectorField& D() const
            {
                return D_;
            }

            //- Return cell displacement field
            volVectorField& D()
            {
                return D_;
            }

            //- Return point mesh
            const pointMesh& pMesh() const
            {
                return pMesh_;
            }

            const volSymmTensorField& sigma() const
            {
                return sigma_;
            }

            const volSymmTensorField& epsilon() const
            {
                return epsilon_;
            }

            //- Return multi-material interface pointer
            virtual const autoPtr<materialInterface>& interface() const
            {
                return interface_;
            }

            //- Return multi-material interface pointer
            virtual autoPtr<materialInterface>& interface()
            {
                return interface_;
            }

            virtual const constitutiveModel& rheology() const
            {
                return rheology_;
            }

            const volScalarField& mu() const
            {
                return mu_;
            }

            const volScalarField& lambda() const
            {
                return lambda_;
            }

            //- Return point displacement field
            virtual const pointVectorField& pointD() const
            {
                return pointD_;
            }

            //- Return multi-material interface pointer
            const leastSquaresVolPointInterpolation& volToPoint() const
            {
                return volToPoint_;
            }

	    //- Return multi-material interface pointer
	    const pointVolInterpolation& pointToVol() const
            {
                return pointToVol_;
            }

            //- Return point displacement field
            virtual vector pointU(label pointID) const;

            //- Patch point displacement increment
            virtual tmp<vectorField> patchPointDisplacementIncrement
            (
                const label patchID
            ) const;

            //- Face zone point displacement increment
            virtual tmp<vectorField> faceZonePointDisplacementIncrement
            (
                const label zoneID
            ) const;

            //- Patch point displacement increment
            virtual tmp<vectorField> patchPointDisplacement
            (
                const label patchID
            ) const;

            //- Face zone point displacement increment
            virtual tmp<vectorField> faceZonePointDisplacement
            (
                const label zoneID
            ) const;

            //- Patch face acceleration
            virtual tmp<vectorField> patchFaceAcceleration
            (
                const label patchID
            ) const;

            //- Face zone velocity
            virtual tmp<vectorField> faceZoneAcceleration
            (
                const label zoneID,
                const label patchID
            ) const;

            //- Face zone velocity
            virtual tmp<vectorField> faceZoneVelocity
            (
                const label zoneID,
                const label patchID
            ) const;

            //- Face zone velocity surface (tangential) gradient
            virtual tmp<tensorField> faceZoneSurfaceGradientOfVelocity
            (
                const label zoneID,
                const label patchID
            ) const;

            //- Face zone current normal
            virtual tmp<vectorField> faceZoneNormal
            (
                const label zoneID,
                const label patchID
            ) const;

            //- Face zone primitive patch in current configuration
            virtual tmp<vectorField>
            currentFaceZonePoints(const label zoneID) const;

        // Edit

            //- Set traction at specified patch
            virtual void setTraction
            (
                const label patchID,
                const vectorField& traction
            );

            //- Set pressure at specified patch
            virtual void setPressure
            (
                const label patchID,
                const scalarField& pressure
            );

            //- Set traction at specified patch
            virtual void setTraction
            (
                const label patchID,
                const label zoneID,
                const vectorField& zoneTraction
            );

            //- Set pressure at specified patch
            virtual void setPressure
            (
                const label patchID,
                const label zoneID,
                const scalarField& zonePressure
            );

            //- Set traction at specified patch
            virtual void setVelocityAndTraction
            (
                const label patchID,
                const vectorField& traction,
                const vectorField& velocity,
                const vectorField& normal
            );

            //- Set traction at specified patch
            virtual void setVelocityAndTraction
            (
                const label patchID,
                const label zoneID,
                const vectorField& faceZoneTraction,
                const vectorField& faceZoneVelocity,
                const vectorField& faceZoneNormal
            );

            //- Predict traction at specified zone
            virtual tmp<vectorField> predictTraction
            (
                const label patchID,
                const label zoneID
            );

            //- Predict pressure at specified zone
            virtual tmp<scalarField> predictPressure
            (
                const label patchID,
                const label zoneID
            );

            //- Evolve the solid solver
            virtual bool evolve();

            //- Evolve the solid solver
            virtual void predict();

            //- Smooth solution
            virtual scalar smooth
            (
                const volVectorField& source,
                label nCorrectors = 0
            );

            //- Smooth solution
            virtual scalar smooth
            (
                const volVectorField& source,
                const volVectorField& refSolution,
                label nCorrectors = 0
            );

            //- Get current hydrosatic pressure
            virtual tmp<volScalarField> hydPressure() const;

            //- Return solution (displacement field)
            virtual const volVectorField& solution() const
            {
                return D_;
            }

            //- Return solution (displacement field)
            virtual volVectorField& solution()
            {
                return D_;
            }

            //- Initialise solution
            virtual void initialise
            (
                const volVectorField& sol,
                bool consistentBoundaryField = false
            );

            //- Return current residual
            virtual tmp<volVectorField> residual
            (
                const volVectorField& source
            );

            //- Return current relative momentum resudual
            scalar residual() const;

            //- Update fields
            virtual void updateFields();

            //- Get current traction force
            virtual tmp<surfaceVectorField> traction() const;

            //- Write using given format, version and compression
            virtual bool writeObject
            (
                IOstream::streamFormat,
                IOstream::versionNumber,
                IOstream::compressionType
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
