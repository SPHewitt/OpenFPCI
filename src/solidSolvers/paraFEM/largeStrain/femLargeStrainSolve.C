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

#include "femLargeStrain.H"

#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGradf.H"

#include "tractionDisplacementFvPatchVectorField.H"
#include "skewCorrectionVectors.H"
#include "multiMaterial.H"
#include "twoDPointCorrector.H"

#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidSolvers
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volVectorField> femLargeStrain::residual
(
    const volVectorField& source
)
{
    tmp<volVectorField> tResidual
    (
        new volVectorField
        (
            IOobject
            (
                "residual(" + D_.name() + ")",
                runTime().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedVector("0", dimForce, vector::zero)
        )
    );

    return tResidual;
}


void femLargeStrain::initialise
(
    const volVectorField& sol,
    bool consistentBoundaryField
)
{
    Info << "Initialise solid solver solution" << endl;
}


scalar femLargeStrain::smooth
(
    const volVectorField& source,
    label nCorrectors
)
{
    Info << "Smoothing solid solver solution" << endl;
    return 0;
}


scalar femLargeStrain::smooth
(
    const volVectorField& source,
    const volVectorField& refSolution,
    label nCorrectors
)
{
    Info << "Smoothing solid solver solution" << endl;
    return 0;
}

scalar femLargeStrain::residual() const
{
    return 0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidSolvers
} // End namespace Foam

// ************************************************************************* //
