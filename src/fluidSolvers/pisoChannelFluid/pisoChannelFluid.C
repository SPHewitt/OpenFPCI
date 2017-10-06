/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "pisoChannelFluid.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "findRefCell.H"
#include "adjustPhi.H"
#include "fluidSolidInterface.H"
#include "fixedGradientFvPatchFields.H"
#include "IFstream.H"
#include "OFstream.H"
#include "argList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fluidSolvers
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pisoChannelFluid, 0);
addToRunTimeSelectionTable(fluidSolver, pisoChannelFluid, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pisoChannelFluid::pisoChannelFluid(const fvMesh& mesh)
:
    fluidSolver(this->typeName, mesh),
    U_
    (
        IOobject
        (
            "U",
            runTime().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    p_
    (
        IOobject
        (
            "p",
            runTime().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    gradp_(fvc::grad(p_)),
    gradU_(fvc::grad(U_)),
    phi_
    (
        IOobject
        (
            "phi",
            runTime().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(U_) & mesh.Sf()
    ),
    laminarTransport_(U_, phi_),
    turbulence_
    (
        incompressible::turbulenceModel::New
        (
            U_, phi_, laminarTransport_
        )
    ),
    rho_
    (
        IOdictionary
        (
            IOobject
            (
                "transportProperties",
                runTime().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).lookup("rho")
    ),
    Ubar_
    (
        IOdictionary
        (
            IOobject
            (
                "transportProperties",
                runTime().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).lookup("Ubar")
    ),
     gradPf_
    (        
        "gradPf",
        dimensionSet(0, 1, -2, 0, 0),
        0.0
    ),
    adjustTimeStep_(false),
    maxCo_(1.0),
    maxDeltaT_(GREAT),
    writePFlux_(false),
    writeTime_("EmptyWord"),
    gradPfOld_
    (
        "gradPfOld",
	dimensionSet(0, 1, -2, 0, 0),
	0.
    )
{
     // Initialise driving force
    IFstream gradPFile
    (
        runTime().path()/runTime().timeName()/"uniform"/"gradPf.raw"
    );

    if(gradPFile.good())
    {
        gradPFile >> gradPf_;
        Info<< "Reading average pressure gradient: " << gradPf_.value()
            << endl;
    }
    else
    {
        Info<< "Initializing with 0 pressure gradient" <<endl
            << endl;
    }; 

    // Intialise read Fields
    adjustTimeStep_ =
        runTime().controlDict().lookupOrDefault("adjustTimeStep", false);

    maxCo_ =
        runTime().controlDict().lookupOrDefault<scalar>("maxCo", 1.0);

    maxDeltaT_ =
        runTime().controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const volVectorField& pisoChannelFluid::U() const
{
    return U_;
}


const volScalarField& pisoChannelFluid::p() const
{
    return p_;
}


//- Patch viscous force (N/m2)
tmp<vectorField> pisoChannelFluid::patchViscousForce(const label patchID) const
{
    tmp<vectorField> tvF
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

    tvF() =
        rho_.value()
       *(
            mesh().boundary()[patchID].nf()
          & turbulence_->devReff()().boundaryField()[patchID]
        );

//     vectorField n = mesh().boundary()[patchID].nf();
//     tvF() -= n*(n&tvF());

    return tvF;
}

//- Patch pressure force (N/m2)
tmp<scalarField> pisoChannelFluid::patchPressureForce(const label patchID) const
{
    tmp<scalarField> tpF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    tpF() = rho_.value()*p().boundaryField()[patchID];

    return tpF;
}

//- Patch viscous force (N/m2)
tmp<vectorField> pisoChannelFluid::faceZoneViscousForce
(
    const label zoneID,
    const label patchID
) const
{
    vectorField pVF = patchViscousForce(patchID);

    tmp<vectorField> tvF
    (
        new vectorField(mesh().faceZones()[zoneID].size(), vector::zero)
    );
    vectorField& vF = tvF();

    const label patchStart =
        mesh().boundaryMesh()[patchID].start();

    forAll(pVF, i)
    {
        vF[mesh().faceZones()[zoneID].whichFace(patchStart + i)] =
            pVF[i];
    }

    // Parallel data exchange: collect pressure field on all processors
    reduce(vF, sumOp<vectorField>());


    return tvF;
}

//- Patch pressure force (N/m2)
tmp<scalarField> pisoChannelFluid::faceZonePressureForce
(
    const label zoneID,
    const label patchID
) const
{
    scalarField pPF = patchPressureForce(patchID);

    tmp<scalarField> tpF
    (
        new scalarField(mesh().faceZones()[zoneID].size(), 0)
    );
    scalarField& pF = tpF();

    const label patchStart =
        mesh().boundaryMesh()[patchID].start();

    forAll(pPF, i)
    {
        pF[mesh().faceZones()[zoneID].whichFace(patchStart + i)] =
            pPF[i];
    }

    // Parallel data exchange: collect pressure field on all processors
    reduce(pF, sumOp<scalarField>());

    return tpF;
}

tmp<scalarField> pisoChannelFluid::faceZoneMuEff
(
    const label zoneID,
    const label patchID
) const
{
    scalarField pMuEff =
        rho_.value()*turbulence_->nuEff()().boundaryField()[patchID];

    tmp<scalarField> tMuEff
    (
        new scalarField(mesh().faceZones()[zoneID].size(), 0)
    );
    scalarField& muEff = tMuEff();

    const label patchStart =
        mesh().boundaryMesh()[patchID].start();

    forAll(pMuEff, i)
    {
        muEff[mesh().faceZones()[zoneID].whichFace(patchStart + i)] =
            pMuEff[i];
    }

    // Parallel data exchange: collect pressure field on all processors
    reduce(muEff, sumOp<scalarField>());

    return tMuEff;
}

void pisoChannelFluid::evolve()
{
    Info << "Evolving Fluid Solver: "<< this->type() << endl;

    const fvMesh& mesh = fluidSolver::mesh();

    int nCorr(readInt(fluidProperties().lookup("nCorrectors")));

    int nNonOrthCorr =
        readInt(fluidProperties().lookup("nNonOrthogonalCorrectors"));

    // Prepare for the pressure solution
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p_, fluidProperties(), pRefCell, pRefValue);

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi_, U_);

    #include "courantNoFsi.H"

    dimensionedScalar magUbar = mag(Ubar_);
    vector flowDirection = (Ubar_/magUbar).value();

    fvVectorMatrix UEqn
    (
        fvm::ddt(U_)
      + fvm::div(phi_, U_)
      + turbulence_->divDevReff()
      ==
      flowDirection*gradPf_
    );

    solve(UEqn == -gradp_);

    // --- PISO loop
    volScalarField rUA = 1.0/UEqn.A();

    for (int corr=0; corr<nCorr; corr++)
    {
        U_ = rUA*UEqn.H();
        phi_ = (fvc::interpolate(U_) & mesh.Sf());
            // + fvc::ddtPhiCorr(rUA, U_, phi_);

        // Non-orthogonal pressure corrector loop
        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix pEqn
            (
                fvm::laplacian(rUA, p_) == fvc::div(phi_)
            );

            pEqn.setReference(pRefCell, pRefValue);

            if
            (
                corr == nCorr-1 && nonOrth == nNonOrthCorr
            )
            {
                pEqn.solve(mesh.solutionDict().solver("pFinal"));
            }
            else
            {
                pEqn.solve();
            }

            if (nonOrth == nNonOrthCorr)
            {
                phi_ -= pEqn.flux();
            }
        }

        #include "continuityErrsFsi.H"

        // Make the fluxes relative to the mesh motion
        fvc::makeRelative(phi_, U_);

        gradp_ = fvc::grad(p_);
        U_ -= rUA*gradp_;
        U_.correctBoundaryConditions();

        //gradU_ = fvc::grad(U_);
    }

    // Correct Driving force for a constant mass flow rate
    // Extract the velocity in the flow direction
    dimensionedScalar magUbarStar =
    (flowDirection & U_)().weightedAverage(mesh.V());

    // Calculate the pressure gradient increment needed to
    // adjust the average flow-rate to the correct value
    dimensionedScalar gragPplus =
        (magUbar - magUbarStar)/rUA.weightedAverage(mesh.V());   

    U_ += flowDirection*rUA*gragPplus;

    gradPf_ += gragPplus;

    Info<< "Uncorrected Ubar = " << magUbarStar.value() << tab
    << "pressure gradient = " << gradPf_.value() << endl;

    // Write happens dt after runtime.write()
    #include "writeGradPf.H"


    turbulence_->correct();

    // Make the fluxes absolute to the mesh motion
    fvc::makeAbsolute(phi_, U_);
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolvers
} // End namespace Foam

// ************************************************************************* //
