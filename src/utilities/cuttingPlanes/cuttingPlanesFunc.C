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

Author
    Sam Hewitt, University of Manchester.  All rights reserved

\*----------------------------------------------------------------------------*/
#include "fvc.H"
#include "cuttingPlanesFunc.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "OStringStream.H"
#include "IStringStream.H"
#include "IOmanip.H"
#include "plane.H"
#include "cuttingPlane.H"
#include "cellSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cuttingPlanesFunc, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        cuttingPlanesFunc,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::cuttingPlanesFunc::writeData()
{

    // Reference to Mesh
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(polyMesh::defaultRegion);

    // Create Cutting plane
    plane pl1(point_,dir_);
    cuttingPlane cutPlane(pl1,mesh);

    // Create cellSet of cut cells
    const labelList& cutCells = cutPlane.cutCells();

    // Lookup fields to write
    const volScalarField& p =
    mesh.lookupObject<volScalarField>("p");

    const volVectorField& U =
    mesh.lookupObject<volVectorField>("U");
    
    const volVectorField UMean = 
    mesh.lookupObject<volVectorField>("UMean");

    const volScalarField pMean = 
    mesh.lookupObject<volScalarField>("pMean");

    const volSymmTensorField UPrime2Mean = 
    mesh.lookupObject<volSymmTensorField>("UPrime2Mean");
   
    // Calculate vorticity
    volVectorField vorticity = fvc::curl(U);
    
    if(cutCells.size()>0)
    { 
       fileName outFile(planeName_+"_"+time_.timeName()+".csv");
    
       OFstream os(time_.path()/"postProcessing"/outFile);
       os << "# x, y, z, p, pmean, u(x,y,z), Umean(x,y,z), ";
       os << "UU, VV, WW, UV, VW, UW, Vort(x,y,z)" << endl;

       forAll(cutCells,i)          
       {
           os << mesh.C()[cutCells[i]].x() << ",";
           os << mesh.C()[cutCells[i]].y() << ",";
           os << mesh.C()[cutCells[i]].z() << ",";

           os << p[cutCells[i]] << ",";
           os << pMean[cutCells[i]]<< ",";
           
           os << U[cutCells[i]].x() << ",";
           os << U[cutCells[i]].y() << ",";
           os << U[cutCells[i]].z() << ",";
           
           os << UMean[cutCells[i]].x() << ",";
           os << UMean[cutCells[i]].y() << ",";
           os << UMean[cutCells[i]].z() << ",";
           
           os << UPrime2Mean[cutCells[i]].xx() << ",";
           os << UPrime2Mean[cutCells[i]].yy() << ",";
           os << UPrime2Mean[cutCells[i]].zz() << ",";
           os << UPrime2Mean[cutCells[i]].xy() << ",";
           os << UPrime2Mean[cutCells[i]].yz() << ",";
           os << UPrime2Mean[cutCells[i]].xz() << ",";

           os << vorticity[cutCells[i]].x() << ",";
           os << vorticity[cutCells[i]].y() << ",";
           os << vorticity[cutCells[i]].z();
           os << endl;
       }

    } // if cutCells.size() > 0
    
    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cuttingPlanesFunc::cuttingPlanesFunc
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    point_(dict.lookup("point")),
    dir_(dict.lookup("dir")),
    planeName_(dict.lookup("name")),
    startTime_(readScalar(dict.lookup("startTime")))
{
    Info << "Creating " << this->name() << " function object." << endl;

    // Create postProcessing dir
    fileName outputDir;
    outputDir = time_.path()/"postProcessing";

    if(!isDir(outputDir))
    {
        // Create directory if does not exist.
        mkDir(outputDir);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::cuttingPlanesFunc::start()
{
    return writeData();
}


bool Foam::cuttingPlanesFunc::execute()
{
    if(time_.write() && time_.value()>=startTime_)
    {
        return writeData();
    }
    else
    {
        return true;
    }
}


bool Foam::cuttingPlanesFunc::read(const dictionary& dict)
{
    if(dict.found("startTime"))
    {
        startTime_ = readScalar(dict.lookup("startTime"));
    }

    return true;
}

// ************************************************************************* //
