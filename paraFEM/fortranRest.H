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
    fortranRest.H

Description
    Class formats the OpenFOAM boundary information into  format read
    by ParaFEM 

Author
    Sam Hewitt.  All rights reserved.

SourceFiles
    fortranRest.C

\*---------------------------------------------------------------------------*/

#ifndef fortranRest_H
#define fortranRest_H

//#include "solidSolver.H"
//#include "volFields.H"
//#include "surfaceFields.H"
//#include "pointFields.H"
//#include "constitutiveModel.H"
//#include "leastSquaresVolPointInterpolation.H"
//#include "pointVolInterpolation.H"
//#include "TLMaterialInterface.H"
//#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace std
{

/*---------------------------------------------------------------------------*\
                        Class fortranRest Declaration
\*---------------------------------------------------------------------------*/

class fortranRest
{
    // Private data
  
        //- List of Node Numbers
        int* Node_;

        //- X restraint
        int* x_;

        //- Y restraint
        int* y_;

        //- Z restraint
        int* z_;

        //- number of restrained nodes
        int size_;

        //- Counter
        int  counter_; 

    // Private Member Functions

        //- Disallow default bitwise copy construct
        fortranRest(const fortranRest&);

        //- Disallow default bitwise assignment
        void operator=(const fortranRest&);


public:

    // Constructors

        //- Construct from components
        fortranRest(int size);


    // Destructor
        ~fortranRest();

    // Member Functions

        // Access

            //- Add a Node and its restraints to object 
            void addNode(int i, int x, int y,int z);

            //- Edit the Node 
	    void editNode(int i, int x, int y, int z);   
	
	    //- Print the object to fortran array
            void print(int *array);
	    
	    //- Print object in ENSI Gold Format
	    void printENSI(int* array);

	    //- Print Object from ENSI 
	    void printFromENSI(int* array,int arraySize);

	    //- get Counter Value
	    int getCounter();

	    //- Get Num of rest Nodes
	    int getNumRestNodes();

};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace std

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
