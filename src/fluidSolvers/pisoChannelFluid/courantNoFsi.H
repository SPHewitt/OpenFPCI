    // CourantNo
    {
        scalar CoNum = 0.0;
        scalar meanCoNum = 0.0;
        scalar velMag = 0.0;

        if (mesh.nInternalFaces())
        {
            surfaceScalarField SfUfbyDelta =
                mesh.surfaceInterpolation::deltaCoeffs()*mag(phi_);

            CoNum = max(SfUfbyDelta/mesh.magSf())
                .value()*runTime().deltaT().value();

            meanCoNum = (sum(SfUfbyDelta)/sum(mesh.magSf()))
                .value()*runTime().deltaT().value();

            velMag = max(mag(phi_)/mesh.magSf()).value();
        }

        Info<< "Courant Number mean: " << meanCoNum
            << " max: " << CoNum
            << " velocity magnitude: " << velMag << endl;

        // Adjustable Time Step    
        //#include "setDeltaTFsi.H"
        //Info<< "Time = " << runTime().timeName() << nl << endl;
    }
