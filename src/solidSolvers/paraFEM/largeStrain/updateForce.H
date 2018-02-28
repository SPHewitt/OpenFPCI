// Date (1st August 2017)

if (rbfUpdate_)
{
    Info << "Update force using rbf interpolation" << endl;

    // Loop over each patch in solid field D 
    label gi = 0;
    forAll(D_.boundaryField(), patchI)
    {
        // If Patch is defined as "tractionDisplacement" in 0 directory
        if
        (
            isA<tractionDisplacementFvPatchVectorField>
            (
                D_.boundaryField()[patchI]
            )
        )
	{
            tractionDisplacementFvPatchVectorField& patchD =
                refCast<tractionDisplacementFvPatchVectorField>
                (
                    D_.boundaryField()[patchI]
                );
               
//	    const labelList& mp = mesh().boundaryMesh()[patchI].meshPoints();
	    const vectorField& lp = mesh().boundaryMesh()[patchI].localPoints();
	    const vectorField& fc = mesh().boundaryMesh()[patchI].faceCentres();

	    vectorField fTraction = patchD.traction();
	    vectorField pt(lp.size(), vector::zero);

	    const surfaceScalarField& magSf = mesh().magSf();
	    vectorField fForce(patchD.size(),vector::zero);

 	    // Create Face forces from Face tractions and face Areas 
        forAll(patchD,faceI)
        {
            fForce[faceI] = fTraction[faceI];
            fForce[faceI] *= magSf.boundaryField()[patchI][faceI];
        }

	    matrix faceX(fc.size(), 3);
	    matrix faceY(fc.size(), 3);
	    matrix pointX(lp.size(), 3);
	    matrix pointY(lp.size(), 3);

        forAll(fc, faceI)
        {
            faceX(faceI, 0) = fc[faceI].x();
            faceX(faceI, 1) = fc[faceI].y();
            faceX(faceI, 2) = fc[faceI].z();

            faceY(faceI, 0) = fForce[faceI].x();
            faceY(faceI, 1) = fForce[faceI].y();
            faceY(faceI, 2) = fForce[faceI].z();
        }

	    forAll(lp, pointI)
	    {
            pointX(pointI, 0) = lp[pointI].x();
            pointX(pointI, 1) = lp[pointI].y();
            pointX(pointI, 2) = lp[pointI].z();
	    }

	    {
	        std::shared_ptr<RBFFunctionInterface> rbfFunction;
	        rbfFunction = 
                    std::shared_ptr<RBFFunctionInterface>
                    (
                        new TPSFunction()
                    );
	        std::shared_ptr<RBFInterpolation> rbf =
                    std::shared_ptr<RBFInterpolation>
                    (
                       new RBFInterpolation( rbfFunction ) 
                    );

            rbf->compute( faceX, pointX );
	        rbf->interpolate( faceY, pointY );
	    }
 

	    vectorField pTraction(lp.size(), vector::zero);
	    forAll(pTraction, pointI)
	    {
            pTraction[pointI] = 
                vector
                (
                    pointY(pointI, 0),
                    pointY(pointI, 1),
                    pointY(pointI, 2)
                );
	    }

        if(twoDimensional_)
        {
            pTraction /= 2;
        }

        Info << "Total Point Force: " << gSum(pTraction) << endl;
             
            forAll(pTraction, pointI)
            {
                // Populate Force Vector force(FX1 FY1 FZ1.....FXN FYN FZN)
                fext_OF_[gi*ndim+0] = pTraction[pointI].x();
                fext_OF_[gi*ndim+1] = pTraction[pointI].y();
                fext_OF_[gi*ndim+2] = pTraction[pointI].z();
                gi++;

	        }
        }
    }
}
else
{
    Info << "Update force using primitive patch interpolation" << endl;
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
            tractionDisplacementFvPatchVectorField& patchD =
                refCast<tractionDisplacementFvPatchVectorField>
                (
                    D_.boundaryField()[patchI]
                );

 	        const surfaceScalarField& magSf = mesh().magSf();

	        vectorField fTraction = patchD.traction();
	        vectorField fForce(patchD.size(),vector::zero);

 	    // Create Face forces from Face tractions and face Areas 
            forAll(patchD,faceI)
            {
                fForce[faceI] = fTraction[faceI];
		        fForce[faceI] *= magSf.boundaryField()[patchI][faceI];

            }	

	        primitivePatchInterpolation patchInterpolator
            (
                mesh().boundaryMesh()[patchI]
            );	

            vectorField pForce =
            patchInterpolator.faceToPointInterpolate
            (
                fForce
            );

            if(twoDimensional_)
            {
                pForce /= 2;
            }
            
            Info << "Total Point Force: " << gSum(pForce) << endl;

            gi = 0;
            forAll(pForce,pointI)
            {
                    fext_OF_[gi*ndim+0] = pForce[pointI].x();
                    fext_OF_[gi*ndim+1] = pForce[pointI].y();
                    fext_OF_[gi*ndim+2] = pForce[pointI].z();
                    gi++;
            }
	    

        }  // End If (if PAtch = forcePAtch)

    }  // End forAll (forAll Patches)

}  // End Else (rbf vs PrimPatch)
