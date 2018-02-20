{
    rUA_ = 1.0/UEqn.A();

    U_ = rUA_*UEqn.H();

    phi_ = (fvc::interpolate(U_) & mesh.Sf());
         //+ (fvc::ddtPhiCorr(rUA_, U_, phi_));

    //adjustPhi(phi_, U_, p_);	

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rUA_, p_) == fvc::div(phi_)
        );

        pEqn.setReference(pRefCell, pRefValue);

        if
        (
            corr == nCorr-1
         && nonOrth == nNonOrthCorr
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

    #include "checkContinuityErrors.H"

    // Explicitly relax pressure for momentum corrector
    if (oCorr != nOuterCorr - 1)
    {
        p_.relax();
    }

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi_, U_);

    gradp_ = fvc::grad(p_);

    U_ -= rUA_*gradp_;
    U_.correctBoundaryConditions();

    gradU_ = fvc::grad(U_);
}
