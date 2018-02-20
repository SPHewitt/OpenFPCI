fvVectorMatrix UEqn
(
    fvm::ddt(U_)
  + fvm::div(phi_, U_)
  + turbulence_->divDevReff()
);

/*
    if (oCorr == nOuterCorr - 1)
    {
        if (mesh.solutionDict().relax("UFinal"))
        {
            UEqn.relax(mesh.solutionDict().relaxationFactor("UFinal"));
        }
        else
        {
            UEqn.relax(1);
        }
    }
    else
    {
        UEqn.relax();
    }
*/

if (momentumPredictor)
{
    solve(UEqn == -fvc::grad(p_));
}

