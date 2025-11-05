    MRF.correctBoundaryVelocity(U);

    #include "computePlanetaryForce.H"

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(rho, U)
     ==
        (
            (facePlanetaryForce) ?
            fvOptions(rho,U)
            :
            F_planetary + fvOptions(rho,U)
        )
    );

    #include "UEqnAddPorosity.H"

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        surfaceScalarField UForceFace
        (
            (
                mixture.surfaceTensionForce()
              - ghf*fvc::snGrad(rho)
              - fvc::snGrad(p_rgh)
            ) * mesh.magSf()

        );

        if (facePlanetaryForce)
        {
            UForceFace += phiFP;
        }

        solve (UEqn == fvc::reconstruct(UForceFace));

        fvOptions.correct(U);
    }
