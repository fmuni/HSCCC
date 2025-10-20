    MRF.correctBoundaryVelocity(U);

    #include "computePlanetaryForce.H"

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(rho, U)
     ==
        fvOptions(rho, U)
    );

    if(!facePlanetaryForce)
    {
        UEqn += F_planetary;
    }

    #include "UEqnAddPorosity.H"

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        if(!facePlanetaryForce)
        {
            solve
            (
                UEqn
             ==
                fvc::reconstruct
                (
                    (
                        mixture.surfaceTensionForce()
                      - ghf*fvc::snGrad(rho)
                      - fvc::snGrad(p_rgh)
                    ) * mesh.magSf()
                )
            );
        }
        else
        {
            solve
            (
                UEqn
             ==
                fvc::reconstruct
                (
                    (
                        mixture.surfaceTensionForce()
                      - ghf*fvc::snGrad(rho)
                      - fvc::snGrad(p_rgh)
                    ) * mesh.magSf()
                    + phiFP
                )
            );
        }
        fvOptions.correct(U);
    }
