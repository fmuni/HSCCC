    MRF.correctBoundaryVelocity(U);


    //Implicit discretization of the planetary force
    // volScalarField magFcoeff = mag(F_planetary)/dimensionedScalar("magU",dimVelocity,1.0);
    //
    // forAll(magFcoeff, cellI)
    // {
    //     if ( mag(U[cellI]) > 1e-5 )
    //     {
    //         magFcoeff[cellI] /= mag(U[cellI]);
    //     }
    //     else
    //     {
    //        magFcoeff[cellI] = 0.;
    //     }
    // }

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(rho, U)
//      + fvm::Sp(magFcoeff,U)
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
            UForceFace += phiFPrelax;
        }

        solve (UEqn == fvc::reconstruct(UForceFace));

        fvOptions.correct(U);
    }
