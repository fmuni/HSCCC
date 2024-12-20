forAll(cPtr, ci)
{
    volScalarField& c = cPtr[ci];
    const dimensionedScalar D1 = D1Ptr[ci];
    const dimensionedScalar D2 = D2Ptr[ci];
    const dimensionedScalar H = HPtr[ci];

    surfaceScalarField alpha1f = fvc::interpolate(alpha1);
    surfaceScalarField alpha2f = 1. - alpha1f;

    //- Harmonic average of the diffusion coefficients
    surfaceScalarField Dh( D1*D2 / ( alpha1f*D1 + alpha2f*D2 ) );

    //- Drift flux
    surfaceScalarField phiH
    (
        Dh * ( 1. - 1./H) / ( alpha1f - alpha2f/H ) * fvc::snGrad(alpha1) * mesh.magSf()
    );

    fvScalarMatrix cEqn
    (
        fvm::ddt(c)
      + fvm::div(phi,c)
      + fvm::div(phiH,c,"div(phiH,c)")
      - fvm::laplacian(Dh,c,"laplacian(D,c)")
    );

    cEqn.relax();
    cEqn.solve();
    c.correctBoundaryConditions();
}