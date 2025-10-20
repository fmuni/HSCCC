/**
 * Computes the planetary force applied to the NS equations
 */

scalar t = runTime.value();
scalar omega = omegaPtr->value(t) * constant::mathematical::twoPi / 60.0;

if ( mag(omega) > 1e-12 )
{
    vector rotDit(0,0,1);

    scalar R = Rs.value();

    forAll(F_planetary, cellI)
    {
        vector Xr = mesh.C()[cellI] - Xp.value();

        F_planetary[cellI] = computePlanetaryAcceleration(Xr.x(), Xr.y(), R, omega, t)*rho[cellI];
    }

    if( facePlanetaryForce )
    {
        surfaceScalarField rhof = fvc::interpolate(rho);

        forAll(phiFP, faceI)
        {
            vector Xr = mesh.Cf()[faceI] - Xp.value();

            vector acc = computePlanetaryAcceleration(Xr.x(), Xr.y(), R, omega, t);

            phiFP[faceI] = rhof[faceI]*(acc & mesh.Sf()[faceI]);
        }
    }

    // Update patches as well
    forAll(F_planetary.boundaryField(),patchI)
    {
        fvPatchField<vector>& Fps = F_planetary.boundaryFieldRef()[patchI];
        const fvPatchField<vector>& Cf = mesh.C().boundaryField()[patchI];
        const fvPatchField<scalar>& rhof = rho.boundaryField()[patchI];
        fvsPatchField<scalar>& phiPFs = phiFP.boundaryFieldRef()[patchI];

        forAll(Fps,faceI)
        {
            vector Xr =Cf[faceI] - Xp.value();
            vector acc = computePlanetaryAcceleration(Xr.x(), Xr.y(), R, omega, t);

            Fps[faceI] = acc*rhof[faceI];

            if(facePlanetaryForce)
            {
                phiPFs[faceI] = Fps[faceI] & Fps.patch().Sf()[faceI];
            }
        }
    }
}
