/**
 * Computes the planetary force applied to the NS equations
 */

if ( mag(OmegaZero).value() > 1e-12 )
{

    // Update omega
    dimensionedVector Omega = OmegaZero;

    scalar t = runTime.value();

    if ( t < rumpTime.value() )
    {
        Omega *= t/rumpTime.value();
    }

    vector rotDit = (OmegaZero / mag(OmegaZero)).value();

    scalar R = Rs.value();
    scalar omega = mag(Omega).value() * constant::mathematical::twoPi / 60.0;

    forAll(F_planetary, cellI)
    {
        vector Xr = mesh.C()[cellI] - Xp.value();
        scalar r = mag(Xr&vector(1,1,0));

        scalar ang_phi = std::atan2(Xr[1],Xr[0]);

        scalar beta = r/R;

        // Compute global acceleration
        scalar axR = -R*sqr(omega)*(Foam::cos(omega*t) + 4.0*beta*Foam::cos(2.0*omega*t + ang_phi)  );
        scalar ayR = -R*sqr(omega)*(Foam::sin(omega*t) + 4.0*beta*Foam::sin(2.0*omega*t + ang_phi)  );

        // Compute local acceleration
        scalar axr = axR*Foam::cos(omega*2.0*t) + ayR*Foam::sin(2.0*omega*t);
        scalar ayr = -axR*Foam::sin(omega*2.0*t) + ayR*Foam::cos(2.0*omega*t);

        F_planetary[cellI][0] = axr*rho[cellI];
        F_planetary[cellI][1] = ayr*rho[cellI];
        F_planetary[cellI][2] = 0.;
    }

    if( facePlanetaryForce )
    {
        surfaceScalarField rhof = fvc::interpolate(rho);

        forAll(phiFP, faceI)
        {
            vector Xr = mesh.Cf()[faceI] - Xp.value();
            scalar r = mag(Xr&vector(1,1,0));

            scalar ang_phi = std::atan2(Xr[1],Xr[0]);

            scalar beta = r/R;

            // Compute global acceleration
            scalar axR = -R*sqr(omega)*(Foam::cos(omega*t) + 4.0*beta*Foam::cos(2.0*omega*t + ang_phi)  );
            scalar ayR = -R*sqr(omega)*(Foam::sin(omega*t) + 4.0*beta*Foam::sin(2.0*omega*t + ang_phi)  );

            // Compute local acceleration
            scalar axr = axR*Foam::cos(omega*2.0*t) + ayR*Foam::sin(2.0*omega*t);
            scalar ayr = -axR*Foam::sin(omega*2.0*t) + ayR*Foam::cos(2.0*omega*t);

            phiFP[faceI] = rhof[faceI]*(axr*mesh.Sf()[faceI][0] + ayr*mesh.Sf()[faceI][1]);
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
            scalar r = mag(Xr&vector(1,1,0));

            scalar ang_phi = std::atan2(Xr[1],Xr[0]);

            scalar beta = r/R;

            // Compute global acceleration
            scalar axR = -R*sqr(omega)*(Foam::cos(omega*t) + 4.0*beta*Foam::cos(2.0*omega*t + ang_phi)  );
            scalar ayR = -R*sqr(omega)*(Foam::sin(omega*t) + 4.0*beta*Foam::sin(2.0*omega*t + ang_phi)  );

            // Compute local acceleration
            scalar axr = axR*Foam::cos(omega*2.0*t) + ayR*Foam::sin(2.0*omega*t);
            scalar ayr = -axR*Foam::sin(omega*2.0*t) + ayR*Foam::cos(2.0*omega*t);

            Fps[faceI][0] = axr*rhof[faceI];
            Fps[faceI][1] = ayr*rhof[faceI];
            Fps[faceI][2] = 0.;

            if(facePlanetaryForce)
            {
                phiPFs[faceI] = Fps[faceI] & Fps.patch().Sf()[faceI];
            }
        }
    }
}
