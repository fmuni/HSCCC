{
    boolList interfaceCells = advector.surf().interfaceCell();

    dimensionedScalar SMALL
    (
        "SMALL",
        dimVol/dimTime,
        1e-32
    );

    alphaPhi = advector.alphaPhi();

    // area of the face immersed in the liquid / total area of the face
    alphaf = alphaPhi/stabilise(phi, SMALL);

    // area of the face immersed in the gas / total area of the face
    alpha2f = 1 - alphaf;

    // Face interpolate density
    surfaceScalarField rhof = alphaf*rho1 + (1-alphaf)*rho2;
    volScalarField alphaOld = alpha1.oldTime();
    volScalarField alpha2Old = 1 - alphaOld;
    surfaceScalarField alpha2Phi = phi.oldTime() - alphaPhi;
    volVectorField interfaceVelocity = (interfaceCenters - interfaceCenters_0) / mesh.time().deltaT();

    // Get the surface interpolation scheme
    IStringStream schemeDataVanLeer("vanLeer01");
    tmp<surfaceInterpolationScheme<scalar>> HOdivSchemeTmp = surfaceInterpolationScheme<scalar>::New(
        mesh,
        phi,
        schemeDataVanLeer
    );

    const surfaceInterpolationScheme<scalar>& HOdivScheme = HOdivSchemeTmp();

    IStringStream schemeDataUpwind("upwind");
    tmp<surfaceInterpolationScheme<scalar>> LOdivSchemeTmp = surfaceInterpolationScheme<scalar>::New(
        mesh,
        phi,
        schemeDataUpwind
    );

    const surfaceInterpolationScheme<scalar>& LOdivScheme = LOdivSchemeTmp();

    forAll(c1Ptr, ci)
    {
        /* =====================================================================
            Interpolate concentrations
        =======================================================================*/

        volScalarField& C1 = c1Ptr[ci];
        volScalarField& C2 = c2Ptr[ci];
        volScalarField& C11 = c11Ptr[ci];
        volScalarField& C22 = c22Ptr[ci];

        Info << "Advecting fields " << C1.name() << ", " << C2.name() << endl;

        const dimensionedScalar D1 = D1Ptr[ci];
        const dimensionedScalar D2 = D2Ptr[ci];
        const dimensionedScalar K = KPtr[ci];
        const dimensionedScalar kl = klPtr[ci];

        //- Interpolation weights
        const surfaceScalarField LOC1Weights = LOdivScheme.weights(C1);
        const surfaceScalarField LOC2Weights = LOdivScheme.weights(C2);
        const surfaceScalarField HOC1Weights = HOdivScheme.weights(C1);
        const surfaceScalarField HOC2Weights = HOdivScheme.weights(C2);        

        //- Interpolate at faces
        surfaceScalarField Cf1 = LOdivScheme.interpolate(C1,LOC1Weights);
        surfaceScalarField Cf2 = LOdivScheme.interpolate(C2,LOC2Weights);
        surfaceScalarField HOCf1 = HOdivScheme.interpolate(C1,HOC1Weights);
        surfaceScalarField HOCf2 = HOdivScheme.interpolate(C2,HOC2Weights);

        //- Compute interpolated values far from interface
        bulkInterpolate(C1, interfaceCells_0, HOCf1, Cf1);
        bulkInterpolate(C2, interfaceCells_0, HOCf2, Cf2);

        //- Compute interpolated values at interfaces 
        interfaceInterpolate(C1, alphaOld, interfaceCells_0, alphaf, alphaPhi, Cf1);
        interfaceInterpolate(C2, alpha2Old, interfaceCells_0, alpha2f, alpha2Phi, Cf2);

        // Parallel exchanges between processors
        syncProcPatches(Cf1, interfaceCells_0);
        syncProcPatches(Cf2, interfaceCells_0);

        /* =====================================================================
            First advection
        =======================================================================*/

        //- Advect
        volScalarField C1Star = C1 - fvc::surfaceIntegrate(phi*Cf1) * runTime.deltaT();
        volScalarField C2Star = C2 - fvc::surfaceIntegrate(phi*Cf2) * runTime.deltaT();

        //- Compute temporary phase averages
        
        forAll(C11, cellI)
        {
            if ( alpha1[cellI]  < 1e-3 )
            {
                C11[cellI] = 0.;
            }
            else
            {
                C11[cellI] = C1Star[cellI] / alpha1[cellI];
            }
        }

        forAll(C22, cellI)
        {
            if ( alpha2[cellI] < 1e-3 )
            {
                C22[cellI] = 0.;
            }
            else
            {
                C22[cellI] = C2Star[cellI] / alpha2[cellI];
            }
        }

        /* =====================================================================
            Flux limiter
        =======================================================================*/
        
        //- The saturation should be limited between 0 an 1
        // TODO: Check if the interface transfer leads to different values
        volScalarField saturation1 = C11;
        volScalarField saturation2 = C22;

        //- Compute the error by interpolating the saturations using the modified upwind schemes
        surfaceScalarField advectionError1 = alphaf;
        phaseAverageInterpolate(saturation1, alphaOld, interfaceCells_0, alphaf, alphaPhi, advectionError1);
        syncProcPatches(advectionError1, interfaceCells_0);

        surfaceScalarField advectionError2 = alpha2f;
        phaseAverageInterpolate(saturation2, alpha2Old, interfaceCells_0, alpha2f, alpha2Phi, advectionError2);
        syncProcPatches(advectionError2, interfaceCells_0); 
        
        surfaceScalarField advectionCorrection1 = advectionError1;
        surfaceScalarField advectionCorrection2 = advectionError2;
 
        //- Compute corrections to advection
        forAll(advectionError1, faceI)
        {
            if (advectionError1[faceI] < 1e-08)
            {
                // This means that C11 < 0, which is not possible. Limit flux to zero.
                advectionCorrection1[faceI] = 0;
            }
            else if (advectionError1[faceI] > alphaf[faceI] )
            {
                // This will always reduce the flux because it means that C11/Csat > 1
                advectionCorrection1[faceI] = alphaf[faceI] / advectionError1[faceI];
            }
            else
            {
                // All other faces
                advectionCorrection1[faceI] = 1;
            }
        }

        forAll(advectionError2, faceI)
        {
            if (advectionError2[faceI] < 1e-08)
            {
                // This means that C11 < 0, which is not possible. Limit flux to zero.
                advectionCorrection2[faceI] = 0;
            }
            else if (advectionError2[faceI] > alpha2f[faceI] )
            {
                // This will always reduce the flux because it means that C11/Csat > 1
                advectionCorrection2[faceI] = alpha2f[faceI] / advectionError2[faceI];
            }
            else
            {
                advectionCorrection2[faceI] = 1;
            }
        }

        //- Update face phase average values
        Cf1 *= advectionCorrection1;
        Cf2 *= advectionCorrection2;

        /* =====================================================================
            Final advection and interface mass transfer
        =======================================================================*/
    
        //- Advect
        C1 = C1 - fvc::surfaceIntegrate(phi*Cf1) * runTime.deltaT();
        C2 = C2 - fvc::surfaceIntegrate(phi*Cf2) * runTime.deltaT();

        //- Recompute phase average values
        forAll(C11, cellI)
        {
            if ( alpha1[cellI]  < 1e-12 )
            {
                C11[cellI] = 0.;
            }
            else
            {
                C11[cellI] = C1[cellI] / alpha1[cellI];
            }
        }

        forAll(C22, cellI)
        {
            if ( alpha2[cellI] < 1e-12 )
            {
                C22[cellI] = 0.;
            }
            else
            {
                C22[cellI] = C2[cellI] / alpha2[cellI];
            }
        }


        //- Use phase averages to recompute values
        //- This can remove some ill behaved cells but makes the method non-conservative
        // C1 = C11 * alpha1;
        // C2 = C22 * alpha2;

        C1.correctBoundaryConditions();
        C2.correctBoundaryConditions();
    }
}
