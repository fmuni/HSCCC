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
    alphaf = alphaPhi/stabilise(phi.oldTime(), SMALL);

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
        volScalarField& C1 = c1Ptr[ci];
        volScalarField& C2 = c1Ptr[ci];

        const dimensionedScalar D1 = D1Ptr[ci];
        const dimensionedScalar D2 = D2Ptr[ci];
        const dimensionedScalar K = KPtr[ci];
        const dimensionedScalar kl = klPtr[ci];

        //- Concentrations averaged over volumes
        volScalarField c1_ave = C1 / ( alpha1.oldTime() + 1e-12);
        volScalarField c2_ave = C2 / ( alpha2.oldTime() + 1e-12);

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
    }
}
