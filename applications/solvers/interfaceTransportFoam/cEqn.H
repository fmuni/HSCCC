{
    boolList interfaceCells = advector.surf().interfaceCell();
    alphaPhi = advector.alphaPhi();

    surfaceScalarField alpha2Phi = phi.oldTime() - alphaPhi;
    //volVectorField interfaceVelocity = (interfaceCenters - interfaceCenters_0) / mesh.time().deltaT();

    tmp<surfaceInterpolationScheme<scalar>> tdivScheme1
    (
        surfaceInterpolationScheme<scalar>::New
        (
            mesh,
            alphaPhi,
            mesh.divScheme("div(phi,C)")
        )
    );

    tmp<surfaceInterpolationScheme<scalar>> tdivScheme2
    (
        surfaceInterpolationScheme<scalar>::New
        (
            mesh,
            alpha2Phi,
            mesh.divScheme("div(phi,C)")
        )
    );

    surfaceInterpolationScheme<scalar>& divScheme1 = tdivScheme1();
    surfaceInterpolationScheme<scalar>& divScheme2 = tdivScheme2();

    forAll(c1Ptr, ci)
    {

        volScalarField& C1 = c1Ptr[ci];
        volScalarField& C2 = c2Ptr[ci];
        volScalarField& C11 = c11Ptr[ci];
        volScalarField& C22 = c22Ptr[ci];

        dimensionedScalar K = KPtr[ci];
        dimensionedScalar kl = klPtr[ci];

        volScalarField Lambda = mag(advector.surf().normal())*kl;

        /* =====================================================================
            Assemble equations
        =======================================================================*/

        // forAll(C11, cellI)
        // {
        //     if ( alpha1[cellI]  < 1e-3 )
        //     {
        //         C11[cellI] = 0.;
        //     }
        //     else
        //     {
        //         C11[cellI] = C1[cellI] / alpha1[cellI];
        //     }
        // }

        // forAll(C22, cellI)
        // {
        //     if ( alpha2[cellI] < 1e-3 )
        //     {
        //         C22[cellI] = 0.;
        //     }
        //     else
        //     {
        //         C22[cellI] = C2[cellI] / alpha2[cellI];
        //     }
        // }

        // C11.correctBoundaryConditions();
        // C22.correctBoundaryConditions();

        fvScalarMatrix C1Eqn
        (
            fvm::ddt(C1)
        //   + fvc::div(alphaPhi, C11,"div(phi,C)")
        );


        fvScalarMatrix C2Eqn
        (
            fvm::ddt(C2)
        //   + fvc::div(alpha2Phi, C22,"div(phi,C)")
        );




        /* =====================================================================
            Add interphase transfer -- partial elimination
        =======================================================================*/
        
        scalarField& a1 = C1Eqn.diag();
        scalarField& s1 = C1Eqn.source();
        scalarField& a2 = C2Eqn.diag();
        scalarField& s2 = C2Eqn.source();

        const volScalarField H1 = C1Eqn.H();
        const volScalarField H2 = C2Eqn.H();

        forAll(interfaceCells, Icell)
        {
            label cellI = interfaceCells[Icell];

            scalar a1o = a1[cellI];

            scalar L = Lambda[cellI];
            scalar KL = L * K.value();

            a1[cellI] += L / ( a2[cellI] + KL ) * a2[cellI];
            s1[cellI] += H2[cellI] * KL / ( a2[cellI] + KL  );

            a2[cellI] += KL / ( a1o + L ) * a1o;
            s2[cellI] += H1[cellI] * L / ( a1o + L  );
        }

        C1Eqn.solve();
        C2Eqn.solve();

    }
}
