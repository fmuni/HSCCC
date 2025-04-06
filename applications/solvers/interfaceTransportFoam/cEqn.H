{
    //boolList interfaceCells = advector.surf().interfaceCell();
    alphaPhi = advector.alphaPhi();

    surfaceScalarField alpha2Phi = phi.oldTime() - alphaPhi;
    //volVectorField interfaceVelocity = (interfaceCenters - interfaceCenters_0) / mesh.time().deltaT();

    forAll(c1Ptr, ci)
    {
        /* =====================================================================
            Interpolate concentrations
        =======================================================================*/

        volScalarField& C1 = c1Ptr[ci];
        volScalarField& C2 = c2Ptr[ci];
        volScalarField& C11 = c11Ptr[ci];
        volScalarField& C22 = c22Ptr[ci];

        forAll(C11, cellI)
        {
            if ( alpha1[cellI]  < 1e-3 )
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
            if ( alpha2[cellI] < 1e-3 )
            {
                C22[cellI] = 0.;
            }
            else
            {
                C22[cellI] = C2[cellI] / alpha2[cellI];
            }
        }

        C11.correctBoundaryConditions();
        C22.correctBoundaryConditions();

        fvScalarMatrix C1Eqn
        (
            fvm::ddt(C1)
          + fvc::div(alphaPhi, C11,"div(phi,C)")
        );

        C1Eqn.solve();

        fvScalarMatrix C2Eqn
        (
            fvm::ddt(C2)
          + fvc::div(alpha2Phi, C22,"div(phi,C)")
        );

        C2Eqn.solve();

    }
}
