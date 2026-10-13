if (specieNames.size() > 0)
{
    const scalar alphaTolerance = 1e-4;
    boolList interfaceCells = advector.surf().interfaceCell();
    alphaPhi = advector.alphaPhi();

    surfaceScalarField alpha2Phi = phi - alphaPhi;
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

    const surfaceInterpolationScheme<scalar>& divScheme1 = tdivScheme1();
    const surfaceInterpolationScheme<scalar>& divScheme2 = tdivScheme2();


    forAll(c1Ptr, ci)
    {

        volScalarField& C1 = c1Ptr[ci];
        volScalarField& C2 = c2Ptr[ci];
        volScalarField& C11 = c11Ptr[ci];
        volScalarField& C22 = c22Ptr[ci];

        dimensionedScalar K = KPtr[ci];
        dimensionedScalar kl = klPtr[ci];

        dimensionedScalar D1 = D1Ptr[ci];
        dimensionedScalar D2 = D2Ptr[ci];

        volScalarField Lambda = mag(advector.surf().normal())*kl;

        C1.correctBoundaryConditions();
        C2.correctBoundaryConditions();

        tmp<fvScalarMatrix> tC1Adv
        (
            new fvScalarMatrix
            (
                C1,
                C1.dimensions()/dimTime*dimVol
            )
        );

        fvScalarMatrix& C1Adv = tC1Adv.ref();

        tmp<fvScalarMatrix> tC2Adv
        (
            new fvScalarMatrix
            (
                C2,
                C2.dimensions()/dimTime*dimVol
            )
        );

        fvScalarMatrix& C2Adv= tC2Adv.ref();

        /* =====================================================================
            Advection operator for C1
        =======================================================================*/

        tmp<surfaceScalarField> tweights1 = tdivScheme1().weights(C1);
        const surfaceScalarField& weights1 = tweights1();


        forAll(C1Adv.lower(),ownerI)
        {
            scalar coeff = alpha1[mesh.owner()[ownerI]];

            if( coeff < alphaTolerance )
            {
                coeff = 0;
            }
            else
            {
                coeff = 1.0/coeff;
            }

            C1Adv.lower()[ownerI] = -coeff*weights1.primitiveField()[ownerI]*alphaPhi.primitiveField()[ownerI];

        }

        forAll(C1Adv.upper(),neighI)
        {
            scalar coeff = alpha1[mesh.neighbour()[neighI]];

            if( coeff < alphaTolerance )
            {
                coeff = 0;
            }
            else
            {
                coeff = 1.0/coeff;
            }

            C1Adv.upper()[neighI] =  coeff* ( 1.0 - weights1.primitiveField()[neighI])*alphaPhi.primitiveField()[neighI];

        }

        C1Adv.negSumDiag();

        // Compute phase average

        forAll(C11, cellI)
        {
            if ( alpha1[cellI]  < alphaTolerance )
            {
                C11[cellI] = 0.;
            }
            else
            {
                C11[cellI] = C1[cellI] / alpha1[cellI];
            }
        }

        C11.correctBoundaryConditions();
        forAll(C1.boundaryField(), patchi)
        {
            fvPatchField<scalar>& psf = C11.boundaryFieldRef()[patchi];
            const fvPatchField<scalar>& C1f = C1.boundaryField()[patchi];
            const fvPatchField<scalar>& a1f = alpha1.boundaryField()[patchi];

            forAll(psf,faceI)
            {
                if ( a1f[faceI]  < alphaTolerance )
                {
                    psf[faceI] = 0.;
                }
                else
                {
                    psf[faceI] = C1f[faceI] / a1f[faceI];
                }                
            }

        }

        forAll(C1.boundaryField(), patchi)
        {
            const fvPatchField<scalar>& psf = C1.boundaryField()[patchi];

            const fvsPatchScalarField& patchFlux = alphaPhi.boundaryField()[patchi];
            const fvsPatchScalarField& pw = weights1.boundaryField()[patchi];
            const fvPatchField<scalar>& a1f = alpha1.boundaryField()[patchi];
            const auto& icoeffs = psf.valueInternalCoeffs(pw);
            const auto& bcoeffs = psf.valueBoundaryCoeffs(pw);

            forAll(psf,faceI)
            {
                if ( a1f[faceI]  < alphaTolerance )
                {
                    C1Adv.internalCoeffs()[patchi][faceI] += 0.;
                    C1Adv.boundaryCoeffs()[patchi][faceI] += 0.;
                }
                else
                {
                    C1Adv.internalCoeffs()[patchi][faceI] += patchFlux[faceI]*icoeffs()[faceI]/a1f[faceI];
                    C1Adv.boundaryCoeffs()[patchi][faceI] -= patchFlux[faceI]*bcoeffs()[faceI]/a1f[faceI];
                }                
            }
        }

        if (divScheme1.corrected())
        {
            C1Adv += fvc::surfaceIntegrate(alphaPhi*divScheme1.correction(C11));
        }

        /* =====================================================================
            Advection operator for C2
        =======================================================================*/

        tmp<surfaceScalarField> tweights2 = tdivScheme2().weights(C2);
        const surfaceScalarField& weights2 = tweights2();

        forAll(C2Adv.lower(),ownerI)
        {
            scalar coeff = alpha2[mesh.owner()[ownerI]];

            if( coeff < alphaTolerance )
            {
                coeff = 0;
            }
            else
            {
                coeff = 1.0/coeff;
            }

            C2Adv.lower()[ownerI] = -coeff*weights2.primitiveField()[ownerI]*alpha2Phi.primitiveField()[ownerI];

        }

        forAll(C2Adv.upper(),neighI)
        {
            scalar coeff = alpha2[mesh.neighbour()[neighI]];

            if( coeff < alphaTolerance )
            {
                coeff = 0;
            }
            else
            {
                coeff = 1.0/coeff;
            }

            C2Adv.upper()[neighI] =  coeff* ( 1.0 - weights2.primitiveField()[neighI])*alpha2Phi.primitiveField()[neighI];

        }

        C2Adv.negSumDiag();

        // Compute phase average

        forAll(C22, cellI)
        {
            if ( alpha2[cellI]  < alphaTolerance )
            {
                C22[cellI] = 0.;
            }
            else
            {
                C22[cellI] = C2[cellI] / alpha2[cellI];
            }
        }

        C22.correctBoundaryConditions();
        forAll(C2.boundaryField(), patchi)
        {
            fvPatchField<scalar>& psf = C22.boundaryFieldRef()[patchi];
            const fvPatchField<scalar>& C2f = C2.boundaryField()[patchi];
            const fvsPatchField<scalar>& a2f = alpha2f.boundaryField()[patchi];

            forAll(psf,faceI)
            {
                if ( a2f[faceI]  < alphaTolerance )
                {
                    psf[faceI] = 0.;
                }
                else
                {
                    psf[faceI] = C2f[faceI] / a2f[faceI];
                }                
            }

        }

        forAll(C2.boundaryField(), patchi)
        {
            const fvPatchField<scalar>& psf = C2.boundaryField()[patchi];

            const fvsPatchScalarField& patchFlux = alpha2Phi.boundaryField()[patchi];
            const fvsPatchScalarField& pw = weights2.boundaryField()[patchi];
            const fvPatchField<scalar>& a2f = alpha2.boundaryField()[patchi];
            const auto& icoeffs = psf.valueInternalCoeffs(pw);
            const auto& bcoeffs = psf.valueBoundaryCoeffs(pw);

            forAll(psf,faceI)
            {
                if ( a2f[faceI]  < alphaTolerance )
                {
                    C2Adv.internalCoeffs()[patchi][faceI] += 0.;
                    C2Adv.boundaryCoeffs()[patchi][faceI] += 0.;
                }
                else
                {
                    C2Adv.internalCoeffs()[patchi][faceI] += patchFlux[faceI]*icoeffs()[faceI]/a2f[faceI];
                    C2Adv.boundaryCoeffs()[patchi][faceI] -= patchFlux[faceI]*bcoeffs()[faceI]/a2f[faceI];
                }                
            }
        }

        if (divScheme2.corrected())
        {
            C2Adv += fvc::surfaceIntegrate(alpha2Phi*divScheme2.correction(C22));
        }
        
        /* =====================================================================
            Assemble Matrices
        =======================================================================*/

        surfaceScalarField D1f = mixture.nuf();
        surfaceScalarField D2f = D1f;

        forAll(mesh.neighbour(), faceI)
        {
            label cellO = mesh.owner()[faceI];
            label cellN = mesh.neighbour()[faceI];
            
            if (alpha1[cellO] < alphaTolerance || alpha1[cellN] < alphaTolerance)
            {
                D1f[faceI] = 0.;
            }
            else
            {
                D1f[faceI] = D1.value();
            }

            if (alpha2[cellO] < alphaTolerance || alpha2[cellN] < alphaTolerance)
            {
                D2f[faceI] = 0.;
            }
            else
            {
                D2f[faceI] = D2.value();
            }
        }

        fvScalarMatrix C1Eqn
        (
            fvm::ddt(C1) 
        //  + C1Adv 
            + fvc::div(alphaPhi, C11, "div(alphaPhi,C)")
        //  - fvm::laplacian(D1f, C1)
        );
        
        fvScalarMatrix C2Eqn
        (
            fvm::ddt(C2) 
        //  + C2Adv
            + fvc::div(alpha2Phi, C22, "div(alphaPhi,C)") 
        //  - fvm::laplacian(D2f, C2)
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

        forAll(Lambda, cellI)
        {
            if(alpha1[cellI] > 1e-2 && alpha2[cellI] > 1e-2)
            {
                a1[cellI] += Lambda[cellI];// / alpha1[cellI];
                s1[cellI] += Lambda[cellI]*K.value()*C2[cellI]; 
                a2[cellI] += K.value()*Lambda[cellI];// / alpha2[cellI];
                s2[cellI] += Lambda[cellI]*C1[cellI];  

            }
                // scalar a1o = a1[cellI];

            // scalar L = Lambda[cellI];
            // scalar KL = L * K.value();

            // a1[cellI] += L / ( a2[cellI] + KL ) * a2[cellI];
            // s1[cellI] += H2[cellI] * KL / ( a2[cellI] + KL  );

            // a2[cellI] += KL / ( a1o + L ) * a1o;
            // s2[cellI] += H1[cellI] * L / ( a1o + L  );
        }

        /* =====================================================================
            Solve the equations
        =======================================================================*/

        // C1Eqn.relax();
        // C2Eqn.relax();

        C1Eqn.solve();
        C2Eqn.solve();

    }
}
// {
//     //boolList interfaceCells = advector.surf().interfaceCell();
//     alphaPhi = advector.alphaPhi();

//     surfaceScalarField alpha2Phi = phi.oldTime() - alphaPhi;
//     //volVectorField interfaceVelocity = (interfaceCenters - interfaceCenters_0) / mesh.time().deltaT();

//     forAll(c1Ptr, ci)
//     {
//         /* =====================================================================
//             Interpolate concentrations
//         =======================================================================*/

//         volScalarField& C1 = c1Ptr[ci];
//         volScalarField& C2 = c2Ptr[ci];
//         volScalarField& C11 = c11Ptr[ci];
//         volScalarField& C22 = c22Ptr[ci];

//         forAll(C11, cellI)
//         {
//             if ( alpha1[cellI]  < alphaTolerance )
//             {
//                 C11[cellI] = 0.;
//             }
//             else
//             {
//                 C11[cellI] = C1[cellI] / alpha1[cellI];
//             }
//         }

//         forAll(C22, cellI)
//         {
//             if ( alpha2[cellI] < alphaTolerance )
//             {
//                 C22[cellI] = 0.;
//             }
//             else
//             {
//                 C22[cellI] = C2[cellI] / alpha2[cellI];
//             }
//         }

//         C11.correctBoundaryConditions();
//         C22.correctBoundaryConditions();

//         fvScalarMatrix C1Eqn
//         (
//             fvm::ddt(C1)
//           + fvc::div(alphaPhi, C11,"div(phi,C)")
//         );

//         C1Eqn.solve();

//         fvScalarMatrix C2Eqn
//         (
//             fvm::ddt(C2)
//           + fvc::div(alpha2Phi, C22,"div(phi,C)")
//         );

//         C2Eqn.solve();

//     }
// }