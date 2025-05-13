if (specieNames.size() > 0)
{
    boolList interfaceCells = advector.surf().interfaceCell();
    alphaPhi = advector.alphaPhi();

    surfaceScalarField alpha2Phi = phi.oldTime() - alphaPhi;
    //volVectorField interfaceVelocity = (interfaceCenters - interfaceCenters_0) / mesh.time().deltaT();

    // tmp<surfaceInterpolationScheme<scalar>> tdivScheme1
    // (
    //     surfaceInterpolationScheme<scalar>::New
    //     (
    //         mesh,
    //         alphaPhi,
    //         mesh.divScheme("div(phi,C)")
    //     )
    // );

    // tmp<surfaceInterpolationScheme<scalar>> tdivScheme2
    // (
    //     surfaceInterpolationScheme<scalar>::New
    //     (
    //         mesh,
    //         alpha2Phi,
    //         mesh.divScheme("div(phi,C)")
    //     )
    // );

    // const surfaceInterpolationScheme<scalar>& divScheme1 = tdivScheme1();
    // const surfaceInterpolationScheme<scalar>& divScheme2 = tdivScheme2();


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
            Material derivative for C1
        =======================================================================*/

        // tmp<surfaceScalarField> tweights1 = tdivScheme1().weights(C1);
        // const surfaceScalarField& weights1 = tweights1();

        // fvScalarMatrix C1Eqn( fvm::ddt(C1) );

        // forAll(C1Eqn.lower(),ownerI)
        // {
        //     scalar coeff = alpha1[mesh.owner()[ownerI]];

        //     if( coeff < 1e-3 )
        //     {
        //         coeff = 0;
        //     }
        //     else
        //     {
        //         coeff = 1.0/coeff;
        //     }

        //     C1Eqn.lower()[ownerI] = -coeff*weights1.primitiveField()[ownerI]*alphaPhi.primitiveField()[ownerI];

        // }

        // forAll(C1Eqn.upper(),neighI)
        // {
        //     scalar coeff = alpha1[mesh.neighbour()[neighI]];

        //     if( coeff < 1e-3 )
        //     {
        //         coeff = 0;
        //     }
        //     else
        //     {
        //         coeff = 1.0/coeff;
        //     }

        //     C1Eqn.upper()[neighI] =  coeff* ( 1.0 - weights1.primitiveField()[neighI])*alphaPhi.primitiveField()[neighI];

        // }

        // C1Eqn.negSumDiag();

        // // Compute phase average

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

        // C11.correctBoundaryConditions();
        // // forAll(C1.boundaryField(), patchi)
        // // {
        // //     fvPatchField<scalar>& psf = C11.boundaryFieldRef()[patchi];
        // //     const fvPatchField<scalar>& C1f = C1.boundaryField()[patchi];
        // //     const fvPatchField<scalar>& a1f = alpha1.boundaryField()[patchi];

        // //     forAll(psf,faceI)
        // //     {
        // //         if ( alphaf[faceI]  < 1e-4 )
        // //         {
        // //             psf[faceI] = 0.;
        // //         }
        // //         else
        // //         {
        // //             psf[faceI] = C1f[faceI] / a1f[faceI];
        // //         }                
        // //     }

        // // }

        // forAll(C1.boundaryField(), patchi)
        // {
        //     const fvPatchField<scalar>& psf = C11.boundaryField()[patchi];
        //     const fvsPatchScalarField& patchFlux = alphaPhi.boundaryField()[patchi];
        //     const fvsPatchScalarField& pw = weights1.boundaryField()[patchi];

        //     C1Eqn.internalCoeffs()[patchi] = patchFlux*psf.valueInternalCoeffs(pw);
        //     C1Eqn.boundaryCoeffs()[patchi] = -patchFlux*psf.valueBoundaryCoeffs(pw);
        // }

        // if (divScheme1.corrected())
        // {
        //     C1Eqn += fvc::surfaceIntegrate(alphaPhi*divScheme1.correction(C11));
        // }

        // /* =====================================================================
        //     Material derivative for C2
        // =======================================================================*/

        // tmp<surfaceScalarField> tweights2 = tdivScheme2().weights(C2);
        // const surfaceScalarField& weights2 = tweights2();

        // fvScalarMatrix C2Eqn( fvm::ddt(C2) );

        // forAll(C2Eqn.lower(),ownerI)
        // {
        //     scalar coeff = alpha2[mesh.owner()[ownerI]];

        //     if( coeff < 1e-3 )
        //     {
        //         coeff = 0;
        //     }
        //     else
        //     {
        //         coeff = 1.0/coeff;
        //     }

        //     C2Eqn.lower()[ownerI] = -coeff*weights2.primitiveField()[ownerI]*alpha2Phi.primitiveField()[ownerI];

        // }

        // forAll(C2Eqn.upper(),neighI)
        // {
        //     scalar coeff = alpha2[mesh.neighbour()[neighI]];

        //     if( coeff < 1e-3 )
        //     {
        //         coeff = 0;
        //     }
        //     else
        //     {
        //         coeff = 1.0/coeff;
        //     }

        //     C2Eqn.upper()[neighI] =  coeff* ( 1.0 - weights2.primitiveField()[neighI])*alpha2Phi.primitiveField()[neighI];

        // }

        // C2Eqn.negSumDiag();

        // // Compute phase average

        // forAll(C22, cellI)
        // {
        //     if ( alpha2[cellI]  < 1e-4 )
        //     {
        //         C22[cellI] = 0.;
        //     }
        //     else
        //     {
        //         C22[cellI] = C2[cellI] / alpha2[cellI];
        //     }
        // }

        // C22.correctBoundaryConditions();

        // forAll(C2.boundaryField(), patchi)
        // {
        //     const fvPatchField<scalar>& psf = C22.boundaryField()[patchi];
        //     const fvsPatchScalarField& patchFlux = alpha2Phi.boundaryField()[patchi];
        //     const fvsPatchScalarField& pw = weights2.boundaryField()[patchi];

        //     C2Eqn.internalCoeffs()[patchi] = patchFlux*psf.valueInternalCoeffs(pw);
        //     C2Eqn.boundaryCoeffs()[patchi] = -patchFlux*psf.valueBoundaryCoeffs(pw);
        // }

        // if (divScheme2.corrected())
        // {
        //     C2Eqn += fvc::surfaceIntegrate(alpha2Phi*divScheme2.correction(C22));
        // }


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


        fvScalarMatrix C2Eqn
        (
            fvm::ddt(C2)
          + fvc::div(alpha2Phi, C22,"div(phi,C)")
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

            if(alpha1[cellI] > 1e-3 && alpha2[cellI] > 1e-3)
            {
                a1[cellI] += Lambda[cellI] / alpha1[cellI];
                s1[cellI] += K.value()*Lambda[cellI]*C22[cellI]; 
                a2[cellI] += K.value()*Lambda[cellI] / alpha2[cellI];
                s2[cellI] += Lambda[cellI]*C11[cellI]; 

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
//             if ( alpha1[cellI]  < 1e-3 )
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
//             if ( alpha2[cellI] < 1e-3 )
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