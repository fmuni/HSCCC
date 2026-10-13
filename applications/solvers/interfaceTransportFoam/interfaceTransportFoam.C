/*---------------------------------------------------------------------------*\
Application
    interfaceTransportIsoFoam

Author
    Federico Municchi. NREL.

Description
    Solver based on interIsoFoam for solving problems where two immiscible
    carrier fluids exchange species with different solubility across the interface.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "isoAdvection.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "incompressibleInterPhaseTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "dynamicRefineFvMesh.H"
#include <math.h>
#include "computePlanetaryAcceleration.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver based on interIsoFoam for solving problems where two immiscible"
        " carrier fluids exchange species with different solubility. The differece.\n"
        "n dissolution is modeled as an Henry jump condition. "
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"
    #include "createAlphaFluxes.H"

    #include "porousCourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "porousCourantNo.H"
        #include "porousAlphaCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        interfaceCenters_0 = interfaceCenters;
        interfaceNormals_0 = interfaceNormals;
        boolList interfaceCells_0 = advector.surf().interfaceCell();

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                if (isA<dynamicRefineFvMesh>(mesh))
                {
                    advector.surf().reconstruct();
                }

                mesh.update();

                if (mesh.changing())
                {
                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    if (isA<dynamicRefineFvMesh>(mesh))
                    {
                        advector.surf().mapAlphaField();
                        alpha2 = 1.0 - alpha1;
                        alpha2.correctBoundaryConditions();
                        rho == alpha1*rho1 + alpha2*rho2;
                        rho.correctBoundaryConditions();
                        rho.oldTime() = rho;
                        alpha2.oldTime() = alpha2;
                    }

                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);

                        mixture.correct();
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            #include "alphaControls.H"

            if (algebraicVOF)
            {
                #include "alphaEqnSubCycle_MULES.H"
            }
            else
            {
                #include "alphaEqnSubCycle.H"
            }

            interfaceCenters = advector.surf().centre();
            interfaceNormals = advector.surf().normal();

            mixture.correct();

            if (pimple.frozenFlow())
            {
                continue;
            }

            #include "computePlanetaryForce.H"

            phiFPrelax = phiFPrelax*(1. - omegap) + phiFP * omegap;

            #include "UEqn.H"
            #include "cEqn.H"

            // --- Pressure corrector loop
            label corrPISO = 1;
            while (corrPISO <= pimple.nCorrPISO())
            {
                #include "pEqn.H"

                if (mag(sumLocalContErr) < contErrTolPISO)
                {
                    Info<< "PISO converged" << endl;
                    break;
                }

                corrPISO++;

            }

            cumulativeContErr += globalContErr;
            Info<< "Cumulative continuity errors= " << cumulativeContErr
                << endl;

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
