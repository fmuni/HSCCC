#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();


// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);


typedef incompressibleInterPhaseTransportModel
    <
        immiscibleIncompressibleTwoPhaseMixture
    > transportModelType;

autoPtr<transportModelType> turbulence
(
    new transportModelType
    (
        rho, U, phi, rhoPhi, mixture
    )
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

#include "createMRF.H"
#include "createFvOptions.H"
#include "createPorosity.H"

isoAdvection advector(alpha1, phi, U);

//- List of species
List<word> specieNames(mixture.get<List<word>>("species"));

//- Pointers to species fields, diffusion coefficients, and Henry constants
PtrList<volScalarField>     cPtr(specieNames.size());
PtrList<dimensionedScalar>  D1Ptr(specieNames.size()); // phase 1
PtrList<dimensionedScalar>  D2Ptr(specieNames.size()); // phase 2
PtrList<dimensionedScalar>  HPtr(specieNames.size());

forAll(cPtr,ci)
{
    Info << "Reading specie " << specieNames[ci] << endl;

    //- Read specie dictionary
    const dictionary cDict(mixture.subDict(specieNames[ci]));

    //- Create specie field
    cPtr.set
    (
        ci,
        new volScalarField
        (
            IOobject
            (
                specieNames[ci],
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );

    if ( cDict.found("initValue." + mixture.phase1Name()) )
    {
        //- Init according to dictionary
        scalar initPhase1 = cDict.get<scalar>("initValue." + mixture.phase1Name());
        scalar initPhase2 = cDict.get<scalar>("initValue." + mixture.phase2Name());
    
        volScalarField& c = cPtr[ci];
        forAll(c, cellI)
        {
            c[cellI] = alpha1[cellI]*initPhase1 + alpha2[cellI]*initPhase2;
        }
        c.correctBoundaryConditions();
    }
    
    //- Read diffusion coefficients
    D1Ptr.set
    (
        ci,
        new dimensionedScalar
        (
            "D." + mixture.phase1Name(),
            dimViscosity,
            cDict.get<scalar>("D." + mixture.phase1Name())
        )
    );

    D2Ptr.set
    (
        ci,
        new dimensionedScalar
        (
            "D." + mixture.phase2Name(),
            dimViscosity,
            cDict.get<scalar>("D." + mixture.phase2Name())
        )
    );

    //- Read Henry constant
    HPtr.set
    (
        ci,
        new dimensionedScalar
        (
            "H",
            dimless,
            cDict.get<scalar>("H")
        )
    );

}

