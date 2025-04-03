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

surfaceScalarField alphaf
(
    IOobject
    (
        "alphaf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    0*fvc::interpolate(alpha1)
);

surfaceScalarField alpha2f
(
    IOobject
    (
        "alpha2f",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    1 - alphaf
);

surfaceScalarField alphaPhi
(
    IOobject
    (
        "alphaPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    phi
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

volVectorField interfaceNormals_0
(
    IOobject
    (
        "interfaceNormals_0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    advector.surf().normal()
);

volVectorField interfaceNormals
(
    IOobject
    (
        "interfaceNormals",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    advector.surf().normal()
);

volVectorField interfaceCenters_0
(
    IOobject
    (
        "interfaceCenters_0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    advector.surf().centre()
);

volVectorField interfaceCenters
(
    IOobject
    (
        "interfaceCenters",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    advector.surf().centre()
);


//- List of species
List<word> specieNames(mixture.get<List<word>>("species"));

//- Pointers to species fields, diffusion coefficients, and Henry constants
PtrList<volScalarField>     c1Ptr(specieNames.size()); // Species in liquid 1
PtrList<volScalarField>     c2Ptr(specieNames.size()); // Species in liquid 2
PtrList<volScalarField>     c11Ptr(specieNames.size()); // Species in liquid 1 averaged in phase 1
PtrList<volScalarField>     c22Ptr(specieNames.size()); // Species in liquid 2 averaged in phase 2
PtrList<dimensionedScalar>  D1Ptr(specieNames.size()); // Diffusion coefficient 1
PtrList<dimensionedScalar>  D2Ptr(specieNames.size()); // Diffusion coefficient 2
PtrList<dimensionedScalar>  KPtr(specieNames.size()); // Partition coefficient
PtrList<dimensionedScalar> klPtr(specieNames.size()); // Transfer rate

forAll(c1Ptr,ci)
{
    Info << "Reading specie " << specieNames[ci] << endl;

    //- Read specie dictionary
    const dictionary cDict(mixture.subDict(specieNames[ci]));

    //- Create specie field
    c1Ptr.set
    (
        ci,
        new volScalarField
        (
            IOobject
            (
                specieNames[ci] + "." + mixture.phase1Name(),
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );

    c2Ptr.set
    (
        ci,
        new volScalarField
        (
            IOobject
            (
                specieNames[ci] + "." + mixture.phase2Name(),
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );

    c11Ptr.set
    (
        ci,
        new volScalarField
        (
            IOobject
            (
                specieNames[ci] + "_ave." + mixture.phase1Name(),
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            c1Ptr[ci]
        )
    );

    c22Ptr.set
    (
        ci,
        new volScalarField
        (
            IOobject
            (
                specieNames[ci] + "_ave." + mixture.phase2Name(),
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            c2Ptr[ci]
        )
    );

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

    //- Read partition coefficient
    KPtr.set
    (
        ci,
        new dimensionedScalar
        (
            "H",
            dimless,
            cDict.get<scalar>("K")
        )
    );


    //- Read transfer rate
    klPtr.set
    (
        ci,
        new dimensionedScalar
        (
            "kl",
            dimLength/dimTime,
            cDict.get<scalar>("kl")
        )
    );

}

