volVectorField g
(
    IOobject
    (
        "g",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("g",dimVelocity/dimTime, Zero)
);

Info<< "Reading planetaryMotionDict\n" << endl;

IOdictionary planetaryMotionDict
(
    IOobject
    (
        "planetaryMotionDict",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// Read solar radius
dimensionedScalar R("R", dimLength, planetaryMotionDict);

// Read motion frequency
dimensionedScalar omega("omega", dimless/dimTime, planetaryMotionDict);

// Read phase
dimensionedScalar pm_phase("phase", dimless, planetaryMotionDict);

Info << "Reading gravity\n";
const meshObjects::gravity& g_vec = meshObjects::gravity::New(runTime);

// Read ramp time
dimensionedScalar rampTime("rampTime",dimTime,planetaryMotionDict);


const int oldLocal = volVectorField::Boundary::localConsistency;
volVectorField::Boundary::localConsistency = 0;

surfaceScalarField ghf("ghf", fvc::interpolate(g) & (mesh.Cf() ));

volScalarField gh("gh",g & mesh.C());

volVectorField::Boundary::localConsistency = oldLocal;

updatePlanetaryMotion (
	g, 
	gh, 
	ghf, 
	R.value(), 
	omega.value(), 
	pm_phase.value(), 
	g_vec.value(), 
	runTime.time().value(),
	rampTime.value()
);




