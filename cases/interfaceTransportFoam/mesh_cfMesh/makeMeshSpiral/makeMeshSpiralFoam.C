#include "fvCFD.H"
#include "dynamicFvMesh.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"


    Info<< "Reading makeMeshSpiralDict\n" << endl;

    IOdictionary mdict
    (
        IOobject
        (
            "blockMeshDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
    
    scalar R(readScalar(mdict.lookup("R"))); //- Radius of the coil
    scalar pitch(readScalar(mdict.lookup("p"))); //- Coil pitch
    scalar twoPi = 2.0*3.14159265359;
    scalar pOverTwoPi = pitch/twoPi;
    
    pointField newPts(mesh.points());
    
    forAll(newPts,pI)
    {
        //- The mesh is supposed to be oriented in x
        //  Plane of coil rotation is x-y
        
        point& p = newPts[pI];
        point pNew = p;
        
        scalar theta = p.z()/(R);
        
        pNew[1] = (R + p[0])*Foam::cos(theta);
        pNew[0] = (R + p[0])*Foam::sin(theta);
        pNew[2] = p[1] + theta*pOverTwoPi;
        p = pNew;
        
    }
    
    //- Moving points
    mesh.movePoints(newPts);
    mesh.update();
    mesh.removeFiles();
    if (!mesh.write())
    {
        FatalErrorInFunction
            << "Failed writing polyMesh."
            << exit(FatalError);
    }



    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
