
#include "processorPolyPatch.H"

namespace Foam
{

void bulkInterpolate
(
    const volScalarField& C,
    const boolList& interfaceCells,
    const surfaceScalarField& Cf_2ndOrder,
    surfaceScalarField& Cf
)
{
    const fvMesh& mesh = C.mesh();

    forAll(C, cellI)
    {
        if (interfaceCells[cellI] == 0)
        {
            labelList faces = mesh.cells()[cellI];
            forAll(faces, faceI)
            {
                if (mesh.isInternalFace(faces[faceI]))
                {
                    Cf[faces[faceI]] = Cf_2ndOrder[faces[faceI]];
                }
            }
        }
    }
}

Foam::scalar interpolateOnFaceI
(
    const scalar& alphaf,
    const scalar& COwn,
    const scalar& CNei,
    const scalar& alphaOwn,
    const scalar& alphaNei,
    const scalar& alphaPhiFaceI
)
{
    const scalar alphaTol = 1e-012;
    scalar Cf = 0;
    if (alphaf > alphaTol)
    {
        if (alphaPhiFaceI > 0 && alphaOwn > 0)
        {
            Cf = alphaf * COwn / alphaOwn;
        }
        if (alphaPhiFaceI < 0 && alphaNei > 0)
        {
            Cf = alphaf * CNei / alphaNei;
        }
    }
    else
    {
        Cf = 0;
    }
    return Cf;
}

void interfaceInterpolate
(
    const volScalarField& C,
    const volScalarField& alpha,
    const boolList& interfaceCells,
    const surfaceScalarField& alphaf,
    const surfaceScalarField& alphaPhi,
    surfaceScalarField& Cf
)
{
    const fvMesh& mesh = Cf.mesh();
    const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    forAll(interfaceCells, cellI)
    {
        labelList interfaces = mesh.cells()[cellI];
        if (interfaceCells[cellI] == 1)
        {
            forAll(interfaces, faceI)
            {
                if (mesh.isInternalFace(interfaces[faceI]))
                {
                    const scalar COwn = C[owner[interfaces[faceI]]];
                    const scalar CNei = C[neighbour[interfaces[faceI]]];
                    const scalar alphaOwn = alpha[owner[interfaces[faceI]]];
                    const scalar alphaNei = alpha[neighbour[interfaces[faceI]]];
                    const scalar alphafaceI = alphaf[interfaces[faceI]];
                    const scalar alphaPhiFaceI = alphaPhi[interfaces[faceI]];
                    Cf.primitiveFieldRef()[interfaces[faceI]] = interpolateOnFaceI(alphafaceI, COwn, CNei, alphaOwn, alphaNei, alphaPhiFaceI);
                }
                else
                {
                    const label patchI = boundaryMesh.whichPatch(interfaces[faceI]);        // index of the patch to which faceI belongs to
                    const label bfaceI = boundaryMesh[patchI].whichFace(interfaces[faceI]); // - mesh.nInternalFaces(); // local faceI index
                    if (Cf.boundaryField()[patchI].coupled())
                    {
                        const scalar COwn = C.boundaryField()[patchI].patchInternalField()()[bfaceI];
                        const scalar CNei = C.boundaryField()[patchI].patchNeighbourField()()[bfaceI];
                        const scalar alphaOwn = alpha.boundaryField()[patchI].patchInternalField()()[bfaceI];
                        const scalar alphaNei = alpha.boundaryField()[patchI].patchNeighbourField()()[bfaceI];
                        const scalar alphafaceI = alphaf.boundaryField()[patchI][bfaceI];
                        const scalar alphaPhiFaceI = alphaPhi.boundaryField()[patchI][bfaceI];
                        Cf.boundaryFieldRef()[patchI][bfaceI] = interpolateOnFaceI(alphafaceI, COwn, CNei, alphaOwn, alphaNei, alphaPhiFaceI);
                    }
                    // else // to be implemented
                    // {
                    //     Cf.boundaryFieldRef()[patch][bfaceI] = alphaf.boundaryField()[patch][bfaceI]*C.primitiveField()[faceCells[bFaceI]]/alpha.primitiveField()[cellI];
                    // }
                }
            }
        }
    }
}



Foam::scalar interpolateFaceConcentrationOnFaceI
(
    const scalar& alphaf,
    const scalar& COwn,
    const scalar& CNei,
    const scalar& alphaOwn,
    const scalar& alphaNei,
    const scalar& alphaPhiFaceI
)
{
    const scalar alphaTol = 1e-012;
    scalar Cf = 0;
    if (alphaf > alphaTol)
    {
        if (alphaPhiFaceI > 0 && alphaOwn > 0)
        {
            Cf = alphaf * COwn;
        }
        if (alphaPhiFaceI < 0 && alphaNei > 0)
        {
            Cf = alphaf * CNei;
        }
    }
    else
    {
        Cf = 0;
    }
    return Cf;
}


void syncProcPatches
(
    surfaceScalarField& Cf,
    const boolList& interfaceCells
)
{
    // references to the mesh and patches
    const fvMesh& mesh = Cf.mesh();
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();
    const fvPatchList& patches = mesh.boundary();
    // list containing the patches indexes
    DynamicList<label> procPatchLabels(pbm.size());
    // list for the communication of the surfaceScalarField
    List<DynamicList<label>> surfaceCellFacesOnProcPatches(0);
    List<DynamicList<scalar>> pFaceValues(0);
    // list used to verify if the patch needs information from the other processor
    List<DynamicList<bool>> requiresNeighbData(0);
    surfaceCellFacesOnProcPatches.resize(pbm.size());
    pFaceValues.resize(pbm.size());
    requiresNeighbData.resize(pbm.size());
    // Append all processor patch labels to the list
    forAll(patches, patchI)
    {
        const fvPatch& thisPatch = patches[patchI];
        const polyPatch& thisPolyPatch = thisPatch.patch();
        if (isA<processorPolyPatch>(thisPolyPatch) && thisPolyPatch.size() > 0)
        {
            procPatchLabels.append(patchI);
            forAll(thisPolyPatch, faceI)
            {
                surfaceCellFacesOnProcPatches[patchI].append(thisPolyPatch.start() + faceI);
                const label cellI = thisPatch.faceCells()[faceI];
                if (interfaceCells[cellI] == 0)
                {
                    requiresNeighbData[patchI].append(1);
                }
                else
                {
                    requiresNeighbData[patchI].append(0);
                }
            }
        }
    }

    if (Pstream::parRun())
    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // Send
        forAll(procPatchLabels, i)
        {
            const label patchI = procPatchLabels[i];
            const processorPolyPatch &procPatch =
                refCast<const processorPolyPatch>(pbm[patchI]);
            UOPstream toNbr(procPatch.neighbProcNo(), pBufs);
            forAll(pbm[patchI], faceI)
            {
                pFaceValues[patchI].append(Cf.boundaryField()[patchI][faceI]);
            }
            const List<label>& surfCellFacesOnProcPatch =
                surfaceCellFacesOnProcPatches[patchI];
            const List<scalar>& faceValues =
                pFaceValues[patchI];
            toNbr << surfCellFacesOnProcPatch << faceValues;
        }

        pBufs.finishedSends();

        // Receive and combine
        forAll(procPatchLabels, patchLabelI)
        {
            const label patchI = procPatchLabels[patchLabelI];
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(pbm[patchI]);
            List<label> faceIDs;
            List<scalar> nbrCfs;
            UIPstream fromNeighb(procPatch.neighbProcNo(), pBufs);
            fromNeighb >> faceIDs >> nbrCfs;
            scalarField& localFaceValues = Cf.boundaryFieldRef()[patchI];
            forAll(faceIDs, i)
            {
                if (requiresNeighbData[patchI][i] == 1)
                {
                    localFaceValues[i] = nbrCfs[i];
                }
            }
        }
    }
    surfaceCellFacesOnProcPatches.clear();
    pFaceValues.clear();
    requiresNeighbData.clear();
}

void phaseAverageInterpolate
(
    const volScalarField& C,
    const volScalarField& alpha,
    const boolList& interfaceCells,
    const surfaceScalarField& alphaf,
    const surfaceScalarField& alphaPhi,
    surfaceScalarField& Cf
)
{
    const fvMesh& mesh = Cf.mesh();
    const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    forAll(interfaceCells, cellI)
    {
        labelList interfaces = mesh.cells()[cellI];
        if (interfaceCells[cellI] == 1)
        {
            forAll(interfaces, faceI)
            {
                if (mesh.isInternalFace(interfaces[faceI]))
                {
                    const scalar COwn = C[owner[interfaces[faceI]]];
                    const scalar CNei = C[neighbour[interfaces[faceI]]];
                    const scalar alphaOwn = alpha[owner[interfaces[faceI]]];
                    const scalar alphaNei = alpha[neighbour[interfaces[faceI]]];
                    const scalar alphafaceI = alphaf[interfaces[faceI]];
                    const scalar alphaPhiFaceI = alphaPhi[interfaces[faceI]];
                    Cf.primitiveFieldRef()[interfaces[faceI]] = interpolateFaceConcentrationOnFaceI(alphafaceI, COwn, CNei, alphaOwn, alphaNei, alphaPhiFaceI);
                }
                else
                {
                    const label patchI = boundaryMesh.whichPatch(interfaces[faceI]);        // index of the patch to which faceI belongs to
                    const label bfaceI = boundaryMesh[patchI].whichFace(interfaces[faceI]); // - mesh.nInternalFaces(); // local faceI index
                    if (Cf.boundaryField()[patchI].coupled())
                    {
                        const scalar COwn = C.boundaryField()[patchI].patchInternalField()()[bfaceI];
                        const scalar CNei = C.boundaryField()[patchI].patchNeighbourField()()[bfaceI];
                        const scalar alphaOwn = alpha.boundaryField()[patchI].patchInternalField()()[bfaceI];
                        const scalar alphaNei = alpha.boundaryField()[patchI].patchNeighbourField()()[bfaceI];
                        const scalar alphafaceI = alphaf.boundaryField()[patchI][bfaceI];
                        const scalar alphaPhiFaceI = alphaPhi.boundaryField()[patchI][bfaceI];
                        Cf.boundaryFieldRef()[patchI][bfaceI] = interpolateFaceConcentrationOnFaceI(alphafaceI, COwn, CNei, alphaOwn, alphaNei, alphaPhiFaceI);
                    }
                    // else // to be implemented
                    // {
                    //     Cf.boundaryFieldRef()[patch][bfaceI] = alphaf.boundaryField()[patch][bfaceI]*C.primitiveField()[faceCells[bFaceI]]/alpha.primitiveField()[cellI];
                    // }
                }
            }
        }
    }
}

}


