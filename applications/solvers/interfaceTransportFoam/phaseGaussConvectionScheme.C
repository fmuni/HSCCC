/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "phaseGaussConvectionScheme.H"
#include "fvcSurfaceIntegrate.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
template<class Type>
GeometricField<Type, fvPatchField, volMesh> 
phaseGaussConvectionScheme<Type>::phaseAverage
(
    const GeometricField<Type, fvPatchField, volMesh>& vf, 
    const scalar tolerance
) const
{
    GeometricField<Type, fvPatchField, volMesh>  vfave = vf;
    forAll(vfave, cellI)
    {
        if ( alpha_[cellI]  < tolerance )
        {
            vfave[cellI] = 0.;
        }
        else
        {
            vfave[cellI] = vf[cellI] / alpha_[cellI];
        }
    }

    vfave.correctBoundaryConditions;

    return vfave;
}

template<class Type>
const surfaceInterpolationScheme<Type>&
phaseGaussConvectionScheme<Type>::interpScheme() const
{
    return tinterpScheme_();
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
phaseGaussConvectionScheme<Type>::interpolate
(
    const surfaceScalarField&,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{

    GeometricField<Type, fvPatchField, volMesh> vfave = phaseAverage(vf,alpha_,scalar(1e-4));

    return tinterpScheme_().interpolate(vfave)*alphaf_;
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
phaseGaussConvectionScheme<Type>::flux
(
    const surfaceScalarField& alphaPhi,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    GeometricField<Type, fvPatchField, volMesh> vfave = phaseAverage(vf,alpha_,scalar(1e-4));
    return alphaPhi*tinterpScheme_().interpolate(vfave);
}


template<class Type>
tmp<fvMatrix<Type>>
phaseGaussConvectionScheme<Type>::fvmDiv
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{


    tmp<surfaceScalarField> tweights = tinterpScheme_().weights(vf);
    const surfaceScalarField& weights = tweights();

    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            faceFlux.dimensions()*vf.dimensions()
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    forAll(fvm.lower(),ownerI)
    {
        scalar coeff = alpha_[alpha_.mesh().owner()[ownerI]];

        if( coeff < 1e-3 )
        {
            coeff = 0;
        }
        else
        {
            coeff = 1.0/coeff;
        }

        fvm.lower()[ownerI] = -coeff*weights.primitiveField()[ownerI]*faceFlux.primitiveField()[ownerI];

    }

    forAll(fvm.upper(),neighI)
    {
        scalar coeff = alpha_[alpha_.mesh().neighbour()[neighI]];

        if( coeff < 1e-3 )
        {
            coeff = 0;
        }
        else
        {
            coeff = 1.0/coeff;
        }

        fvm.upper()[neighI] =  coeff* ( 1.0 - weights.primitiveField()[neighI])*faceFlux.primitiveField()[neighI];

    }

    fvm.negSumDiag();

    GeometricField<Type, fvPatchField, volMesh> vfave = phaseAverage(vf,alpha_,scalar(1e-4));
    forAll(vf.boundaryField(), patchi)
    {
        const fvPatchField<Type>& psf = vfave.boundaryField()[patchi];
        const fvsPatchScalarField& patchFlux = faceFlux.boundaryField()[patchi];
        const fvsPatchScalarField& pw = weights.boundaryField()[patchi];

        fvm.internalCoeffs()[patchi] = patchFlux*psf.valueInternalCoeffs(pw);
        fvm.boundaryCoeffs()[patchi] = -patchFlux*psf.valueBoundaryCoeffs(pw);
    }

    if (tinterpScheme_().corrected())
    {
        fvm += fvc::surfaceIntegrate(faceFlux*tinterpScheme_().correction(vfave));
    }

    return tfvm;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
phaseGaussConvectionScheme<Type>::fvcDiv
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tConvection
    (
        fvc::surfaceIntegrate(flux(faceFlux, vf))
    );

    tConvection.ref().rename
    (
        "convection(" + faceFlux.name() + ',' + vf.name() + ')'
    );

    return tConvection;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
