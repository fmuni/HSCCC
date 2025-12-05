/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "pressureMeanOutlet.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::pressureMeanOutlet::pressureMeanOutlet
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    outletInletFvPatchField<scalar>(p, iF),
    meanValue_()
{}



Foam::pressureMeanOutlet::pressureMeanOutlet
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    outletInletFvPatchField<scalar>(p, iF),
    meanValue_(Function1<scalar>::New("meanValue", dict, &this->db()))
{
    this->phiName_ = dict.getOrDefault<word>("phi", "phi");

    this->readValueEntry(dict, IOobjectOption::MUST_READ);

    this->refValue() = *this;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}



Foam::pressureMeanOutlet::pressureMeanOutlet
(
    const pressureMeanOutlet& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    outletInletFvPatchField<scalar>(ptf, p, iF, mapper),
    meanValue_(ptf.meanValue_.clone())
{}



Foam::pressureMeanOutlet::pressureMeanOutlet
(
    const pressureMeanOutlet& ptf
)
:
    outletInletFvPatchField<scalar>(ptf),
    meanValue_(ptf.meanValue_.clone())
{}



Foam::pressureMeanOutlet::pressureMeanOutlet
(
    const pressureMeanOutlet& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    outletInletFvPatchField<scalar>(ptf, iF),
    meanValue_(ptf.meanValue_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::pressureMeanOutlet::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // This prevents instabilities in the boundary condition.
    // Essentially, it makes sure the pressure BC is computed using the old
    // pressure field. This introduces some additional diffusion in the solution
    // because of the first order temporal approximation, but it allows convergence
    // in the sequence of pressure correctors.
    label timeIndex = this->db().time().timeIndex();
    if ( timeIndex == currTimeIndex_ )
    {
        return;
    }
    else
    {
        currTimeIndex_ = timeIndex;
    }

    Info << "Computing new pressure BC\n";
    const scalar t = this->db().time().timeOutputValue();
    scalar meanValue = meanValue_->value(t);

    Field<scalar> newValues(this->patchInternalField());

    scalar meanValuePsi =
        gSum(this->patch().magSf()*newValues)
       /gSum(this->patch().magSf());


    newValues += (meanValue - meanValuePsi);

    this->refValue() = newValues;

    outletInletFvPatchField<scalar>::updateCoeffs();
}



void Foam::pressureMeanOutlet::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    meanValue_->writeData(os);
    fvPatchField<scalar>::writeValueEntry(os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        pressureMeanOutlet
    );
}

// ************************************************************************* //


// ************************************************************************* //
