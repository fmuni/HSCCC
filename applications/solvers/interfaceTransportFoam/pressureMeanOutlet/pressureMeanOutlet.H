/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenFOAM Foundation
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

Class
    Foam::pressureMeanOutlet

Group
    grpInletBoundaryConditions

Description
    This boundary condition extrapolates field to the patch using the near-cell
    values and adjusts the distribution to match the specified, optionally
    time-varying, mean value.  This extrapolated field is applied as a
    fixedValue for outflow faces but zeroGradient is applied to inflow faces.

    This boundary condition can be applied to pressure when inletOutlet is
    applied to the velocity so that a zeroGradient condition is applied to the
    pressure at inflow faces where the velocity is specified to avoid an
    unphysical over-specification of the set of boundary conditions.

Usage
    \table
        Property     | Description             | Required    | Default value
        meanValue    | mean value Function1    | yes         |
        phi          | Flux field name         | no          | phi
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        scalar            fixedMeanOutletInlet;
        meanValue       1.0;
    }
    \endverbatim

See also
    Foam::fixedMeanFvPatchField
    Foam::outletInletFvPatchField
    Foam::Function1scalars

SourceFiles
    pressureMeanOutlet.C

\*---------------------------------------------------------------------------*/

#ifndef pressureMeanOutlet_H
#define pressureMeanOutlet_H

#include "outletInletFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class pressureMeanOutlet Declaration
\*---------------------------------------------------------------------------*/

class pressureMeanOutlet
:
    public outletInletFvPatchField<scalar>
{
    // Private data

        //- MeanValue value the field is adjusted to maintain
        autoPtr<Function1<scalar>> meanValue_;

        //- Current time step id to check if updated
        label currTimeIndex_ = -1;


public:

    //- Runtime scalar information
    TypeName("pressureMeanOutlet");


    // Constructors

        //- Construct from patch and internal field
        pressureMeanOutlet
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        pressureMeanOutlet
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given pressureMeanOutlet
        //  onto a new patch
        pressureMeanOutlet
        (
            const pressureMeanOutlet&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        pressureMeanOutlet
        (
            const pressureMeanOutlet&
        );

        //- Construct as copy setting internal field reference
        pressureMeanOutlet
        (
            const pressureMeanOutlet&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return fvPatchField<scalar>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return fvPatchField<scalar>::Clone(*this, iF);
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
