/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "kineticallyLimitedPartitioning.H"
#include "constants.H"
#include "cutCellIso.H"
#include "volPointInterpolation.H"

using namespace Foam::constant;

// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
void Foam::speciesTransferModels::
kineticallyLimitedPartitioning<Thermo, OtherThermo>::updateInterface
(
    const volScalarField& T
)
{
    const fvMesh& mesh = this->mesh_;

    const volScalarField& alpha = this->pair().from();

    scalarField ap
    (
        volPointInterpolation::New(mesh).interpolate(alpha)
    );

    cutCellIso cutCell(mesh, ap);

    forAll(interfaceArea_, celli)
    {
        const label status = cutCell.calcSubCell(celli, isoAlpha_);
        interfaceArea_[celli] = 0;
        if (status == 0) // cell is cut
        {
            interfaceArea_[celli] = mag(cutCell.faceArea())/mesh.V()[celli];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::speciesTransferModels::kineticallyLimitedPartitioning<Thermo, OtherThermo>
::kineticallyLimitedPartitioning
(
    const dictionary& dict,
    const phasePair& pair
)
:
    InterfaceCompositionModel<Thermo, OtherThermo>(dict, pair),
    isoAlpha_(dict.getOrDefault<scalar>("isoAlpha", 0.5)),
    kL_("kL", dimless/dimTime/dimArea, dict),
    K_("K", dimless, dict),
    interfaceArea_
    (
        IOobject
        (
            "interfaceArea",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimless/dimLength, Zero)
    ),
    mDotc_
    (
        IOobject
        (
            "mDotc",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimDensity/dimTime, Zero)
    ),
    Tactivate_("Tactivate", dimTemperature, 0, dict)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::speciesTransferModels::kineticallyLimitedPartitioning<Thermo, OtherThermo>
::Kexp
(
    const volScalarField& T
)
{

    const word speciesName(IOobject::member(this->transferSpecie()));

    const volScalarField& Yfrom = this->fromThermo_.composition().Y(speciesName + "." + this->pair().from().name());
    const volScalarField& Yto = this->toThermo_.composition().Y(speciesName + "." + this->pair().to().name());

    updateInterface(T);

    mDotc_ = kL_ * interfaceArea_ * ( Yto - K_ * Yfrom );

    return tmp<volScalarField>::New(mDotc_);
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::speciesTransferModels::kineticallyLimitedPartitioning<Thermo, OtherThermo>::
KSp
(
    label variable,
    const volScalarField& refValue
)
{
    return nullptr;
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::speciesTransferModels::kineticallyLimitedPartitioning<Thermo, OtherThermo>::
KSu
(
    label variable,
    const volScalarField& refValue
)
{
    return nullptr;
}


//************************************************************************ //
