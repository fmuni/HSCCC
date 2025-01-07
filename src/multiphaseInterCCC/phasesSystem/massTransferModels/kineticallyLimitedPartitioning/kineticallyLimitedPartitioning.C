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
    kL_("kL", dimLength/dimTime, dict),
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
    Tactivate_("Tactivate", dimTemperature, 0, dict),
    Sexp_(this->transferSpecie().size())
{
    wordList species = this->transferSpecie();
    forAll(Sexp_,specieI)
    {
        Sexp_.set
        (
            specieI,
            volScalarField::New
            (
                "Sexp." + species[specieI],
                IOobject::NO_REGISTER,
                this->mesh_,
                dimensionedScalar(dimDensity/dimTime, Zero)
            )
        );
    }
    
    //- Check if species exist
    const PtrList<volScalarField>& Yfrom = this->fromThermo_.composition().Y();
    const PtrList<volScalarField>& Yto = this->toThermo_.composition().Y();
    forAll(species,specieI)
    {
        bool found = false;

        forAll(Yfrom,Yi)
        {
            if(Yfrom[Yi].name() == species[specieI] + "." + this->pair().from().name())
            {
                found = true;
            }
        }

        forAll(Yto,Yi)
        {
            if(Yto[Yi].name() == species[specieI] + "." + this->pair().to().name())
            {
                found = true;
            }
        }
        
        if(!found)
        {
            FatalErrorInFunction << "Species " << species[specieI] << " must be present in both phases." << exit(FatalError);
        }


    }

    

}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::speciesTransferModels::kineticallyLimitedPartitioning<Thermo, OtherThermo>
::Kexp
(
    const volScalarField& T
)
{

    const wordList speciesNames(this->transferSpecie());
    updateInterface(T);
    mDotc_ *= 0.;
    const volScalarField from
    (
        clamp(this->pair().from(), zero_one{})
    );
    
    forAll(speciesNames, specieI)
    {
        word speciesName = speciesNames[specieI];
        word speciesNameFrom = speciesName + "." + this->pair().from().name(); 
        word speciesNameTo = speciesName + "." + this->pair().to().name();  
        const PtrList<volScalarField>& Ysfrom = this->fromThermo_.composition().Y();
        const PtrList<volScalarField>& Ysto = this->toThermo_.composition().Y();    
        
        label idFrom = -1;
        label idTo = -1;
        
        forAll(this->fromThermo_.composition().Y(),Yi)
        {
            if(Ysfrom[Yi].name() == speciesName + "." + this->pair().from().name())
            {
                idFrom = Yi;
                break;
            }
        }
        forAll(this->toThermo_.composition().Y(),Yi)
        {
            if(Ysto[Yi].name() == speciesName + "." + this->pair().to().name())
            {
                idTo = Yi;
                break;
            }
        }

        if(idFrom == -1 || idTo == -1)
        {
            FatalErrorInFunction << "Species " << speciesName << " must be present in both phases." << exit(FatalError);
        }        

        const volScalarField& Yfrom = Ysfrom[idFrom];
        const volScalarField& Yto = Ysto[idTo];

        // Sexp is the mass source term for the species.
        // Here "to" is the mobile phase and "from" is the stationary phase. The partition coefficient is 
        // defined as K = Yto/Ytfrom at equilibrium.
        // The first line in the formula comes from the conversion from surface integral to volumetric term
        // The second line is the model equation.
        // The last line ensures that phase from does not acquire mass from phase to.
        Sexp_[specieI] = this->pair().from().rho() * from
                        * kL_ * interfaceArea_ * ( Yto - K_ * Yfrom ) 
                        * pos( K_ * Yfrom - Yto );
                        
        mDotc_ +=  Sexp_[specieI]; // Multiply by rho to convert to mass source term
        

    }

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

template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::speciesTransferModels::kineticallyLimitedPartitioning<Thermo, OtherThermo>::
vdotSpecie(word specieName) const
{

    const wordList speciesNames(this->transferSpecie());
    
    forAll(speciesNames, specieI)
    {
        // Ensure consistency between phases
        if(specieName == speciesNames[specieI] + "." + this->pair().from().name())
        {
            return tmp<volScalarField>::New(Sexp_[specieI]/this->pair().from().rho() );
        }
        
        if(specieName == speciesNames[specieI] + "." + this->pair().to().name())
        {
            return tmp<volScalarField>::New(-Sexp_[specieI]/this->pair().from().rho());
        }                
    }
    
    return nullptr;

}


//************************************************************************ //
