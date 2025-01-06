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

Class
    Foam::meltingEvaporationModels::kineticallyLimitedPartitioning

Description
    Non-equilibrium model where the specie is partitioned between the 
    two phases.

    The interface condition is given by	
    \f[
        \hat{m} = kL \rho_{to} D_{to} ( Y_{to} - K Y_{from}  )
    \f]

    where:
    \vartable
        \hat{m}         | mass flux rate [kg/s/m2]
        \rho_{to}       | density of the phase traveling to
        D_{to}          | diffusion coefficient of the phase traveling to
        kL              | non-equilibrium limiting factor
        K 		| partition coefficient
        Y_{to}          | mass fraction at the phase to
        Y_{from}        | mass fraction at the phase from
    \endvartable

Usage
    Example usage:
    \verbatim
    massTransferModel
    (
        (aqueous to organic)
        {
            // Mandatory entries
            type                kineticallyLimitedPartitioning;
            species             chem.organic;
            kL                  1;
            K			0.5;
            
            // Optional entries
            isoAlpha            <scalar>;
            Tactivate           <scalar>;

            // Inherited entries
            ...
        }
    );
    \endverbatim

    where:
    \table
      Property     | Description                        | Type | Reqd | Deflt
      type         | Type name: kineticallyLimitedPartitioning | word | yes  | -
      saturationPressure | Saturation model             | dict | yes  | -
      isoAlpha    | Iso-alpha for interface            | scalar | no  | 0.5
      kL           | limiting factor                  | scalar | yes | -
      K  | partition coefficient             | scalar | yes  | -
    \endtable

    The inherited entries are elaborated in:
      - \link InterfaceCompositionModel.H \endlink
      - \link saturationModel.H \endlink

SourceFiles
    kineticallyLimitedPartitioning.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_meltingEvaporationModels_kineticallyLimitedPartitioning_H
#define Foam_meltingEvaporationModels_kineticallyLimitedPartitioning_H

#include "InterfaceCompositionModel.H"
#include "saturationModel.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * * * *//

namespace Foam
{

// Forward Declarations
class phasePair;

namespace speciesTransferModels
{

/*---------------------------------------------------------------------------*\
                    Class kineticallyLimitedPartitioning Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo, class OtherThermo>
class kineticallyLimitedPartitioning
:
    public InterfaceCompositionModel<Thermo, OtherThermo>
{
    // Private Data

        //- Saturation pressure model
        autoPtr<saturationModel> saturationModelPtr_;

        //- Interface Iso-value
        scalar isoAlpha_;

        //- Limiting factor
        dimensionedScalar kL_;

        //- Partition coefficient
        dimensionedScalar K_;

        //- Interface area
        volScalarField interfaceArea_;

        //- Mass source
        volScalarField mDotc_;

        //- Activation temperature
        const dimensionedScalar Tactivate_;


    // Private Member Functions

        //- Update interface
        void updateInterface(const volScalarField& T);
        
        // Explicit source terms for each specie
        PtrList<volScalarField> Sexp_;        


public:

    //- Runtime type information
    TypeName("kineticallyLimitedPartitioning");


    // Constructors

        //- Construct from components
        kineticallyLimitedPartitioning
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~kineticallyLimitedPartitioning() = default;


    // Member Functions

        //- Explicit total mass transfer coefficient
        virtual tmp<volScalarField> Kexp
        (
            const volScalarField& field
        );

        //- Implicit mass transfer coefficient
        virtual tmp<volScalarField> KSp
        (
            label modelVariable,
            const volScalarField& field
        );

        //- Explicit mass transfer coefficient
        virtual tmp<volScalarField> KSu
        (
            label modelVariable,
            const volScalarField& field
        );

        //- Return Tactivate
        virtual const dimensionedScalar& Tactivate() const noexcept
        {
            return Tactivate_;
        }

        //- Add/subtract alpha*div(U) as a source term
        //- for alpha, substituting div(U) = mDot(1/rho1 - 1/rho2)
        virtual bool includeDivU() const noexcept
        {
            return true;
        }
        
	virtual tmp<volScalarField> vdotSpecie(word specieName) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace meltingEvaporationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "kineticallyLimitedPartitioning.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
