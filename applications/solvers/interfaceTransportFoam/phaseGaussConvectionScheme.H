
#ifndef phasephaseGaussConvectionScheme_H
#define phasephaseGaussConvectionScheme_H

#include "convectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class phaseGaussConvectionScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class phaseGaussConvectionScheme
:
    public fv::convectionScheme<Type>
{
    // Private data

        tmp<surfaceInterpolationScheme<Type>> tinterpScheme_;

        const volScalarField&       alpha_;
        const volScalarField&       alphaf_;


    // Private Member Functions

        //- No copy construct
        phaseGaussConvectionScheme(const phaseGaussConvectionScheme&) = delete;

        //- No copy assignment
        void operator=(const phaseGaussConvectionScheme&) = delete;

        GeometricField<Type, fvPatchField, volMesh> phaseAverage
        (
            const GeometricField<Type, fvPatchField, volMesh>& ,
            const scalar
        ) const; 


public:

    //- Runtime type information
    TypeName("phaseGauss");


    // Constructors

        //- Constructors
        phaseGaussConvectionScheme
        (
            const fvMesh& mesh,
            const volScalarField& alpha,          
            const surfaceScalarField& alphaf,
            const surfaceScalarField& alphaPhi,
            const tmp<surfaceInterpolationScheme<Type>>& scheme
        )
        :
            convectionScheme<Type>(mesh, alphaPhi),
            tinterpScheme_(scheme),
            alpha_(alpha),
            alphaf_(alphaf)
        {}

        //- Construct from flux and Istream
        phaseGaussConvectionScheme
        (
            const fvMesh& mesh,
            const volScalarField& alpha,          
            const surfaceScalarField& alphaf,
            const surfaceScalarField& alphaPhi,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            convectionScheme<Type>(mesh, faceFlux),
            tinterpScheme_
            (
                surfaceInterpolationScheme<Type>::New(mesh, alphaPhi, is)
            ),
            alpha_(alpha),
            alphaf_(alphaf)
        {}


    // Member Functions

        const surfaceInterpolationScheme<Type>& interpScheme() const;

        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> interpolate
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> flux
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<fvMatrix<Type>> fvmDiv
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDiv
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
