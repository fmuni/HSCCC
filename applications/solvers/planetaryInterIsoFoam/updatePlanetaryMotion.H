#pragma once

void updatePlanetaryMotion(
	volVectorField& g,
	volScalarField& gh,
	surfaceScalarField& ghf,
	scalar R,
	scalar omega,
	scalar phase,
	vector gvec,
	scalar t,
	scalar t0  
)
{
    //- Ramp effect
    omega *= min(t/(t0 + 1e-16),1.0);
    
	// Pre-compute trigonometric expression
	scalar R_omgsqr_cos_omegat = R*Foam::sqr(omega)*Foam::cos(omega*t);
	scalar R_omgsqr_sin_omegat = R*Foam::sqr(omega)*Foam::sin(omega*t);
	scalar _4_R_omgsqr_cos_2omegat_phase = 4.*R*Foam::sqr(omega)*Foam::cos(2.*omega*t + phase);
	scalar _4_R_omgsqr_sin_2omegat_phase = 4.*R*Foam::sqr(omega)*Foam::sin(2.*omega*t + phase);
	scalar cos_2omegat = Foam::cos(2.*omega*t);
	scalar sin_2omegat = Foam::sin(2.*omega*t);

	// Update g
	forAll(g, cellI)
	{
	    vector cc = g.mesh().C()[cellI];
	    scalar r = Foam::sqrt( Foam::sqr(cc.x()) + Foam::sqr(cc.y()) ); 

	    // Acceleration in the Solar frame of reference
	    scalar aRx = -R_omgsqr_cos_omegat + r*_4_R_omgsqr_cos_2omegat_phase;
	    scalar aRy = -R_omgsqr_sin_omegat + r*_4_R_omgsqr_sin_2omegat_phase;
	    
	    // Acceleration in the stationary fluid domain
	    vector acc = gvec;
	    acc.x() += aRx*cos_2omegat + aRy*sin_2omegat;
	    acc.y() += -aRx*sin_2omegat + aRy*cos_2omegat;
	    
	    g[cellI] = acc;
	}
/*	
	// Make sure that boundaries are consistent
	volVectorField::Boundary& gBf = g.boundaryFieldRef();	
	forAll(gBf, patchI)
	{
	    fvPatchField<Vector<double> >& gBfi = const_cast<fvPatchField<Vector<double> >&>(gBf[patchI]);
	    forAll(gBfi, faceI)
	    {
            vector cf = gBf[faceI].patch().Cf()[faceI];
            scalar r = Foam::sqrt( Foam::sqr(cf.x()) + Foam::sqr(cf.y()) );
            
            // Acceleration in the Solar frame of reference
            scalar aRx = -R_omgsqr_cos_omegat + r*_4_R_omgsqr_cos_2omegat_phase;
            scalar aRy = -R_omgsqr_sin_omegat + r*_4_R_omgsqr_sin_2omegat_phase;
            
            // Acceleration in the stationary fluid domain
            vector acc = gvec;
            acc.x() += aRx*cos_2omegat + aRy*sin_2omegat;
            acc.y() += -aRx*sin_2omegat + aRy*cos_2omegat;
            
            // Update ghf
            gBfi[faceI] = acc;	        
	    }
	}
*/	
	// Update gh
	gh = g & g.mesh().C();
	
	// Update ghf
	forAll(ghf, faceI)
	{
	    vector cf = ghf.mesh().Cf()[faceI];
	    scalar r = Foam::sqrt( Foam::sqr(cf.x()) + Foam::sqr(cf.y()) );
	    
	    // Acceleration in the Solar frame of reference
	    scalar aRx = -R_omgsqr_cos_omegat + r*_4_R_omgsqr_cos_2omegat_phase;
	    scalar aRy = -R_omgsqr_sin_omegat + r*_4_R_omgsqr_sin_2omegat_phase;
	    
	    // Acceleration in the stationary fluid domain
	    vector acc = gvec;
	    acc.x() += aRx*cos_2omegat + aRy*sin_2omegat;
	    acc.y() += -aRx*sin_2omegat + aRy*cos_2omegat;
	    
	    // Update ghf
	    ghf[faceI] = acc&cf;
	    
	}
}
