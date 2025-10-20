/* Helper function for computing planetary acceleration */
#pragma once
#include "fvCFD.H"
#include <cmath>

using namespace Foam;

inline vector computePlanetaryAcceleration(
    scalar x,
    scalar y,
    scalar R,
    scalar omega,
    scalar t
)
{

    scalar phi = std::atan2(y,x);
    scalar r = std::sqrt(sqr(x) + sqr(y));

    scalar beta = r/R;

    scalar axR = R*sqr(omega)*( std::cos(omega*t) + 4.*beta*std::cos(2.*omega*t + phi) );
    scalar ayR = R*sqr(omega)*( std::sin(omega*t) + 4.*beta*std::sin(2.*omega*t + phi) );

    scalar axr = axR*std::cos(2.*omega*t) + ayR*std::sin(2.*omega*t);
    scalar ayr = -axR*std::sin(2.*omega*t) + ayR*std::cos(2.*omega*t);

    return vector(axr, ayr, 0.);
}
