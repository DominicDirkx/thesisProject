#ifndef TUDATBUNDLE_CR3BPPERIODICORBITTYPES_H
#define TUDATBUNDLE_CR3BPPERIODICORBITTYPES_H



#include <string>
#include <vector>
#include <map>
#include <iostream>

#include <boost/function.hpp>
#include <boost/bind.hpp>

#include <Eigen/Core>

#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/Astrodynamics/Propagators/stateDerivativeCircularRestrictedThreeBodyProblem.h"

namespace tudat
{

namespace cr3bp
{

enum CR3BPPeriodicOrbitTypes
{
    horizontal_lyapunov_orbit,
    vertical_lyapunov_orbit,
    halo_orbit,
    axial_orbit
};

}

}

#endif  // TUDATBUNDLE_CR3BPPERIODICORBITTYPES_H
