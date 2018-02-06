#ifndef TUDATBUNDLE_RICHARDSONTHIRDORDERAPPROXIMATION_H
#define TUDATBUNDLE_RICHARDSONTHIRDORDERAPPROXIMATION_H



#include <string>

#include "Tudat/Basics/basicTypedefs.h"

#include "cr3bpPeriodicOrbitTypes.h"

std::pair< Eigen::Vector6d, double > richardsonApproximationLibrationPointPeriodicOrbit(
        tudat::cr3bp::CR3BPPeriodicOrbitTypes orbitType, int librationPointNr,
        double amplitude, double n= 1.0 );



#endif  // TUDATBUNDLE_RICHARDSONTHIRDORDERAPPROXIMATION_H
