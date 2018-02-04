#ifndef TUDATBUNDLE_RICHARDSONTHIRDORDERAPPROXIMATION_H
#define TUDATBUNDLE_RICHARDSONTHIRDORDERAPPROXIMATION_H



#include <string>

#include "Tudat/Basics/basicTypedefs.h"

std::pair< Eigen::Vector6d, double > richardsonApproximationLibrationPointPeriodicOrbit( std::string orbitType, int librationPointNr,
                                                   double amplitude, double n= 1.0 );



#endif  // TUDATBUNDLE_RICHARDSONTHIRDORDERAPPROXIMATION_H
