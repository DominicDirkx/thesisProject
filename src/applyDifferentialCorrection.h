#ifndef TUDATBUNDLE_APPLYDIFFERENTIALCORRECTION_H
#define TUDATBUNDLE_APPLYDIFFERENTIALCORRECTION_H

#include <map>
#include <iostream>

#include "Eigen/Core"

#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"

Eigen::VectorXd applyDifferentialCorrection(
        const int librationPointNr, const std::string& orbitType,
        const Eigen::VectorXd& initialStateVector,
        double orbitalPeriod, const double massParameter,
        const boost::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > > integratorSettings,
        double maxPositionDeviationFromPeriodicOrbit,
        double maxVelocityDeviationFromPeriodicOrbit,
        const int maxNumberOfIterations = 1000 );


#endif  // TUDATBUNDLE_APPLYDIFFERENTIALCORRECTION_H
