#ifndef TUDATBUNDLE_APPLYDIFFERENTIALCORRECTION_H
#define TUDATBUNDLE_APPLYDIFFERENTIALCORRECTION_H

#include <map>
#include <iostream>

#include "Eigen/Core"

#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"

#include "cr3bpPeriodicOrbits.h"

Eigen::VectorXd applyDifferentialCorrection(
        const Eigen::VectorXd& initialStateVector,
        double orbitalPeriod,
        const boost::shared_ptr< tudat::cr3bp::CR3BPPeriodicOrbitModel > periodicOrbitModel,
        const boost::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const int currentIteration );


#endif  // TUDATBUNDLE_APPLYDIFFERENTIALCORRECTION_H
