#ifndef TUDATBUNDLE_PROPAGATEORBIT_H
#define TUDATBUNDLE_PROPAGATEORBIT_H

#include <map>

#include <Eigen/Core>

#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"

#include "cr3bpPeriodicOrbitTypes.h"

Eigen::MatrixXd getFullInitialState( const Eigen::Vector6d& initialState );

void writeStateHistoryToFile(
        const std::map< double, Eigen::Vector6d >& stateHistory,
        const int orbitId, const tudat::cr3bp::CR3BPPeriodicOrbitTypes orbitType, const int librationPointNr,
        const int saveEveryNthIntegrationStep, const bool completeInitialConditionsHaloFamily );

std::pair< Eigen::MatrixXd, double > propagateOrbit(
        const Eigen::MatrixXd& stateVectorInclSTM, double massParameter, double currentTime,
        int direction, double initialStepSize = 1.0E-5, double maximumStepSize = 1.0E-4 );

std::pair< Eigen::Matrix< double, 6, 7 >, double >  propagateOrbitToFinalCondition(
        const Eigen::Matrix< double, 6, 7 >& fullInitialState,
        boost::function< Eigen::Matrix< double, 6, 7 >( const double, const Eigen::Matrix< double, 6, 7 >& ) > stateDerivativeFunction,
        const boost::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const double finalTime, std::map< double, Eigen::Vector6d >& stateHistory );

std::pair< Eigen::MatrixXd, double >  propagateOrbitWithStateTransitionMatrixToFinalCondition(
        const Eigen::MatrixXd fullInitialState, const double massParameter, const double finalTime, int direction,
        std::map< double, Eigen::MatrixXd >& stateTransitionMatrixHistory, const int saveFrequency = -1, const double initialTime = 0.0);

#endif  // TUDATBUNDLE_PROPAGATEORBIT_H
