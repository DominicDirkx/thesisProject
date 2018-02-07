#ifndef TUDATBUNDLE_CREATEINITIALCONDITIONS_H
#define TUDATBUNDLE_CREATEINITIALCONDITIONS_H



#include <string>
#include <vector>
#include <map>
#include <iostream>

#include <Eigen/Core>

#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"

#include "cr3bpPeriodicOrbits.h"

void appendResultsVector(
        const double jacobiEnergy, const double orbitalPeriod, const Eigen::VectorXd& initialStateVector,
        const Eigen::MatrixXd& stateVectorInclSTM, std::vector< Eigen::VectorXd >& initialConditions );

void appendDifferentialCorrectionResultsVector(
        const double jacobiEnergyHalfPeriod,  const Eigen::VectorXd& differentialCorrectionResult,
        std::vector< Eigen::VectorXd >& differentialCorrections );

void writeFinalResultsToFiles( const int librationPointNr, const tudat::cr3bp::CR3BPPeriodicOrbitTypes orbitType,
                               std::vector< Eigen::VectorXd > initialConditions,
                               std::vector< Eigen::VectorXd > differentialCorrections );




//! Apply differential correction, and save results for periodic orbit
Eigen::MatrixXd correctPeriodicOrbitInitialState(
                const Eigen::Vector6d& initialStateGuess, double orbitalPeriod,
                const int orbitNumber,
                const boost::shared_ptr< tudat::cr3bp::CR3BPPeriodicOrbitModel > periodicOrbitModel,
                const boost::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > > integratorSettings,
                std::vector< Eigen::VectorXd >& initialConditions,
                std::vector< Eigen::VectorXd >& differentialCorrections );

double getDefaultArcLength(
        const double distanceIncrement,
        const Eigen::Vector6d& currentState );

void createPeriodicOrbitInitialConditionsFromExistingData(
        const boost::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const boost::shared_ptr< tudat::cr3bp::CR3BPPeriodicOrbitModel > periodicOrbitModel,
        std::vector< Eigen::VectorXd >& initialConditions,
        std::vector< Eigen::VectorXd >& differentialCorrections,
        const boost::function< double( const Eigen::Vector6d& ) > pseudoArcLengthFunction );

void createPeriodicOrbitInitialConditions(
        const boost::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const boost::shared_ptr< tudat::cr3bp::CR3BPPeriodicOrbitModel > periodicOrbitModel,
        const boost::function< double( const Eigen::Vector6d& ) > pseudoArcLengthFunction =
        boost::bind( &getDefaultArcLength, 1.0E-4, _1 ) );


#endif  // TUDATBUNDLE_CREATEINITIALCONDITIONS_H
