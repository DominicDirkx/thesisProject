#include <cmath>

#include "Tudat/Astrodynamics/BasicAstrodynamics/celestialBodyConstants.h"
#include "Tudat/Astrodynamics/Gravitation/librationPoint.h"
#include "Tudat/Astrodynamics/Gravitation/jacobiEnergy.h"
#include "Tudat/Astrodynamics/Propagators/integrateEquations.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationTerminationSettings.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationTermination.h"

#include "applyDifferentialCorrection.h"
#include "propagateOrbit.h"
#include "cr3bpPeriodicOrbits.h"

using namespace tudat;

Eigen::VectorXd applyDifferentialCorrectionForPeriodicOrbit(
        const Eigen::VectorXd& initialStateVectorGuess,
        double orbitalPeriod,
        const boost::shared_ptr< tudat::cr3bp::CR3BPPeriodicOrbitModel > periodicOrbitModel,
        const boost::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const int currentIteration )
{
    std::map< double, Eigen::Vector6d > stateHistory;

    // Propagate initial guess
    Eigen::Matrix< double, 6, 7 > initialStateWithSTM = Eigen::Matrix< double, 6, 7 >::Zero( );
    initialStateWithSTM.block( 0, 0, 6, 1 ) = initialStateVectorGuess;
    initialStateWithSTM.block( 0, 1, 6, 6 ).setIdentity( );
    std::pair< Eigen::Matrix< double, 6, 7 >, double > halfPeriodState = propagateOrbitToFinalCondition(
                initialStateWithSTM,
                periodicOrbitModel->getStateDerivativeFunctionWithStateTransition( ),
                integratorSettings, orbitalPeriod / 2.0, stateHistory );

    // Retrieve state and time of propagation
    Eigen::Matrix< double, 6, 7 > stateWithSTMAtHalfPeriod = halfPeriodState.first;
    double timeAtHalfPeriod = halfPeriodState.second;
    Eigen::VectorXd stateVectorAtHalfPeriod = stateWithSTMAtHalfPeriod.block( 0, 0, 6, 1 );

    // Initialize variables
    Eigen::VectorXd differentialCorrection(7);
    Eigen::VectorXd outputVector(15);

    int numberOfIterations = 0;

    // Apply differential correction and propagate to half-period point until converged.
    do
    {
        // Compute differential correction
        differentialCorrection = periodicOrbitModel->computeDifferentialCorrection(
                   stateWithSTMAtHalfPeriod.block( 0, 1, 6, 6 ), stateVectorAtHalfPeriod, timeAtHalfPeriod, numberOfIterations );

        // Update initial state and time
        initialStateWithSTM.block( 0, 0, 6, 1 ) += differentialCorrection.segment( 0, 6 );
        orbitalPeriod  = orbitalPeriod + 2.0 * differentialCorrection( 6 );

        // Propagate corrected state and retrieve results
        std::pair< Eigen::Matrix< double, 6, 7 >, double > halfPeriodState = propagateOrbitToFinalCondition(
                    initialStateWithSTM, periodicOrbitModel->getStateDerivativeFunctionWithStateTransition( ),
                    integratorSettings, orbitalPeriod / 2.0, stateHistory );
        stateWithSTMAtHalfPeriod = halfPeriodState.first;
        timeAtHalfPeriod = halfPeriodState.second;
        stateVectorAtHalfPeriod = stateWithSTMAtHalfPeriod.block( 0, 0, 6, 1 );

        numberOfIterations += 1;
    }
    while ( periodicOrbitModel->continueDifferentialCorrection( stateVectorAtHalfPeriod, numberOfIterations ) );

//    double jacobiEnergyHalfPeriod       = tudat::gravitation::computeJacobiEnergy(massParameter, stateVectorAtHalfPeriod);
//    double jacobiEnergyInitialCondition = tudat::gravitation::computeJacobiEnergy(massParameter, initialStateWithSTM.block( 0, 0, 6, 1 ));

//    std::cout << "\nCorrected initial state vector:" << std::endl << initialStateWithSTM.block( 0, 0, 6, 1 )        << std::endl
//              << "\nwith orbital period: "           << orbitalPeriod                                              << std::endl
//              << "||J(0) - J(T/2|| = "               << std::abs(jacobiEnergyInitialCondition - jacobiEnergyHalfPeriod) << std::endl
//              << "||T/2 - t|| = "                    << std::abs(orbitalPeriod/2.0 - timeAtHalfPeriod) << "\n"               << std::endl;

    // The output vector consists of:
    // 1. Corrected initial state vector, including orbital period
    // 2. Half period state vector, including timeAtHalfPeriod of integration
    // 3. numberOfIterations
    outputVector.segment(0,6)    = initialStateWithSTM.block( 0, 0, 6, 1 );
    outputVector(6)              = orbitalPeriod;
    outputVector.segment(7,6)    = stateVectorAtHalfPeriod;
    outputVector(13)             = timeAtHalfPeriod;
    outputVector(14)             = numberOfIterations;

    std::cout<<orbitalPeriod - 2.0 * timeAtHalfPeriod<<" "<<timeAtHalfPeriod<<std::endl;

    return outputVector;
}
