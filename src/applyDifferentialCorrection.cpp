#include <cmath>

#include "Tudat/Astrodynamics/BasicAstrodynamics/celestialBodyConstants.h"
#include "Tudat/Astrodynamics/Gravitation/librationPoint.h"
#include "Tudat/Astrodynamics/Gravitation/jacobiEnergy.h"
#include "Tudat/Astrodynamics/Propagators/integrateEquations.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationTerminationSettings.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationTermination.h"

#include "applyDifferentialCorrection.h"
#include "computeDifferentialCorrection.h"
#include "propagateOrbit.h"
#include "stateDerivativeModel.h"
#include "cr3bpPeriodicOrbits.h"

using namespace tudat;

Eigen::VectorXd applyDifferentialCorrection(
        const Eigen::VectorXd& initialStateVector,
        double orbitalPeriod,
        const boost::shared_ptr< tudat::cr3bp::CR3BPPeriodicOrbitModel > periodicOrbitModel,
        const boost::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const int currentIteration )
{
    Eigen::MatrixXd initialStateVectorInclSTM = Eigen::MatrixXd::Zero( 6, 7 );

    initialStateVectorInclSTM.block( 0, 0, 6, 1 ) = initialStateVector;
    initialStateVectorInclSTM.block( 0, 1, 6, 6 ).setIdentity( );

    std::map< double, Eigen::Vector6d > stateHistory;
    std::pair< Eigen::MatrixXd, double > halfPeriodState = propagateOrbitToFinalCondition(
                initialStateVectorInclSTM,
                periodicOrbitModel->getStateDerivativeFunctionWithStateTransition( ),
                integratorSettings, orbitalPeriod / 2.0, stateHistory );

    Eigen::MatrixXd stateVectorInclSTM = halfPeriodState.first;
    double currentTime = halfPeriodState.second;
    Eigen::VectorXd stateVectorOnly = stateVectorInclSTM.block( 0, 0, 6, 1 );

    // Initialize variables
    Eigen::VectorXd differentialCorrection(7);
    Eigen::VectorXd outputVector(15);

    int numberOfIterations = 0;
    // Apply differential correction and propagate to half-period point until converged.
    while ( periodicOrbitModel->continueDifferentialCorrection( stateVectorOnly, numberOfIterations ) )
    {

        differentialCorrection = periodicOrbitModel->computeDifferentialCorrection(
                   stateVectorInclSTM.block( 0, 1, 6, 6 ), stateVectorOnly, currentTime, numberOfIterations );
        initialStateVectorInclSTM.block( 0, 0, 6, 1 ) += differentialCorrection.segment( 0, 6 ) / 1.0;        

        orbitalPeriod  = orbitalPeriod + 2.0 * differentialCorrection( 6 ) / 1.0;

        std::pair< Eigen::MatrixXd, double > halfPeriodState = propagateOrbitToFinalCondition(
                    initialStateVectorInclSTM, periodicOrbitModel->getStateDerivativeFunctionWithStateTransition( ),
                    integratorSettings, orbitalPeriod / 2.0, stateHistory );

        stateVectorInclSTM = halfPeriodState.first;
        currentTime = halfPeriodState.second;
        stateVectorOnly = stateVectorInclSTM.block( 0, 0, 6, 1 );

        numberOfIterations += 1;
    }
//    std::cout<<"Number of iterations: "<<numberOfIterations<<std::endl;

//    double jacobiEnergyHalfPeriod       = tudat::gravitation::computeJacobiEnergy(massParameter, stateVectorOnly);
//    double jacobiEnergyInitialCondition = tudat::gravitation::computeJacobiEnergy(massParameter, initialStateVectorInclSTM.block( 0, 0, 6, 1 ));

//    std::cout << "\nCorrected initial state vector:" << std::endl << initialStateVectorInclSTM.block( 0, 0, 6, 1 )        << std::endl
//              << "\nwith orbital period: "           << orbitalPeriod                                              << std::endl
//              << "||J(0) - J(T/2|| = "               << std::abs(jacobiEnergyInitialCondition - jacobiEnergyHalfPeriod) << std::endl
//              << "||T/2 - t|| = "                    << std::abs(orbitalPeriod/2.0 - currentTime) << "\n"               << std::endl;

    // The output vector consists of:
    // 1. Corrected initial state vector, including orbital period
    // 2. Half period state vector, including currentTime of integration
    // 3. numberOfIterations
    outputVector.segment(0,6)    = initialStateVectorInclSTM.block( 0, 0, 6, 1 );
    outputVector(6)              = orbitalPeriod;
    outputVector.segment(7,6)    = stateVectorOnly;
    outputVector(13)             = currentTime;
    outputVector(14)             = numberOfIterations;

    return outputVector;
}
