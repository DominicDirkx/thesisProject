#include <cmath>

#include "Tudat/Astrodynamics/BasicAstrodynamics/celestialBodyConstants.h"
#include "Tudat/Astrodynamics/Gravitation/librationPoint.h"
#include "Tudat/Astrodynamics/Gravitation/jacobiEnergy.h"
#include "Tudat/Astrodynamics/Propagators/integrateEquations.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationTerminationSettings.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationTermination.h"
#include "Tudat/Mathematics/NumericalIntegrators/rungeKuttaVariableStepSizeIntegrator.h"

#include "applyDifferentialCorrection.h"
#include "computeDifferentialCorrection.h"
#include "propagateOrbit.h"
#include "stateDerivativeModel.h"

using namespace tudat;

Eigen::VectorXd applyDifferentialCorrection(const int librationPointNr, const std::string& orbitType,
                                            const Eigen::VectorXd& initialStateVector,
                                            double orbitalPeriod, const double massParameter,
                                            double maxPositionDeviationFromPeriodicOrbit,
                                            double maxVelocityDeviationFromPeriodicOrbit,
                                            const int maxNumberOfIterations )
{
    std::cout << "\nApply differential correction:" << std::endl;

    Eigen::MatrixXd initialStateVectorInclSTM = Eigen::MatrixXd::Zero( 6, 7 );

    initialStateVectorInclSTM.block( 0, 0, 6, 1 ) = initialStateVector;
    initialStateVectorInclSTM.block( 0, 1, 6, 6 ).setIdentity( );

    std::map< double, Eigen::Vector6d > stateHistory;

    clock_t tstart2, tend2;
    tstart2 = clock ();
    std::pair< Eigen::MatrixXd, double > halfPeriodState = propagateOrbitToFinalCondition(
                initialStateVectorInclSTM, massParameter, orbitalPeriod / 2.0, 1.0, stateHistory, 1, 0.0 );
    tend2 = clock ();

    std::cout<<"Runtime with original scheme "<<tend2 -tstart2<<std::endl;

    std::map< double, Eigen::MatrixXd > halfPeriodState2;
    {

        std::map< double, double > cummulativeComputationTimeHistory;
        std::map< double, Eigen::VectorXd > dependentVariableHistory;
        boost::function< Eigen::VectorXd( ) > dependentVariableFunction;
        int saveFrequency = 1;

        double minimumStepSize   = std::numeric_limits<double>::epsilon( ); // 2.22044604925031e-16
        const double relativeErrorTolerance = 100.0 * std::numeric_limits<double>::epsilon( ); // 2.22044604925031e-14
        const double absoluteErrorTolerance = 1.0e-24;

        // Create integrator to be used for propagating.
        boost::shared_ptr< numerical_integrators::NumericalIntegrator< double, Eigen::MatrixXd > > orbitIntegrator = boost::make_shared<
                tudat::numerical_integrators::RungeKuttaVariableStepSizeIntegrator< double, Eigen::MatrixXd > >(
                    tudat::numerical_integrators::RungeKuttaCoefficients::get( tudat::numerical_integrators::RungeKuttaCoefficients::rungeKuttaFehlberg78 ),
                    &computeStateDerivative, 0.0, initialStateVectorInclSTM, minimumStepSize, 1.0E-5, relativeErrorTolerance, absoluteErrorTolerance );

        double initialTimeStep = 1.0E-5;
        boost::shared_ptr< propagators::PropagationTerminationCondition > propagationTerminationCondition =
                propagators::createPropagationTerminationConditions(
                    boost::make_shared< propagators::PropagationTimeTerminationSettings >(
                        orbitalPeriod / 2.0, true ), simulation_setup::NamedBodyMap( ), 1.0E-5 );

        clock_t tstart, tend;
        tstart = clock( );
        propagators::integrateEquationsFromIntegrator< Eigen::MatrixXd, double >(
                    orbitIntegrator, initialTimeStep, propagationTerminationCondition, halfPeriodState2,
                    dependentVariableHistory, cummulativeComputationTimeHistory, dependentVariableFunction, saveFrequency );
        tend = clock( );

        std::cout<<"Runtime with Tudat scheme "<<tend -tstart<<std::endl;

    }
    std::cout<<stateHistory.rbegin( )->second.transpose( )<<std::endl;
    std::cout<<halfPeriodState2.rbegin( )->second.transpose( )<<std::endl;


    Eigen::MatrixXd stateVectorInclSTM      = halfPeriodState.first;
    double currentTime             = halfPeriodState.second;
    Eigen::VectorXd stateVectorOnly = stateVectorInclSTM.block( 0, 0, 6, 1 );

    // Initialize variables
    Eigen::VectorXd differentialCorrection(7);
    Eigen::VectorXd outputVector(15);
    double positionDeviationFromPeriodicOrbit;
    double velocityDeviationFromPeriodicOrbit;

    if (orbitType == "axial")
    {
        // Initial condition for axial family should be [x, 0, 0, 0, ydot, zdot]
        positionDeviationFromPeriodicOrbit = sqrt(pow(stateVectorOnly(1), 2) + pow(stateVectorOnly(2), 2));
        velocityDeviationFromPeriodicOrbit = sqrt(pow(stateVectorOnly(3), 2));
    }
    else
    {
        // Initial condition for other families should be [x, 0, y, 0, ydot, 0]
        positionDeviationFromPeriodicOrbit = sqrt(pow(stateVectorOnly(1), 2));
        velocityDeviationFromPeriodicOrbit = sqrt(pow(stateVectorOnly(3), 2) + pow(stateVectorOnly(5), 2));
    }

    std::cout << "\nInitial state vector:\n"                  << initialStateVectorInclSTM.block( 0, 0, 6, 1 )
              << "\nPosition deviation from periodic orbit: " << positionDeviationFromPeriodicOrbit
              << "\nVelocity deviation from periodic orbit: " << velocityDeviationFromPeriodicOrbit
              << "\n\nDifferential correction:"               << std::endl;

    bool deviationFromPeriodicOrbitRelaxed = false;

    int numberOfIterations = 0;
    // Apply differential correction and propagate to half-period point until converged.
    while ( positionDeviationFromPeriodicOrbit > maxPositionDeviationFromPeriodicOrbit or
            velocityDeviationFromPeriodicOrbit > maxVelocityDeviationFromPeriodicOrbit )
    {

        // If the maximum number of iterations has been reached, return a zero vector to stop the numerical continuation
        if ( numberOfIterations > maxNumberOfIterations and deviationFromPeriodicOrbitRelaxed == false )
        {
            // Relax the periodicity constraints after exceeding the maximum number of iterations instead of termination
            maxPositionDeviationFromPeriodicOrbit = 10.0 * maxPositionDeviationFromPeriodicOrbit;
            maxVelocityDeviationFromPeriodicOrbit = 10.0 * maxVelocityDeviationFromPeriodicOrbit;
            deviationFromPeriodicOrbitRelaxed = true;
//            return Eigen::VectorXd::Zero(15);
        }

        // Relax the maximum deviation requirements to compute the horizontal Lyapunov family in L2
        if (deviationFromPeriodicOrbitRelaxed == false and numberOfIterations > 10 and
                orbitType == "horizontal" and librationPointNr == 2)
        {

            maxPositionDeviationFromPeriodicOrbit = 10.0 * maxPositionDeviationFromPeriodicOrbit;
            maxVelocityDeviationFromPeriodicOrbit = 10.0 * maxVelocityDeviationFromPeriodicOrbit;
            deviationFromPeriodicOrbitRelaxed = true;
        }

        // Apply differential correction
        if (numberOfIterations > 10 and orbitType == "axial" and librationPointNr == 2)
        {
            // To compute the full L2 axial family, fix x position after not finding a fully periodic solution after 10 iterations
            differentialCorrection = computeDifferentialCorrection( librationPointNr, orbitType, stateVectorInclSTM, true );
        }
        else
        {
            differentialCorrection = computeDifferentialCorrection( librationPointNr, orbitType, stateVectorInclSTM, false );
        }

        std::cout<<"APPLYING DIFF. CORR: "<<differentialCorrection<<std::endl;

        initialStateVectorInclSTM.block( 0, 0, 6, 1 ) += differentialCorrection.segment( 0, 6 ) / 1.0;
        orbitalPeriod  = orbitalPeriod + 2.0 * differentialCorrection( 6 ) / 1.0;

        std::pair< Eigen::MatrixXd, double > halfPeriodState = propagateOrbitToFinalCondition(
                    initialStateVectorInclSTM, massParameter, orbitalPeriod / 2.0, 1.0, stateHistory, -1, 0.0 );
        stateVectorInclSTM      = halfPeriodState.first;
        currentTime             = halfPeriodState.second;
        stateVectorOnly = stateVectorInclSTM.block( 0, 0, 6, 1 );

        if (orbitType == "axial")
        {
            // Initial condition for axial family should be [x, 0, 0, 0, ydot, zdot]
            positionDeviationFromPeriodicOrbit = sqrt(pow(stateVectorOnly(1), 2) + pow(stateVectorOnly(2), 2));
            velocityDeviationFromPeriodicOrbit = sqrt(pow(stateVectorOnly(3), 2));
        }
        else
        {
            // Initial condition for other families should be [x, 0, y, 0, ydot, 0]
            positionDeviationFromPeriodicOrbit = sqrt(pow(stateVectorOnly(1), 2));
            velocityDeviationFromPeriodicOrbit = sqrt(pow(stateVectorOnly(3), 2) + pow(stateVectorOnly(5), 2));
        }
        numberOfIterations += 1;

        std::cout << "positionDeviationFromPeriodicOrbit: " << positionDeviationFromPeriodicOrbit << std::endl
                  << "velocityDeviationFromPeriodicOrbit: " << velocityDeviationFromPeriodicOrbit << "\n" << std::endl;
    }

    double jacobiEnergyHalfPeriod       = tudat::gravitation::computeJacobiEnergy(massParameter, stateVectorOnly);
    double jacobiEnergyInitialCondition = tudat::gravitation::computeJacobiEnergy(massParameter, initialStateVectorInclSTM.block( 0, 0, 6, 1 ));

    std::cout << "\nCorrected initial state vector:" << std::endl << initialStateVectorInclSTM.block( 0, 0, 6, 1 )        << std::endl
              << "\nwith orbital period: "           << orbitalPeriod                                              << std::endl
              << "||J(0) - J(T/2|| = "               << std::abs(jacobiEnergyInitialCondition - jacobiEnergyHalfPeriod) << std::endl
              << "||T/2 - t|| = "                    << std::abs(orbitalPeriod/2.0 - currentTime) << "\n"               << std::endl;

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
