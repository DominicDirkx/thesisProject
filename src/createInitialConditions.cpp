#include <fstream>
#include <iomanip>
#include <chrono>
#include <ctime>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>
#include <Eigen/Dense>

#include "Tudat/Astrodynamics/BasicAstrodynamics/celestialBodyConstants.h"
#include "Tudat/Astrodynamics/Gravitation/librationPoint.h"
#include "Tudat/Astrodynamics/Gravitation/jacobiEnergy.h"
#include "Tudat/Astrodynamics/Propagators/integrateEquations.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationTerminationSettings.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationTermination.h"
#include "Tudat/Mathematics/NumericalIntegrators/rungeKuttaVariableStepSizeIntegrator.h"

#include "createInitialConditions.h"
#include "applyDifferentialCorrection.h"
#include "propagateOrbit.h"
#include "richardsonThirdOrderApproximation.h"

using namespace tudat;

void appendResultsVector(const double jacobiEnergy, const double orbitalPeriod, const Eigen::VectorXd& initialStateVector,
                         const Eigen::MatrixXd& stateVectorInclSTM, std::vector< Eigen::VectorXd >& initialConditions )
{
    Eigen::VectorXd tempInitialCondition = Eigen::VectorXd( 44 );

    // Add Jacobi energy and orbital period
    tempInitialCondition( 0 ) = jacobiEnergy;
    tempInitialCondition( 1 ) = orbitalPeriod;

    // Add initial condition of periodic solution
    for (int i = 0; i <= 5; i++){
        tempInitialCondition( i + 2 ) = initialStateVector( i );
    }

    // Add Monodromy matrix
    tempInitialCondition( 8 ) = stateVectorInclSTM(0, 1);
    tempInitialCondition( 9 ) = stateVectorInclSTM(0, 2);
    tempInitialCondition( 10 ) = stateVectorInclSTM(0, 3);
    tempInitialCondition( 11 ) = stateVectorInclSTM(0, 4);
    tempInitialCondition( 12 ) = stateVectorInclSTM(0, 5);
    tempInitialCondition( 13 ) = stateVectorInclSTM(0, 6);

    tempInitialCondition( 14 ) = stateVectorInclSTM(1, 1);
    tempInitialCondition( 15 ) = stateVectorInclSTM(1, 2);
    tempInitialCondition( 16 ) = stateVectorInclSTM(1, 3);
    tempInitialCondition( 17 ) = stateVectorInclSTM(1, 4);
    tempInitialCondition( 18 ) = stateVectorInclSTM(1, 5);
    tempInitialCondition( 19 ) = stateVectorInclSTM(1, 6);

    tempInitialCondition( 20 ) = stateVectorInclSTM(2, 1);
    tempInitialCondition( 21 ) = stateVectorInclSTM(2, 2);
    tempInitialCondition( 22 ) = stateVectorInclSTM(2, 3);
    tempInitialCondition( 23 ) = stateVectorInclSTM(2, 4);
    tempInitialCondition( 24 ) = stateVectorInclSTM(2, 5);
    tempInitialCondition( 25 ) = stateVectorInclSTM(2, 6);

    tempInitialCondition( 26 ) = stateVectorInclSTM(3, 1);
    tempInitialCondition( 27 ) = stateVectorInclSTM(3, 2);
    tempInitialCondition( 28 ) = stateVectorInclSTM(3, 3);
    tempInitialCondition( 29 ) = stateVectorInclSTM(3, 4);
    tempInitialCondition( 30 ) = stateVectorInclSTM(3, 5);
    tempInitialCondition( 31 ) = stateVectorInclSTM(3, 6);

    tempInitialCondition( 32 ) = stateVectorInclSTM(4, 1);
    tempInitialCondition( 33 ) = stateVectorInclSTM(4, 2);
    tempInitialCondition( 34 ) = stateVectorInclSTM(4, 3);
    tempInitialCondition( 35 ) = stateVectorInclSTM(4, 4);
    tempInitialCondition( 36 ) = stateVectorInclSTM(4, 5);
    tempInitialCondition( 37 ) = stateVectorInclSTM(4, 6);

    tempInitialCondition( 38 ) = stateVectorInclSTM(5, 1);
    tempInitialCondition( 39 ) = stateVectorInclSTM(5, 2);
    tempInitialCondition( 40 ) = stateVectorInclSTM(5, 3);
    tempInitialCondition( 41 ) = stateVectorInclSTM(5, 4);
    tempInitialCondition( 42 ) = stateVectorInclSTM(5, 5);
    tempInitialCondition( 43 ) = stateVectorInclSTM(5, 6);

    initialConditions.push_back(tempInitialCondition);
}

void appendDifferentialCorrectionResultsVector(
        const double jacobiEnergyHalfPeriod,  const Eigen::VectorXd& differentialCorrectionResult,
        std::vector< Eigen::VectorXd >& differentialCorrections )
{

    Eigen::VectorXd tempDifferentialCorrection = Eigen::VectorXd( 9 );
    tempDifferentialCorrection( 0 ) = differentialCorrectionResult( 14 );  // numberOfIterations
    tempDifferentialCorrection( 1 ) = jacobiEnergyHalfPeriod;  // jacobiEnergyHalfPeriod
    tempDifferentialCorrection( 2 ) = differentialCorrectionResult( 13 );  // currentTime
    for (int i = 7; i <= 12; i++)
    {
        tempDifferentialCorrection( i - 4 ) = differentialCorrectionResult( i );  // halfPeriodStateVector
    }

    
    differentialCorrections.push_back(tempDifferentialCorrection);

}

Eigen::MatrixXd correctPeriodicOrbitInitialState(
        const Eigen::Vector6d& initialStateGuess, double orbitalPeriod,
        const int orbitNumber,
        const boost::shared_ptr< tudat::cr3bp::CR3BPPeriodicOrbitModel > periodicOrbitModel,
        const boost::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > > integratorSettings,
        std::vector< Eigen::VectorXd >& initialConditions,
        std::vector< Eigen::VectorXd >& differentialCorrections )
{
    Eigen::Vector6d initialStateVector = initialStateGuess;


    // Correct state vector guess
    Eigen::VectorXd differentialCorrectionResult = applyDifferentialCorrectionForPeriodicOrbit(
                initialStateVector, orbitalPeriod, periodicOrbitModel, integratorSettings, orbitNumber );
    initialStateVector = differentialCorrectionResult.segment( 0, 6 );
    orbitalPeriod = differentialCorrectionResult( 6 );

    if( orbitalPeriod > 0 )
    {
        // Propagate the initialStateVector for a full period and write output to file.
        std::map< double, Eigen::Vector6d > stateHistory;
        Eigen::MatrixXd stateVectorInclSTM = propagateOrbitToFinalCondition(
                    getFullInitialState( initialStateVector ),
                    periodicOrbitModel->getStateDerivativeFunctionWithStateTransition( ),
                    integratorSettings, orbitalPeriod, stateHistory ).first;

        periodicOrbitModel->addPeriodicOrbit(
                    initialStateVector, stateVectorInclSTM.block( 0, 1, 6, 6 ),
                   orbitalPeriod, stateVectorInclSTM.block( 0, 0, 6, 1 ), differentialCorrectionResult( 14 ) , stateHistory );

//        writeStateHistoryToFile( stateHistory, orbitNumber, periodicOrbitModel->getOrbitType( ),
//                                 periodicOrbitModel->getLibrationPointNumber( ), 1000, false );

////        // Save results
        double jacobiEnergyHalfPeriod = tudat::gravitation::computeJacobiEnergy(
                    periodicOrbitModel->getMassParameter( ), differentialCorrectionResult.segment( 7, 6 ) );
        appendDifferentialCorrectionResultsVector( jacobiEnergyHalfPeriod, differentialCorrectionResult, differentialCorrections );

        double jacobiEnergy = tudat::gravitation::computeJacobiEnergy(
                    periodicOrbitModel->getMassParameter( ), stateVectorInclSTM.block( 0, 0, 6, 1 ));
        appendResultsVector( jacobiEnergy, orbitalPeriod, initialStateVector, stateVectorInclSTM, initialConditions );

        return stateVectorInclSTM;
    }
    else
    {
        return Eigen::MatrixXd::Constant( 6, 7, TUDAT_NAN );
    }


}

void writeFinalResultsToFiles( const int librationPointNr, const tudat::cr3bp::CR3BPPeriodicOrbitTypes orbitType,
                               std::vector< Eigen::VectorXd > initialConditions,
                               std::vector< Eigen::VectorXd > differentialCorrections )
{
    // Prepare file for initial conditions
    remove(("../data/raw/orbits/L" + std::to_string(librationPointNr) + "_" + std::to_string( orbitType ) + "_initial_conditions.txt").c_str());
    std::ofstream textFileInitialConditions;
    textFileInitialConditions.open(("../data/raw/orbits/L" + std::to_string(librationPointNr) + "_" + std::to_string( orbitType ) + "_initial_conditions.txt"));
    textFileInitialConditions.precision(std::numeric_limits<double>::digits10);

    // Prepare file for differential correction
    remove(("../data/raw/orbits/L" + std::to_string(librationPointNr) + "_" + std::to_string( orbitType ) + "_differential_correction.txt").c_str());
    std::ofstream textFileDifferentialCorrection;
    textFileDifferentialCorrection.open(("../data/raw/orbits/L" + std::to_string(librationPointNr) + "_" + std::to_string( orbitType ) + "_differential_correction.txt"));
    textFileDifferentialCorrection.precision(std::numeric_limits<double>::digits10);

    // Write initial conditions to file
    for (unsigned int i=0; i<initialConditions.size(); i++) {

        textFileInitialConditions << std::left << std::scientific                                          << std::setw(25)
                                  << initialConditions[i][0]  << std::setw(25) << initialConditions[i][1]  << std::setw(25)
                                  << initialConditions[i][2]  << std::setw(25) << initialConditions[i][3]  << std::setw(25)
                                  << initialConditions[i][4]  << std::setw(25) << initialConditions[i][5]  << std::setw(25)
                                  << initialConditions[i][6]  << std::setw(25) << initialConditions[i][7]  << std::setw(25)
                                  << initialConditions[i][8]  << std::setw(25) << initialConditions[i][9]  << std::setw(25)
                                  << initialConditions[i][10] << std::setw(25) << initialConditions[i][11] << std::setw(25)
                                  << initialConditions[i][12] << std::setw(25) << initialConditions[i][13] << std::setw(25)
                                  << initialConditions[i][14] << std::setw(25) << initialConditions[i][15] << std::setw(25)
                                  << initialConditions[i][16] << std::setw(25) << initialConditions[i][17] << std::setw(25)
                                  << initialConditions[i][18] << std::setw(25) << initialConditions[i][19] << std::setw(25)
                                  << initialConditions[i][20] << std::setw(25) << initialConditions[i][21] << std::setw(25)
                                  << initialConditions[i][22] << std::setw(25) << initialConditions[i][23] << std::setw(25)
                                  << initialConditions[i][24] << std::setw(25) << initialConditions[i][25] << std::setw(25)
                                  << initialConditions[i][26] << std::setw(25) << initialConditions[i][27] << std::setw(25)
                                  << initialConditions[i][28] << std::setw(25) << initialConditions[i][29] << std::setw(25)
                                  << initialConditions[i][30] << std::setw(25) << initialConditions[i][31] << std::setw(25)
                                  << initialConditions[i][32] << std::setw(25) << initialConditions[i][33] << std::setw(25)
                                  << initialConditions[i][34] << std::setw(25) << initialConditions[i][35] << std::setw(25)
                                  << initialConditions[i][36] << std::setw(25) << initialConditions[i][37] << std::setw(25)
                                  << initialConditions[i][38] << std::setw(25) << initialConditions[i][39] << std::setw(25)
                                  << initialConditions[i][40] << std::setw(25) << initialConditions[i][41] << std::setw(25)
                                  << initialConditions[i][42] << std::setw(25) << initialConditions[i][43] << std::endl;

        textFileDifferentialCorrection << std::left << std::scientific   << std::setw(25)
                                       << differentialCorrections[i][0]  << std::setw(25) << differentialCorrections[i][1 ]  << std::setw(25)
                                       << differentialCorrections[i][2 ]  << std::setw(25) << differentialCorrections[i][3]  << std::setw(25)
                                       << differentialCorrections[i][4]  << std::setw(25) << differentialCorrections[i][5]  << std::setw(25)
                                       << differentialCorrections[i][6]  << std::setw(25) << differentialCorrections[i][7]  << std::setw(25)
                                       << differentialCorrections[i][8]  << std::setw(25) << std::endl;
    }
}


double getDefaultArcLength(
        const double distanceIncrement,
        const Eigen::Vector6d& currentState )
{
    return std::max( distanceIncrement / currentState.segment( 0, 3 ).norm( ), 0.01 );
}

void createPeriodicOrbitInitialConditionsFromExistingData(
        const boost::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const boost::shared_ptr< tudat::cr3bp::CR3BPPeriodicOrbitModel > periodicOrbitModel,
        std::vector< Eigen::VectorXd >& initialConditions,
        std::vector< Eigen::VectorXd >& differentialCorrections,
        const boost::function< double( const Eigen::Vector6d& ) > pseudoArcLengthFunction )

{
    // Set exit parameters of continuation procedure
    int numberOfInitialConditions = 2;
    int maximumNumberOfInitialConditions = 10000;

    // Initialize state vectors and orbital periods
    Eigen::Vector6d initialStateVector = Eigen::Vector6d::Zero( );
    Eigen::MatrixXd stateVectorInclSTM = Eigen::MatrixXd::Zero( 6, 7 );

    // Generate periodic orbits until termination
    double orbitalPeriod  = 0.0, periodIncrement = 0.0, pseudoArcLengthCorrection = 0.0;
    bool continueNumericalContinuation = true;
    Eigen::Vector6d stateIncrement;

    std::clock_t c_start = std::clock();
    while( ( numberOfInitialConditions < maximumNumberOfInitialConditions ) && continueNumericalContinuation)
    {
        // Determine increments to state and time
        stateIncrement = initialConditions[ initialConditions.size( ) - 1 ].segment( 2, 6 ) -
                initialConditions[ initialConditions.size( ) - 2 ].segment( 2, 6 );
        periodIncrement = initialConditions[ initialConditions.size( ) - 1 ]( 1 ) -
                initialConditions[ initialConditions.size( ) - 2 ]( 1 );
        pseudoArcLengthCorrection =
                pseudoArcLengthFunction( stateIncrement );

        // Apply numerical continuation
        initialStateVector = initialConditions[ initialConditions.size( ) - 1 ].segment( 2, 6 ) +
                stateIncrement * pseudoArcLengthCorrection;
        orbitalPeriod = initialConditions[ initialConditions.size( ) - 1 ]( 1 ) +
                periodIncrement * pseudoArcLengthCorrection;
        stateVectorInclSTM = correctPeriodicOrbitInitialState(
                    initialStateVector, orbitalPeriod, numberOfInitialConditions, periodicOrbitModel,
                    integratorSettings, initialConditions, differentialCorrections );

        continueNumericalContinuation =
                periodicOrbitModel->terminateNumericalContinuation(
                    stateVectorInclSTM.block( 0, 1, 6, 6 ), stateVectorInclSTM.block( 0, 0, 6, 1 ), orbitalPeriod, numberOfInitialConditions );

        numberOfInitialConditions += 1;
    }
    std::clock_t c_end = std::clock();

    std::cout<<"Generating initial state, "<<numberOfInitialConditions<<" orbits at. Total generation time = "<< 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC <<" milliseconds "<<std::endl;

    writeFinalResultsToFiles( periodicOrbitModel->getLibrationPointNumber( ), periodicOrbitModel->getOrbitType( ), initialConditions, differentialCorrections );
}

void createPeriodicOrbitInitialConditions(
        const boost::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const boost::shared_ptr< tudat::cr3bp::CR3BPPeriodicOrbitModel > periodicOrbitModel,
        const boost::function< double( const Eigen::Vector6d& ) > pseudoArcLengthFunction )

{

    std::vector< Eigen::VectorXd > initialConditions;
    std::vector< Eigen::VectorXd > differentialCorrections;

    // Perform first two iteration
    std::pair< Eigen::Vector6d, double > richardsonThirdOrderApproximationResultIteration1 =
            periodicOrbitModel->getFirstInitialStateGuess( );
    std::pair< Eigen::Vector6d, double > richardsonThirdOrderApproximationResultIteration2 =
            periodicOrbitModel->getSecondInitialStateGuess( );

    correctPeriodicOrbitInitialState(
                richardsonThirdOrderApproximationResultIteration1.first, richardsonThirdOrderApproximationResultIteration1.second, 0,
                 periodicOrbitModel, integratorSettings, initialConditions, differentialCorrections );

    correctPeriodicOrbitInitialState(
                richardsonThirdOrderApproximationResultIteration2.first, richardsonThirdOrderApproximationResultIteration2.second, 1,
                periodicOrbitModel, integratorSettings, initialConditions, differentialCorrections );


    createPeriodicOrbitInitialConditionsFromExistingData(
                integratorSettings,
                periodicOrbitModel,
                initialConditions,
                differentialCorrections,
                pseudoArcLengthFunction );

}


