#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "Tudat/Astrodynamics/Gravitation/librationPoint.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/celestialBodyConstants.h"

#include "createInitialConditions.h"
#include "computeManifolds.h"
#include "propagateOrbit.h"
//#include "completeInitialConditionsHaloFamily.h"
//#include "createInitialConditionsAxialFamily.h"
#include "connectManifoldsAtTheta.h"
//#include "omp.h"

//
//  TODO, Code refactoring:
//
//  - In correctPeriodicOrbitInitialState function, the appendDifferentialCorrectionResultsVector and appendResultsVector functions
//  should be changed, so that the output is written to a well organized struct, with clear identifiers for all information.
//  File output function should be created for consistency with Koen's format.
//  -  librationPointNr, orbitType, massParameter -> should be passed as a function that creates the state derivative model.
//  -  maximumNumberOfInitialConditions: Should be an input parameter to createPeriodicOrbitInitialConditions, as should settings for
//  continuation
//  DONE: - integrator settings should be an input variable to e.g. createPeriodicOrbitInitialConditions.

double massParameter = tudat::gravitation::circular_restricted_three_body_problem::computeMassParameter( tudat::celestial_body_constants::EARTH_GRAVITATIONAL_PARAMETER, tudat::celestial_body_constants::MOON_GRAVITATIONAL_PARAMETER );

int main (){

    // ================================
    // == Compute initial conditions ==
    // ================================

    double minimumStepSize   = std::numeric_limits<double>::epsilon( ); // 2.22044604925031e-16
    double maximumStepSize   = 100.0;//std::numeric_limits<double>::infinity( ); // 2.22044604925031e-16

    const double relativeErrorTolerance = 1.0E-10;
    const double absoluteErrorTolerance = 1.0E-14;

    boost::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > > integratorSettings =
            boost::make_shared< tudat::numerical_integrators::RungeKuttaVariableStepSizeSettings< > >
            ( tudat::numerical_integrators::rungeKuttaVariableStepSize, 0.0, 1.0E-5,
              tudat::numerical_integrators::RungeKuttaCoefficients::rungeKutta87DormandPrince, minimumStepSize, maximumStepSize,
              relativeErrorTolerance, absoluteErrorTolerance );
    integratorSettings->saveFrequency_ = 1;

    #pragma omp parallel num_threads(6)
    {
        #pragma omp for
        for (unsigned int i=1; i<=6; i++) {
            if (i ==1)
            {
                createPeriodicOrbitInitialConditions( 1, "horizontal", integratorSettings );
            }
            if (i ==2)
            {
                createPeriodicOrbitInitialConditions( 2, "horizontal", integratorSettings );
            }
            if (i ==3)
            {
                createPeriodicOrbitInitialConditions(1, "halo", integratorSettings );
            }
            if (i ==4)
            {
                createPeriodicOrbitInitialConditions(2, "halo", integratorSettings );
            }
            if (i ==5)
            {
                createPeriodicOrbitInitialConditions(1, "vertical", integratorSettings );
            }
            if (i ==6)
            {
                createPeriodicOrbitInitialConditions(2, "vertical", integratorSettings );
            }
        }
    }

//    // ======================================
//    // == Compute manifolds at theta (III) ==
//    // ======================================

//    double desiredJacobiEnergy = 3.1;
//    int thetaStoppingAngleMin = -180;
//    int thetaStoppingAngleMax = 0;
//    int numberOfTrajectoriesPerManifold = 5000;

//    for (int orbitTypeNumber = 0; orbitTypeNumber <= 0; orbitTypeNumber++) {

//        std::string orbitType;
//        if (orbitTypeNumber == 0) {
//            orbitType = "vertical";
//        } if (orbitTypeNumber == 1) {
//            orbitType = "horizontal";
//        } if (orbitTypeNumber == 2) {
//            orbitType = "halo";
//        }

//        std::vector<std::pair<double, Eigen::MatrixXd> > assembledResults;

//        #pragma omp parallel num_threads(20)
//        {
//            #pragma omp for
//            for (int i = thetaStoppingAngleMin; i <= thetaStoppingAngleMax; i++) {
//                double thetaStoppingAngle = static_cast<double>(i);
//                Eigen::MatrixXd minimumImpulseStateVectorsAtPoincare = connectManifoldsAtTheta(orbitType,
//                                                                                               thetaStoppingAngle,
//                                                                                               numberOfTrajectoriesPerManifold,
//                                                                                               desiredJacobiEnergy);
//                assembledResults.push_back(std::make_pair(thetaStoppingAngle, minimumImpulseStateVectorsAtPoincare));
//            }
//        }

//        std::ostringstream desiredJacobiEnergyStr;
//        std::string fileNameString;
//        desiredJacobiEnergyStr << std::setprecision(4) << desiredJacobiEnergy;
//        fileNameString = ("../data/raw/poincare_sections/" + orbitType + "_" + desiredJacobiEnergyStr.str() + "_minimum_impulse_connections.txt");
//        remove(fileNameString.c_str());
//        std::ofstream textFileAssembledResults(fileNameString.c_str());

//        textFileAssembledResults.precision(std::numeric_limits<double>::digits10);

//        for (int i = thetaStoppingAngleMin; i <= thetaStoppingAngleMax; i++) {
//            double thetaStoppingAngle = static_cast<double>(i);

//            for (unsigned int idx = 0; idx < assembledResults.size(); idx++) {
//                if (assembledResults.at(idx).first == thetaStoppingAngle) {
//                    textFileAssembledResults << std::left << std::scientific << std::setw(30) << thetaStoppingAngle
//                                             << std::setw(30)
//                                             << assembledResults.at(idx).second(0, 0) << std::setw(30)
//                                             << assembledResults.at(idx).second(0, 1) << std::setw(30)
//                                             << assembledResults.at(idx).second(0, 2) << std::setw(30)
//                                             << assembledResults.at(idx).second(0, 3) << std::setw(30)
//                                             << assembledResults.at(idx).second(0, 4) << std::setw(30)
//                                             << assembledResults.at(idx).second(0, 5) << std::setw(30)
//                                             << assembledResults.at(idx).second(0, 6) << std::setw(30)
//                                             << assembledResults.at(idx).second(0, 7) << std::setw(30)
//                                             << assembledResults.at(idx).second(1, 0) << std::setw(30)
//                                             << assembledResults.at(idx).second(1, 1) << std::setw(30)
//                                             << assembledResults.at(idx).second(1, 2) << std::setw(30)
//                                             << assembledResults.at(idx).second(1, 3) << std::setw(30)
//                                             << assembledResults.at(idx).second(1, 4) << std::setw(30)
//                                             << assembledResults.at(idx).second(1, 5) << std::setw(30)
//                                             << assembledResults.at(idx).second(1, 6) << std::setw(30)
//                                             << assembledResults.at(idx).second(1, 7) << std::endl;
//                }
//            }
//        }
//        textFileAssembledResults.close();
//    }

//    // ================================
//    // == Compute manifolds ==
//    // ================================
//    #pragma omp parallel num_threads(18)
//    {
//        #pragma omp for
//        for (unsigned int i=0; i<18; i++) {

//            std::string orbitType;
//            int librationPointNr;
//            int orbitIdOne;
//            double desiredJacobiEnergy;

//            if (i == 0){
//                orbitType = "horizontal";
//                librationPointNr = 1;
//                orbitIdOne = 808;
//                desiredJacobiEnergy = 3.05;
//            } if (i == 1){
//                orbitType = "horizontal";
//                librationPointNr = 1;
//                orbitIdOne = 577;
//                desiredJacobiEnergy = 3.1;
//            } if (i == 2){
//                orbitType = "horizontal";
//                librationPointNr = 1;
//                orbitIdOne = 330;
//                desiredJacobiEnergy = 3.15;
//            } if (i == 3){
//                orbitType = "horizontal";
//                librationPointNr = 2;
//                orbitIdOne = 1066;
//                desiredJacobiEnergy = 3.05;
//            } if (i == 4){
//                orbitType = "horizontal";
//                librationPointNr = 2;
//                orbitIdOne = 760;
//                desiredJacobiEnergy = 3.1;
//            } if (i == 5){
//                orbitType = "horizontal";
//                librationPointNr = 2;
//                orbitIdOne = 373;
//                desiredJacobiEnergy = 3.15;
//            } if (i == 6){
//                orbitType = "halo";
//                librationPointNr = 1;
//                orbitIdOne = 1235;
//                desiredJacobiEnergy = 3.05;
//            } if (i == 7){
//                orbitType = "halo";
//                librationPointNr = 1;
//                orbitIdOne = 836;
//                desiredJacobiEnergy = 3.1;
//            } if (i == 8){
//                orbitType = "halo";
//                librationPointNr = 1;
//                orbitIdOne = 358;
//                desiredJacobiEnergy = 3.15;
//            } if (i == 9){
//                orbitType = "halo";
//                librationPointNr = 2;
//                orbitIdOne = 1093;
//                desiredJacobiEnergy = 3.05;
//            } if (i == 10){
//                orbitType = "halo";
//                librationPointNr = 2;
//                orbitIdOne = 651;
//                desiredJacobiEnergy = 3.1;
//            } if (i == 11){
//                orbitType = "halo";
//                librationPointNr = 2;
//                orbitIdOne = 0;
//                desiredJacobiEnergy = 3.15;
//            } if (i == 12){
//                orbitType = "vertical";
//                librationPointNr = 1;
//                orbitIdOne = 1664;
//                desiredJacobiEnergy = 3.05;
//            } if (i == 13){
//                orbitType = "vertical";
//                librationPointNr = 1;
//                orbitIdOne = 1159;
//                desiredJacobiEnergy = 3.1;
//            } if (i == 14){
//                orbitType = "vertical";
//                librationPointNr = 1;
//                orbitIdOne = 600;
//                desiredJacobiEnergy = 3.15;
//            } if (i == 15){
//                orbitType = "vertical";
//                librationPointNr = 2;
//                orbitIdOne = 1878;
//                desiredJacobiEnergy = 3.05;
//            } if (i == 16){
//                orbitType = "vertical";
//                librationPointNr = 2;
//                orbitIdOne = 1275;
//                desiredJacobiEnergy = 3.1;
//            } if (i == 17){
//                orbitType = "vertical";
//                librationPointNr = 2;
//                orbitIdOne = 513;
//                desiredJacobiEnergy = 3.15;
//            }

//            Eigen::VectorXd selectedInitialConditions = readInitialConditionsFromFile(librationPointNr, orbitType, orbitIdOne, orbitIdOne + 1, massParameter);
//            Eigen::VectorXd refinedJacobiEnergyResult = refineOrbitJacobiEnergy(librationPointNr, orbitType, desiredJacobiEnergy,
//                                                                                selectedInitialConditions.segment(1, 6),
//                                                                                selectedInitialConditions(0),
//                                                                                selectedInitialConditions.segment(8, 6),
//                                                                                selectedInitialConditions(7), massParameter);
//            Eigen::VectorXd initialStateVector = refinedJacobiEnergyResult.segment(0, 6);
//            double orbitalPeriod               = refinedJacobiEnergyResult(6);

//            Eigen::MatrixXd fullInitialState = getFullInitialState( initialStateVector );
//            std::map< double, Eigen::Vector6d > stateHistory;
//            std::pair< Eigen::MatrixXd, double > endState = propagateOrbitToFinalCondition(
//                        fullInitialState, massParameter, integratorSettings, orbitalPeriod, 1, stateHistory, 100 );

//            writeStateHistoryToFile( stateHistory, orbitIdOne, orbitType, librationPointNr, 1000, false );

//            // ===============================================================
//            // == Compute manifolds based on precomputed initial conditions ==
//            // ===============================================================
//            computeManifolds(initialStateVector, orbitalPeriod, orbitIdOne, librationPointNr, orbitType);

////                    for (unsigned int librationPointNr = 2; librationPointNr <= 2; librationPointNr++) {

//            // =========================================
//            // == Load precomputed initial conditions ==
//            // =========================================

////            std::ifstream textFileInitialConditions("../data/raw/orbits/L" + std::to_string(librationPointNr) + "_" + orbitType + "_initial_conditions.txt");
////            std::vector<std::vector<double>> initialConditions;
////
////            if (textFileInitialConditions) {
////                std::string line;
////
////                while (std::getline(textFileInitialConditions, line)) {
////                    initialConditions.push_back(std::vector<double>());
////
////                    // Break down the row into column values
////                    std::stringstream split(line);
////                    double value;
////
////                    while (split >> value)
////                        initialConditions.back().push_back(value);
////                }
////            }
////            double orbitalPeriod1;
////            double orbitalPeriod2;
////            Eigen::VectorXd initialStateVector1;
////            Eigen::VectorXd initialStateVector2;

//            // ==============================================================================================
//            // == Complete initial conditions for the halo family, until connection to horizontal Lyapunov ==
//            // ==============================================================================================

////            orbitalPeriod1         = initialConditions[2][1];
////            initialStateVector1    = Eigen::VectorXd::Zero(6);
////            initialStateVector1(0) = initialConditions[2][2];
////            initialStateVector1(1) = initialConditions[2][3];
////            initialStateVector1(2) = initialConditions[2][4];
////            initialStateVector1(3) = initialConditions[2][5];
////            initialStateVector1(4) = initialConditions[2][6];
////            initialStateVector1(5) = initialConditions[2][7];
////
////            orbitalPeriod2         = initialConditions[1][1];
////            initialStateVector2    = Eigen::VectorXd::Zero(6);
////            initialStateVector2(0) = initialConditions[1][2];
////            initialStateVector2(1) = initialConditions[1][3];
////            initialStateVector2(2) = initialConditions[1][4];
////            initialStateVector2(3) = initialConditions[1][5];
////            initialStateVector2(4) = initialConditions[1][6];
////            initialStateVector2(5) = initialConditions[1][7];
////
////            completeInitialConditionsHaloFamily( initialStateVector1, initialStateVector2, orbitalPeriod1, orbitalPeriod2, librationPointNr );

//            // ==================================================================================================================
//            // == Create initial conditions for the axial family, based on the bifurcation from the horizontal Lyapunov family ==
//            // ==================================================================================================================

//            // Create initial conditions axial family
////            int orbitIdForBifurcationToAxial;
////            double offsetForBifurcationToAxial1;
////            double offsetForBifurcationToAxial2;
////            if (librationPointNr == 1) {
////                orbitIdForBifurcationToAxial = 938;  // Indices for L1 bifurcations: [181, 938, 1233]
////                offsetForBifurcationToAxial1 = -0.01;
////                offsetForBifurcationToAxial2 = -0.02;
////            } else {
////                orbitIdForBifurcationToAxial = 1262;  // Indices for L2 bifurcations: [353, 1262, 1514]
////                offsetForBifurcationToAxial1 = -8.170178408168770e-04;
////                offsetForBifurcationToAxial2 = -1.817017673845933e-03;
////            }
////
////            orbitalPeriod1         = initialConditions[orbitIdForBifurcationToAxial][1];
////            initialStateVector1    = Eigen::VectorXd::Zero(6);
////            initialStateVector1(0) = initialConditions[orbitIdForBifurcationToAxial][2];
////            initialStateVector1(1) = initialConditions[orbitIdForBifurcationToAxial][3];
////            initialStateVector1(2) = initialConditions[orbitIdForBifurcationToAxial][4];
////            initialStateVector1(3) = initialConditions[orbitIdForBifurcationToAxial][5];
////            initialStateVector1(4) = initialConditions[orbitIdForBifurcationToAxial][6];
////            initialStateVector1(5) = initialConditions[orbitIdForBifurcationToAxial][7] + offsetForBifurcationToAxial1;
////            Eigen::VectorXd stateVectorInclSTM;
////            stateVectorInclSTM     = writePeriodicOrbitToFile( initialStateVector1, librationPointNr, "axial", 0, orbitalPeriod1, massParameter);
////
////            orbitalPeriod2         = initialConditions[orbitIdForBifurcationToAxial][1];
////            initialStateVector2    = Eigen::VectorXd::Zero(6);
////            initialStateVector2(0) = initialConditions[orbitIdForBifurcationToAxial][2];
////            initialStateVector2(1) = initialConditions[orbitIdForBifurcationToAxial][3];
////            initialStateVector2(2) = initialConditions[orbitIdForBifurcationToAxial][4];
////            initialStateVector2(3) = initialConditions[orbitIdForBifurcationToAxial][5];
////            initialStateVector2(4) = initialConditions[orbitIdForBifurcationToAxial][6];
////            initialStateVector2(5) = initialConditions[orbitIdForBifurcationToAxial][7] + offsetForBifurcationToAxial2;
////
////            stateVectorInclSTM     = writePeriodicOrbitToFile( initialStateVector2, librationPointNr, "axial", 1, orbitalPeriod2, massParameter);
////
////            createInitialConditionsAxialFamily(initialStateVector1, initialStateVector2, orbitalPeriod1, orbitalPeriod2, librationPointNr);



//        }
//    }

    return 0;
}
