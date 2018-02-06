#ifndef TUDATBUNDLE_CR3BPPERIODICORBITS_H
#define TUDATBUNDLE_CR3BPPERIODICORBITS_H



#include <string>
#include <vector>
#include <map>
#include <iostream>

#include <boost/function.hpp>
#include <boost/bind.hpp>

#include <Eigen/Core>

#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/Astrodynamics/Propagators/stateDerivativeCircularRestrictedThreeBodyProblem.h"

#include "cr3bpPeriodicOrbitTypes.h"

namespace tudat
{

namespace cr3bp
{


double getEarthMoonPeriodicOrbitAmplitude( const int librationPointNumber, const CR3BPPeriodicOrbitTypes orbitType, const int guessIteration );

std::pair< Eigen::Vector6d, double > getLibrationPointPeriodicOrbitInitialStateVectorGuess(
        const int librationPointNumber, const CR3BPPeriodicOrbitTypes orbitType, const int guessIteration,
        const boost::function< double( const int librationPointNumber, const std::string& orbitType, const int guessIteration ) > getAmplitude =
        getEarthMoonPeriodicOrbitAmplitude );

class CR3BPPeriodicOrbitModel
{
public:
    CR3BPPeriodicOrbitModel(
            const CR3BPPeriodicOrbitTypes orbitType,
            const int librationPointNumber,
            const double massParameter,
            const double orbitAmplitudeForFirstInitialStateGuess,
            const double orbitAmplitudeForSecondInitialStateGuess,
            const int maximumNumberOfDifferentialCorrectionIterations,
            const double maximumAllowedPositionDeviationFromPeriodicOrbit,
            const double maximumAllowedVelocityDeviationFromPeriodicOrbit,
            const double maximumAllowedEigenvalueDeviation );

    CR3BPPeriodicOrbitModel(
            const CR3BPPeriodicOrbitTypes orbitType,
            const int librationPointNumber,
            const double massParameter,
            const std::pair< Eigen::Vector6d, double > firstInitialStateGuess,
            const std::pair< Eigen::Vector6d, double > secondInitialStateGuess,
            const int maximumNumberOfDifferentialCorrectionIterations,
            const double maximumAllowedPositionDeviationFromPeriodicOrbit,
            const double maximumAllowedVelocityDeviationFromPeriodicOrbit,
            const double maximumAllowedEigenvalueDeviation );

    boost::function< Eigen::Matrix< double, 6, 7 >( const double, const Eigen::Matrix< double, 6, 7 >& ) >
    getStateDerivativeFunctionWithStateTransition( )
    {
        return boost::bind(
                    &propagators::StateDerivativeCircularRestrictedThreeBodyProblem::computeStateDerivativeWithStateTransitionMatrix,
                    stateDerivativeModel_, _1, _2 );
    }

    boost::function< Eigen::Matrix< double, 6, 1 >( const double, const Eigen::Matrix< double, 6, 1 >& ) >
    getStateDerivativeFunction( )
    {
        return boost::bind(
                    &propagators::StateDerivativeCircularRestrictedThreeBodyProblem::computeStateDerivative,
                    stateDerivativeModel_, _1, _2 );
    }

    Eigen::VectorXd computeDifferentialCorrection(
            const Eigen::Matrix6d& stmPartOfStateVectorInMatrixForm,
            const Eigen::Vector6d& cartesianState,
            const double currentTime );


    bool terminateNumericalContinuation(
            const Eigen::Matrix6d& stmPartOfStateVectorInMatrixForm,
            const Eigen::Vector6d& cartesianState,
            const double currentTime );

    bool continueDifferentialCorrection(
            const Eigen::Vector6d stateVector, const int numberOfIterations );

    std::pair< Eigen::Vector6d, double > getFirstInitialStateGuess( )
    {
        return firstInitialStateGuess_;
    }

    std::pair< Eigen::Vector6d, double > getSecondInitialStateGuess( )
    {
        return secondInitialStateGuess_;
    }

    double getMassParameter( )
    {
        return massParameter_;
    }

private:

    bool checkMonodromyEigenvalues( const Eigen::Matrix6d &monodromyMatrix, const bool moduleOneInsteadOfRealOne );


    boost::shared_ptr< propagators::StateDerivativeCircularRestrictedThreeBodyProblem  > stateDerivativeModel_;

    CR3BPPeriodicOrbitTypes orbitType_;

    int librationPointNumber_;

    double massParameter_;


    int maximumNumberOfDifferentialCorrectionIterations_;

    double maximumAllowedPositionDeviationFromPeriodicOrbit_;

    double maximumAllowedVelocityDeviationFromPeriodicOrbit_;

    double maximumAllowedEigenvalueDeviation_;

    std::pair< Eigen::Vector6d, double > firstInitialStateGuess_;

    std::pair< Eigen::Vector6d, double > secondInitialStateGuess_;
};

}

}

#endif  // TUDATBUNDLE_CR3BPPERIODICORBITS_H
