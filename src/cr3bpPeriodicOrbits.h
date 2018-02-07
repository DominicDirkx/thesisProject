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
#include "Tudat/Astrodynamics/BasicAstrodynamics/celestialBodyConstants.h"
#include "Tudat/Astrodynamics/Propagators/stateDerivativeCircularRestrictedThreeBodyProblem.h"
#include "Tudat/Astrodynamics/Gravitation/librationPoint.h"

#include "cr3bpPeriodicOrbitTypes.h"

namespace tudat
{

namespace cr3bp
{

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
            const double maximumAllowedEigenvalueDeviation,
            const int maximumNumberOfNumericalContinuations );

    CR3BPPeriodicOrbitModel(
            const CR3BPPeriodicOrbitTypes orbitType,
            const int librationPointNumber,
            const double massParameter,
            const std::pair< Eigen::Vector6d, double > firstInitialStateGuess,
            const std::pair< Eigen::Vector6d, double > secondInitialStateGuess,
            const int maximumNumberOfDifferentialCorrectionIterations,
            const double maximumAllowedPositionDeviationFromPeriodicOrbit,
            const double maximumAllowedVelocityDeviationFromPeriodicOrbit,
            const double maximumAllowedEigenvalueDeviation,
            const int maximumNumberOfNumericalContinuations );

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
            const double currentTime,
            const int numberOfIterations );


    bool terminateNumericalContinuation(
            const Eigen::Matrix6d& stmPartOfStateVectorInMatrixForm,
            const Eigen::Vector6d& cartesianState,
            const double currentTime,
            const int numberOfOrbitsGenerated );

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

    CR3BPPeriodicOrbitTypes getOrbitType( )
    {
        return orbitType_;
    }

    int getLibrationPointNumber( )
    {
        return librationPointNumber_;
    }

private:

    bool checkMonodromyEigenvalues( const Eigen::Matrix6d &monodromyMatrix );


    boost::shared_ptr< propagators::StateDerivativeCircularRestrictedThreeBodyProblem  > stateDerivativeModel_;

    CR3BPPeriodicOrbitTypes orbitType_;

    int librationPointNumber_;

    double massParameter_;


    int maximumNumberOfDifferentialCorrectionIterations_;

    int maximumNumberOfNumericalContinuations_;

    double maximumAllowedPositionDeviationFromPeriodicOrbit_;

    double maximumAllowedVelocityDeviationFromPeriodicOrbit_;

    double maximumAllowedEigenvalueDeviation_;

    std::pair< Eigen::Vector6d, double > firstInitialStateGuess_;

    std::pair< Eigen::Vector6d, double > secondInitialStateGuess_;
};

double getEarthMoonPeriodicOrbitAmplitude( const int librationPointNumber, const CR3BPPeriodicOrbitTypes orbitType, const int guessIteration );

std::pair< Eigen::Vector6d, double > getLibrationPointPeriodicOrbitInitialStateVectorGuess(
        const int librationPointNumber, const CR3BPPeriodicOrbitTypes orbitType, const int guessIteration,
        const boost::function< double( const int librationPointNumber, const std::string& orbitType, const int guessIteration ) > getAmplitude =
        getEarthMoonPeriodicOrbitAmplitude );

boost::shared_ptr< CR3BPPeriodicOrbitModel > createEarthMoonPeriodicOrbitModel(
        const int librationPointNumber,
        const CR3BPPeriodicOrbitTypes orbitType,
        const double massParameter =
        gravitation::circular_restricted_three_body_problem::computeMassParameter(
            tudat::celestial_body_constants::EARTH_GRAVITATIONAL_PARAMETER, tudat::celestial_body_constants::MOON_GRAVITATIONAL_PARAMETER ),
        const int maximumNumberOfDifferentialCorrectionIterations = 25,
        const double maximumAllowedPositionDeviationFromPeriodicOrbit = 1.0E-12,
        const double maximumAllowedVelocityDeviationFromPeriodicOrbit = 1.0E-12,
        const double maximumAllowedEigenvalueDeviation = 1.0E-3,
        const int maximumNumberOfNumericalContinuations = 1E4 );

}

}

#endif  // TUDATBUNDLE_CR3BPPERIODICORBITS_H
