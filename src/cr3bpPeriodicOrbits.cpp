#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "cr3bpPeriodicOrbits.h"
#include "richardsonThirdOrderApproximation.h"

namespace tudat
{

namespace cr3bp
{

double getEarthMoonPeriodicOrbitAmplitude( const int librationPointNumber, const CR3BPPeriodicOrbitTypes orbitType, const int guessIteration )
{
    if( librationPointNumber < 1 ||librationPointNumber > 2 )
    {
        throw std::runtime_error( "Error when getting Earth-Moon periodic orbit amplitude, libration point " +
                                  std::to_string( librationPointNumber ) + " is not supported." );
    }
    double amplitude = TUDAT_NAN;
    if( guessIteration == 0 )
    {
        if (orbitType == horizontal_lyapunov_orbit )
        {
            if (librationPointNumber == 1 )
            {
                amplitude = 1.0e-3;
            }
            else if (librationPointNumber == 2 )
            {
                amplitude = 1.0e-4;
            }
        }
        else if (orbitType == vertical_lyapunov_orbit )
        {
            if (librationPointNumber == 1)
            {
                amplitude = 1.0e-1;
            }
            else if (librationPointNumber == 2)
            {
                amplitude = 1.0e-1;
            }
        }
        else if (orbitType == halo_orbit )
        {
            if (librationPointNumber == 1)
            {
                amplitude = -1.1e-1;
            }
            else if (librationPointNumber == 2)
            {
                amplitude = 1.5e-1;
            }
        }
    }
    else if( guessIteration == 1 )
    {

        if (orbitType == horizontal_lyapunov_orbit )
        {
            if (librationPointNumber == 1)
            {
                amplitude = 1.0e-4;
            }
            else if (librationPointNumber == 2)
            {
                amplitude = 1.0e-3;
            }
        }
        else if (orbitType == vertical_lyapunov_orbit )
        {
            if (librationPointNumber == 1)
            {
                amplitude = 2.0e-1;
            }
            else if (librationPointNumber == 2)
            {
                amplitude = 2.0e-1;
            }
        }
        else if (orbitType == halo_orbit )
        {
            if (librationPointNumber == 1)
            {
                amplitude = -1.2e-1;
            }
            else if (librationPointNumber == 2)
            {
                amplitude = 1.6e-1;
            }
        }
    }
    else
    {
        throw std::runtime_error( "Error when getting Earth-Moon periodic orbit amplitude, iteration " +
                                  std::to_string( guessIteration ) + " is not supported." );
    }

    return amplitude;
}

std::pair< Eigen::Vector6d, double >  getLibrationPointPeriodicOrbitInitialStateVectorGuess(
        const int librationPointNumber, const CR3BPPeriodicOrbitTypes orbitType, const int guessIteration,
        const boost::function< double( const int librationPointNumber, const CR3BPPeriodicOrbitTypes orbitType, const int guessIteration ) > getAmplitude )
{
    double amplitude = getAmplitude( librationPointNumber, orbitType, guessIteration );
    std::pair< Eigen::Vector6d, double >  richardsonThirdOrderApproximationResult =
            richardsonApproximationLibrationPointPeriodicOrbit( orbitType, librationPointNumber, amplitude);

    return richardsonThirdOrderApproximationResult;
}

CR3BPPeriodicOrbitModel::CR3BPPeriodicOrbitModel(
        const CR3BPPeriodicOrbitTypes orbitType,
        const int librationPointNumber,
        const double massParameter,
        const double orbitAmplitudeForFirstInitialStateGuess,
        const double orbitAmplitudeForSecondInitialStateGuess,
        const int maximumNumberOfDifferentialCorrectionIterations,
        const double maximumAllowedPositionDeviationFromPeriodicOrbit,
        const double maximumAllowedVelocityDeviationFromPeriodicOrbit,
        const double maximumAllowedEigenvalueDeviation ):
    orbitType_( orbitType ), librationPointNumber_( librationPointNumber ), massParameter_( massParameter ),
    maximumNumberOfDifferentialCorrectionIterations_( maximumNumberOfDifferentialCorrectionIterations ),
    maximumAllowedPositionDeviationFromPeriodicOrbit_( maximumAllowedPositionDeviationFromPeriodicOrbit ),
    maximumAllowedVelocityDeviationFromPeriodicOrbit_( maximumAllowedVelocityDeviationFromPeriodicOrbit ),
    maximumAllowedEigenvalueDeviation_( maximumAllowedEigenvalueDeviation )
{
    firstInitialStateGuess_ =
            richardsonApproximationLibrationPointPeriodicOrbit( orbitType_, librationPointNumber_, orbitAmplitudeForFirstInitialStateGuess );
    firstInitialStateGuess_ =
            richardsonApproximationLibrationPointPeriodicOrbit( orbitType_, librationPointNumber_, orbitAmplitudeForSecondInitialStateGuess );
    stateDerivativeModel_ = boost::make_shared< propagators::StateDerivativeCircularRestrictedThreeBodyProblem  >( massParameter_ );
}

CR3BPPeriodicOrbitModel::CR3BPPeriodicOrbitModel(
        const CR3BPPeriodicOrbitTypes orbitType,
        const int librationPointNumber,
        const double massParameter,
        const std::pair< Eigen::Vector6d, double > firstInitialStateGuess,
        const std::pair< Eigen::Vector6d, double > secondInitialStateGuess,
        const int maximumNumberOfDifferentialCorrectionIterations,
        const double maximumAllowedPositionDeviationFromPeriodicOrbit,
        const double maximumAllowedVelocityDeviationFromPeriodicOrbit,
        const double maximumAllowedEigenvalueDeviation ):
    orbitType_( orbitType ), librationPointNumber_( librationPointNumber ), massParameter_( massParameter ),
    firstInitialStateGuess_( firstInitialStateGuess ), secondInitialStateGuess_( secondInitialStateGuess ),
    maximumNumberOfDifferentialCorrectionIterations_( maximumNumberOfDifferentialCorrectionIterations ),
    maximumAllowedPositionDeviationFromPeriodicOrbit_( maximumAllowedPositionDeviationFromPeriodicOrbit ),
    maximumAllowedVelocityDeviationFromPeriodicOrbit_( maximumAllowedVelocityDeviationFromPeriodicOrbit ),
    maximumAllowedEigenvalueDeviation_( maximumAllowedEigenvalueDeviation )
{
    stateDerivativeModel_ = boost::make_shared< propagators::StateDerivativeCircularRestrictedThreeBodyProblem  >( massParameter_ );
}

Eigen::VectorXd CR3BPPeriodicOrbitModel::computeDifferentialCorrection(
        const Eigen::Matrix6d& stmPartOfStateVectorInMatrixForm,
        const Eigen::Vector6d& cartesianState,
        const double currentTime )
{
    Eigen::Vector7d differentialCorrection;
    Eigen::Vector3d corrections;
    Eigen::Matrix3d updateMatrix;
    Eigen::Vector3d multiplicationMatrix;

    // Compute the accelerations and velocities (in X- and Z-direction) on the spacecraft and put them in a 2x1 vector.
    Eigen::VectorXd cartesianAccelerations = stateDerivativeModel_->computeStateDerivative( currentTime, cartesianState );

    // If type is axial, the desired state vector has the form [x, 0, 0, 0, ydot, zdot] and requires a differential correction for {x, ydot, T/2}
    if (orbitType_== axial_orbit )
    {
        // Check which deviation is larger: x-velocity or z-position.
        if ( std::abs(cartesianState(2)) < std::abs(cartesianState(3)) and !xPositionFixed )
        {
            // Correction on {x, ydot, T/2} for constant {zdot}

            // Set the correct multiplication matrix (state at T/2)
            multiplicationMatrix << cartesianState(1), cartesianState(2), cartesianState(3);

            // Compute the update matrix.
            updateMatrix << stmPartOfStateVectorInMatrixForm(1, 0), stmPartOfStateVectorInMatrixForm(1, 4), cartesianState(4),
                    stmPartOfStateVectorInMatrixForm(2, 0), stmPartOfStateVectorInMatrixForm(2, 4), cartesianState(5),
                    stmPartOfStateVectorInMatrixForm(3, 0), stmPartOfStateVectorInMatrixForm(3, 4), cartesianAccelerations(3);

            // Compute the necessary differential correction.
            corrections = updateMatrix.inverse() * multiplicationMatrix;

            // Put corrections in correct format.
            differentialCorrection.setZero( );
            differentialCorrection(0) = -corrections(0);
            differentialCorrection(4) = -corrections(1);
            differentialCorrection(6) = -corrections(2);
        }
        else
        {
            // Correction on {ydot, zdot T/2} for constant {x}

            // Set the correct multiplication matrix (state at T/2)
            multiplicationMatrix << cartesianState(1), cartesianState(2), cartesianState(3);

            // Compute the update matrix.
            updateMatrix << stmPartOfStateVectorInMatrixForm(1, 4), stmPartOfStateVectorInMatrixForm(1, 5), cartesianState(4),
                    stmPartOfStateVectorInMatrixForm(2, 4), stmPartOfStateVectorInMatrixForm(2, 5), cartesianState(5),
                    stmPartOfStateVectorInMatrixForm(3, 4), stmPartOfStateVectorInMatrixForm(3, 5), cartesianAccelerations(3);

            // Compute the necessary differential correction.
            corrections = updateMatrix.inverse() * multiplicationMatrix;

            // Put corrections in correct format.
            differentialCorrection.setZero( );
            differentialCorrection(4) = -corrections(0);
            differentialCorrection(5) = -corrections(1);
            differentialCorrection(6) = -corrections(2);
        }
    }

    // If type is not axial, the desired state vector has the form [x, 0, z, 0, ydot, 0] and requires a differential correction for either {z, ydot, T/2} or {x, ydot, T/2}
    else
    {

        // Check which deviation is larger: x-velocity or z-velocity.
        if ( std::abs(cartesianState(3)) < std::abs(cartesianState(5)) || orbitType_== horizontal_lyapunov_orbit ||
             (orbitType_== halo_orbit and librationPointNumber_ == 2) )
        {
            // Correction on {z, ydot, T/2} for constant {x}

            // Set the correct multiplication matrix (state at T/2)
            multiplicationMatrix << cartesianState(1), cartesianState(3), cartesianState(5);

            // Compute the update matrix.
            updateMatrix         << stmPartOfStateVectorInMatrixForm(1, 2), stmPartOfStateVectorInMatrixForm(1, 4), cartesianState(4),
                    stmPartOfStateVectorInMatrixForm(3, 2), stmPartOfStateVectorInMatrixForm(3, 4), cartesianAccelerations(3),
                    stmPartOfStateVectorInMatrixForm(5, 2), stmPartOfStateVectorInMatrixForm(5, 4), cartesianAccelerations(5);

            // Compute the necessary differential correction.
            corrections = updateMatrix.inverse() * multiplicationMatrix;

            // Put corrections in correct format.
            differentialCorrection.setZero();
            differentialCorrection(2) = -corrections(0);
            differentialCorrection(4) = -corrections(1);
            differentialCorrection(6) = -corrections(2);
        }
        else
        {
            // Correction on {x, ydot, T/2} for constant {z}
            // Set the correct multiplication matrix (state at T/2)
            multiplicationMatrix << cartesianState(1), cartesianState(3), cartesianState(5);

            // Compute the update matrix.
            updateMatrix << stmPartOfStateVectorInMatrixForm(1, 0), stmPartOfStateVectorInMatrixForm(1, 4), cartesianState(4),
                    stmPartOfStateVectorInMatrixForm(3, 0), stmPartOfStateVectorInMatrixForm(3, 4), cartesianAccelerations(3),
                    stmPartOfStateVectorInMatrixForm(5, 0), stmPartOfStateVectorInMatrixForm(5, 4), cartesianAccelerations(5);

            // Compute the necessary differential correction.
            corrections = updateMatrix.inverse() * multiplicationMatrix;

            // Put corrections in correct format.
            differentialCorrection.setZero();
            differentialCorrection(0) = -corrections(0);
            differentialCorrection(4) = -corrections(1);
            differentialCorrection(6) = -corrections(2);
        }
    }

    // Return differential correction.
    return differentialCorrection;
}


bool CR3BPPeriodicOrbitModel::terminateNumericalContinuation(
        const Eigen::Matrix6d& stmPartOfStateVectorInMatrixForm,
        const Eigen::Vector6d& cartesianState,
        const double currentTime )
{
    // Check termination conditions
    bool continueNumericalContinuation = true;
    if( stmPartOfStateVectorInMatrixForm != stmPartOfStateVectorInMatrixForm )
    {
        continueNumericalContinuation = false;
        std::cout << "\n\nNUMERICAL CONTINUATION STOPPED DUE TO NaN STATE TRANSITION MATROX \n\n" << std::endl;
    }
    else if ( differentialCorrections.at( differentialCorrections.size( ) - 1 ).segment(0, 6) == Eigen::VectorXd::Zero(6) )
    {
        continueNumericalContinuation = false;
        std::cout << "\n\nNUMERICAL CONTINUATION STOPPED DUE TO EXCEEDING MAXIMUM NUMBER OF ITERATIONS\n\n" << std::endl;
    }
    else
    {        // Check eigenvalue condition (at least one pair equalling a real one)
        continueNumericalContinuation = checkMonodromyEigenvalues( stmPartOfStateVectorInMatrixForm );
    }
    return continueNumericalContinuation;
}

bool CR3BPPeriodicOrbitModel::continueDifferentialCorrection(
        const Eigen::Vector6d stateVector, const int numberOfIterations )
{
    double maximumPositionDeviationToUse, maximumVelocityDeviationToUse;

    // If the maximum number of iterations has been reached, return a zero vector to stop the numerical continuation
    if ( numberOfIterations > maximumNumberOfDifferentialCorrectionIterations_ )
    {
        // Relax the periodicity constraints after exceeding the maximum number of iterations instead of termination
        maximumPositionDeviationToUse = 10.0 * maximumAllowedPositionDeviationFromPeriodicOrbit_;
        maximumVelocityDeviationToUse = 10.0 * maximumAllowedVelocityDeviationFromPeriodicOrbit_;
    }
    // Relax the maximum deviation requirements to compute the horizontal Lyapunov family in L2
    else if ( numberOfIterations > 10 && orbitType_ == horizontal_lyapunov_orbit && librationPointNumber_ == 2)
    {
        maximumPositionDeviationToUse = 10.0 * maximumAllowedPositionDeviationFromPeriodicOrbit_;
        maximumVelocityDeviationToUse = 10.0 * maximumAllowedVelocityDeviationFromPeriodicOrbit_;
    }
    else
    {
        maximumPositionDeviationToUse = maximumAllowedPositionDeviationFromPeriodicOrbit_;
        maximumVelocityDeviationToUse = maximumAllowedVelocityDeviationFromPeriodicOrbit_;
    }

    double positionDeviationFromPeriodicOrbit, velocityDeviationFromPeriodicOrbit;
    if ( orbitType_ == axial_orbit )
    {
        // Initial condition for axial family should be [x, 0, 0, 0, ydot, zdot]
        positionDeviationFromPeriodicOrbit = sqrt(pow(stateVector(1), 2) + pow(stateVector(2), 2));
        velocityDeviationFromPeriodicOrbit = sqrt(pow(stateVector(3), 2));
    }
    else
    {
        // Initial condition for other families should be [x, 0, y, 0, ydot, 0]
        positionDeviationFromPeriodicOrbit = sqrt(pow(stateVector(1), 2));
        velocityDeviationFromPeriodicOrbit = sqrt(pow(stateVector(3), 2) + pow(stateVector(5), 2));
    }

    return std::make_pair( ( positionDeviationFromPeriodicOrbit > maximumPositionDeviationToUse ) ||
                           ( velocityDeviationFromPeriodicOrbit > maximumVelocityDeviationToUse ) );
}


bool CR3BPPeriodicOrbitModel::checkMonodromyEigenvalues( const Eigen::Matrix6d &monodromyMatrix )
{
    bool moduleOneInsteadOfRealOne;

    // Exception for the horizontal Lyapunov family in Earth-Moon L2: eigenvalue may be of module one instead of a real one to compute a more extensive family
    if ( ( librationPointNumber_ == 2 ) && ( orbitType_ == horizontal_lyapunov_orbit ) )
    {
        moduleOneInsteadOfRealOne = true;
    }
    else
    {
        moduleOneInsteadOfRealOne = false;
    }

    // Initialize variables
    bool eigenvalueRealOne = false;

    // Reshape the STM for one period to matrix form and compute the eigenvalues
    Eigen::EigenSolver<Eigen::Matrix6d> eig( monodromyMatrix );

    // Determine whether the monodromy matrix contains at least one eigenvalue of real one within the maxEigenvalueDeviation
    for (int i = 0; i <= 5; i++)
    {
        if (std::abs(eig.eigenvalues().imag()(i)) < maxEigenvalueDeviation)
        {
            if (std::abs(eig.eigenvalues().real()(i) - 1.0) < maxEigenvalueDeviation)
            {
                eigenvalueRealOne = true;
            }
        }
    }

    // Optional argument to generalize the test from real one to module one
    if (moduleOneInsteadOfRealOne == true)
    {
        for (int i = 0; i <= 5; i++)
        {
            if (std::abs(std::abs(eig.eigenvalues()(i)) - 1.0 ) < maxEigenvalueDeviation)
            {
                eigenvalueRealOne = true;
            }
        }
    }

    return eigenvalueRealOne;
}

}
}

}
