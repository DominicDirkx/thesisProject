#ifndef TUDATBUNDLE_CREATEINITIALCONDITIONS_H
#define TUDATBUNDLE_CREATEINITIALCONDITIONS_H



#include <string>
#include <vector>
#include <map>
#include <iostream>

#include <Eigen/Core>

#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/Astrodynamics/Propagators/stateDerivativeCircularRestrictedThreeBodyProblem.h"

namespace tudat
{

namespace cr3bp
{

enum CR3BPPeriodicOrbitTypes
{
    horizontal_lyapunov_orbit,
    vertical_lyapunov_orbit,
    halo_orbit,
    axial_orbit
};

class CR3BPPeriodicOrbitModel
{
public:
    CR3BPPeriodicOrbitModel(
            const CR3BPPeriodicOrbitTypes orbitType, const int librationPointNumber,
            const double massParameter );

    boost::function< Eigen::Matrix< double, 6, 7 >( const double, const Eigen::Matrix< double, 6, 7 >& ) >
    getStateDerivativeFunctionWithStateTransition( )
    {
        return boost::bind( &propagators::StateDerivativeCircularRestrictedThreeBodyProblem::computeStateDerivativeWithStateTransitionMatrix,
                            boost::make_shared< propagators::StateDerivativeCircularRestrictedThreeBodyProblem  >( massParameter_ ), _1, _2 );
    }

    boost::function< Eigen::Matrix< double, 6, 1 >( const double, const Eigen::Matrix< double, 6, 1 >& ) >
    getStateDerivativeFunction( )
    {
        return boost::bind( &propagators::StateDerivativeCircularRestrictedThreeBodyProblem::computeStateDerivative,
                            boost::make_shared< propagators::StateDerivativeCircularRestrictedThreeBodyProblem  >( massParameter_ ), _1, _2 );
    }

    bool isDifferentialCorrectionConverged( )
    {

    }

    bool terminateNumericalContinuation( )
    {

    }

    std::pair< double, double > getStateDeviation( )
    {

    }

private:

    CR3BPPeriodicOrbitTypes orbitType_;

    int librationPointNumber_;

    double massParameter_;

};

}

}

#endif  // TUDATBUNDLE_CREATEINITIALCONDITIONS_H
