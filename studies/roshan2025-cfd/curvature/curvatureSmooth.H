// SPDX-License-Identifier: GPL-3.0-or-later
// Curvature-only colour-function smoothing; transported alpha is never changed.
// Eq. (10), Antritter et al. (2020), Progress in Computational Fluid Dynamics
// 20(2), 59-70: two cell/face averaging passes with C_kappa = 0.5.
#ifndef curvatureSmooth_H
#define curvatureSmooth_H

#include "gaussGrad.H"
#include "fvcAverage.H"
#include "surfaceInterpolate.H"

namespace Foam
{
namespace fv
{
template<class Type>
class curvatureSmooth : public gradScheme<Type>
{
    label passes_;
    scalar weight_;

public:
    TypeName("curvatureSmooth");

    curvatureSmooth(const fvMesh& mesh, Istream& is)
    : gradScheme<Type>(mesh), passes_(readLabel(is)), weight_(readScalar(is))
    {
        if (passes_ < 0 || passes_ > 8 || weight_ < 0 || weight_ > 1)
        {
            FatalIOErrorInFunction(is)
                << "Require 0 <= passes <= 8 and 0 <= weight <= 1"
                << exit(FatalIOError);
        }
    }

    virtual tmp<VolField<typename outerProduct<vector, Type>::type>> calcGrad
    (
        const VolField<Type>& field,
        const word& name
    ) const
    {
        VolField<Type> smoothed
        (
            IOobject
            (
                "curvatureSmoothed_" + field.name(),
                field.instance(),
                field.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            field
        );
        for (label pass = 0; pass < passes_; ++pass)
        {
            smoothed = (1 - weight_)*smoothed
                + weight_*fvc::average(fvc::interpolate(smoothed));
            smoothed.correctBoundaryConditions();
        }
        return gaussGrad<Type>(this->mesh()).calcGrad(smoothed, name);
    }
};
}
}
#endif
