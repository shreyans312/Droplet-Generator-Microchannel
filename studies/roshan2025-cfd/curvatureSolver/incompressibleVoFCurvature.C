// SPDX-License-Identifier: GPL-3.0-or-later
// CSF colour-function sharpening from Antritter et al. (2020), Eq. (11).
// The phase-advection and density fields stay unchanged. The surface force
// uses the original face-normal discretisation shared with the pressure solve.
#include "incompressibleVoF.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceInterpolate.H"
#include "fvcSnGrad.H"

namespace Foam
{
namespace solvers
{
class incompressibleVoFCurvature : public incompressibleVoF
{
    scalar sharpening_;

protected:
    virtual tmp<surfaceScalarField> surfaceTensionForce() const
    {
        if (sharpening_ == 0)
        {
            return incompressibleVoF::surfaceTensionForce();
        }
        const volScalarField sharpened
        (
            IOobject
            (
                "surfaceForceAlpha", alpha1.instance(), alpha1.mesh(),
                IOobject::NO_READ, IOobject::NO_WRITE, false
            ),
            min(max((alpha1 - sharpening_/2)/(1 - sharpening_), scalar(0)), scalar(1))
        );
        return fvc::interpolate(interface.sigmaK())*fvc::snGrad(sharpened);
    }

public:
    TypeName("incompressibleVoFCurvature");
    incompressibleVoFCurvature(fvMesh& mesh)
    : incompressibleVoF(mesh),
      sharpening_(mesh.time().controlDict().lookupOrDefault<scalar>("surfaceTensionSharpening", .5))
    {
        if (sharpening_ < 0 || sharpening_ >= 1)
        {
            FatalErrorInFunction << "surfaceTensionSharpening must be in [0,1)" << exit(FatalError);
        }
        Info<< "Curvature CSF force sharpening = " << sharpening_ << endl;
    }
};
defineTypeNameAndDebug(incompressibleVoFCurvature, 0);
addToRunTimeSelectionTable(solver, incompressibleVoFCurvature, fvMesh);
}
}
