#ifndef PERIODIC_BECKE_PERIODIC_BECKE_H_
#define PERIODIC_BECKE_PERIODIC_BECKE_H_

#include "../periodic_types.h" // For LatticeVector
namespace PeriodicBox { struct GridPoint; }

namespace PeriodicBox
{
    void BeckeWeights(
            const int n_point,
            GridPoint* points, // The .w_total field will get modified
            const int n_atom,
            const double* atom_xyz,
            const double* atom_radius,
            const LatticeVector unit_cell
            );

    // void BeckeGrad(
    //         int nPt,
    //         const GridPoint* points,
    //         int nAtoms,
    //         const double* xyz,
    //         const double* Radius,
    //         const double* func,
    //         double* Grad
    //         );
}

#endif
