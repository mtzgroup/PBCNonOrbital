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

    void BeckeWeightGradient(
        const int n_point,
        const GridPoint* points,
        const int n_atom,
        const double* atom_xyz,
        const double* atom_radius,
        const LatticeVector unit_cell,
        const double* epsilon_xc,
        double* gradient
    );
}

#endif
