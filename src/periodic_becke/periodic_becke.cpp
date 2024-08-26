#include "periodic_becke.h"

#include "../helper.h"
#include "../gpubox.h"

namespace PeriodicBox
{
    void beckeWeightLauncher(GPUBox* gpu, const int i_thread, const int n_point_total, GridPoint* points,
                             const int n_atom, const double* atom_xyz, const double* atom_radius, const LatticeVector unit_cell);
    void BeckeWeights(const int n_point, GridPoint* points, const int n_atom, const double* atom_xyz, const double* atom_radius, const LatticeVector unit_cell)
    {
        PARALLEL_FOR (int i_thread = 0; i_thread < GPUBox::NumGPUstreams(); i_thread++) {
            GPUBox* gpu = GPUBox::Checkout();
            beckeWeightLauncher(gpu, i_thread, n_point, points, n_atom, atom_xyz, atom_radius, unit_cell);
            GPUBox::Checkin(gpu);
        }
    }

    void beckeGradientLauncher(GPUBox* gpu, const int i_thread, const int n_point_total, const GridPoint* points,
                               const int n_atom, const double* atom_xyz, const double* atom_radius, const LatticeVector unit_cell,
                               const double* epsilon_xc, double* gradient);
    void BeckeWeightGradient(const int n_point, const GridPoint* points,
                            const int n_atom, const double* atom_xyz, const double* atom_radius,
                            const LatticeVector unit_cell, const double* epsilon_xc, double* gradient)
    {
        VecPool gradient_copies(GPUBox::NumGPUstreams(), 3 * n_atom);
        gradient_copies.zero();

        PARALLEL_FOR (int i_thread = 0; i_thread < GPUBox::NumGPUstreams(); i_thread++) {
            GPUBox* gpu = GPUBox::Checkout();
            double* gradient_copy = gradient_copies.checkout();
            beckeGradientLauncher(gpu, i_thread, n_point, points, n_atom, atom_xyz, atom_radius, unit_cell, epsilon_xc, gradient_copy);
            gradient_copies.checkin(gradient_copy);
            GPUBox::Checkin(gpu);
        }

        gradient_copies.reduce(gradient);
    }
}
