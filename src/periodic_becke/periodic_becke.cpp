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

    // void beckeGradThrd(GPUBox* gpu, int tid, int nPt, const GridPoint* Points, int nAtoms, const double* xyz, const double* Radius, const double* func, double* Grad);
    // void BeckeGrad(int nPt, const GridPoint* Points, int nAtoms, const double* xyz, const double* Radius, const double* func, double* Gradient)
    // {
    //     VecPool GradPool(GPUBox::NumGPUstreams(), 3*nAtoms);
    //     GradPool.zero();

    //     PARALLEL_FOR (int tid = 0; tid < GPUBox::NumGPUstreams(); tid++){
    //         GPUBox* gpu = GPUBox::Checkout();
    //         double* t_grad = GradPool.checkout();
    //         beckeGradThrd(gpu, tid, nPt, Points, nAtoms, xyz, Radius, func, t_grad);
    //         GradPool.checkin(t_grad);
    //         GPUBox::Checkin(gpu);
    //     }

    //     GradPool.reduce(Gradient);
    // }
}
