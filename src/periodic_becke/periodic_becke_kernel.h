#ifndef PERIODIC_BECKE_PERIODIC_GRID_KERNEL_H_
#define PERIODIC_BECKE_PERIODIC_GRID_KERNEL_H_

// #define DFT_BOX_THRE 0
// #define MAX_WT_GRID_DIM 65535

// #define WGRAD_PTCTR_BLOCK   256
// #define WGRAD_DOSUMS_BLOCK  256
// #define WGRAD_BLOCK_YDIM_SP   1
// #define WGRAD_BLOCK_XDIM_SP 128
// #define WGRAD_BLOCK_YDIM_DP   1
// #define WGRAD_BLOCK_XDIM_DP 128

#include "../periodic_types.h" // For LatticeVector

#include <vector_types.h> // For float4

class GPUBox;

namespace PeriodicBox
{
    const double default_switch_function_threshold = 1.0e-12;
    // Henry 20240329: The following cutoff radius cannot be derived from any inequality. We got this numerical value (9 Angstrom) from
    //                 Towler, M. D.; Zupan, A.; Causà, M., Density functional theory in periodic systems using local Gaussian basis sets. Computer Physics Communications 1996, 98 (1), 181-205.
    const double default_image_cutoff_radius = 17.0; // Bohr

    void set_atom_and_interatomic_quantities(GPUBox* gpu, double*& d_interatomic_quantities, float4*& d_atoms,
                                             const int n_atom, const double* atom_xyz, const double* atom_radius);
    void weights(const int n_grid, const int n_block,
                 const int n_point, const double* d_point_x, const double* d_point_y, const double* d_point_z,
                 const float4* d_atoms, const int n_atom, const int* d_point_atom_center_index,
                 double* d_weight, const double* d_interatomic_quantities,
                 const double switch_function_threshold, const double image_cutoff_radius,
                 const LatticeVector unit_cell);
    // void Wgrad_ptctr_caller(const int GRID, const int BLOCK, const int numPoints, const int numAtoms, const int grad_pitch, const int* ctrs,
    //                         double* xGrad, double* yGrad, double* zGrad);
    // void Wgrad_dosums_caller(const int GRID, const int BLOCK, const int numPoints, double* final_grad,
    //                         const double* xGrad, const double* yGrad, const double* zGrad);
    // void Wgrad_cache_caller(const dim3 GRID, const dim3 BLOCK, const int numPoints, const int numAtoms, const int grad_pitch,
    //                         const double* xpts, const double* ypts, const double* zpts, const float4* c_atom,
    //                         double* ra_cache, double* p_cache,
    //                         const double2* cache, const double w_thre);
    // void Wgrad_kernel_caller(const dim3 GRID, const dim3 BLOCK, const int numPoints, const int numAtoms, const int grad_pitch, 
    //                         const double* xpts, const double* ypts, const double* zpts, const double* wpts,
    //                         const float4* c_atom, const int* ctrs,
    //                         double* xGrad, double* yGrad, double* zGrad,
    //                         const double* ra_cache, const double* p_cache, const double2* cache);
}

#endif
