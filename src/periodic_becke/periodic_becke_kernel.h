#ifndef PERIODIC_BECKE_PERIODIC_GRID_KERNEL_H_
#define PERIODIC_BECKE_PERIODIC_GRID_KERNEL_H_

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

    const int weight_gradient_block_x_dimension = 128;
    const int weight_gradient_block_y_dimension = 1;
    const int weight_gradient_max_grid_dimension = 65535;
    const int weight_gradient_center_atom_block_dimension = 256;
    const int weight_gradient_sum_over_point_dimension = 256;

    void weight_gradient_cache(const dim3 n_grid, const dim3 n_block, const int n_point, const int n_atom, const int n_point_per_grid,
                               const double* d_point_x, const double* d_point_y, const double* d_point_z, const float4* d_atoms,
                               double* d_p_cache,
                               const double* d_interatomic_quantities, const double switch_function_threshold,
                               const double image_cutoff_radius, const LatticeVector unit_cell);
    void weight_gradient_compute(const dim3 n_grid, const dim3 n_block, const int n_point, const int n_atom, const int n_point_per_grid, 
                                 const double* d_point_x, const double* d_point_y, const double* d_point_z, const double* wpts,
                                 const float4* d_atoms, const int* d_i_atom_for_point,
                                 double* d_gradient_cache_x, double* d_gradient_cache_y, double* d_gradient_cache_z,
                                 const double* d_p_cache, const double* d_interatomic_quantities,
                                 const double image_cutoff_radius, const LatticeVector unit_cell);
    void weight_gradient_center_atom(const int n_grid, const int n_block, const int n_point, const int n_atom, const int n_point_per_grid, const int* d_i_atom_for_point,
                                     double* d_gradient_cache_x, double* d_gradient_cache_y, double* d_gradient_cache_z);
    void weight_gradient_sum_over_point(const int n_grid, const int n_block, const int n_point, double* final_gradient,
                                        const double* d_gradient_cache_x, const double* d_gradient_cache_y, const double* d_gradient_cache_z);
}

#endif
