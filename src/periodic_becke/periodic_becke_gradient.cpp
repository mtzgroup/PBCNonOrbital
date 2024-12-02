
#include "periodic_becke_kernel.h"
#include "../periodic_grid.h"
#include "../periodic_types.h" // For LatticeVector

#include "../helper.h" // For CUERR
#include "../gpubox.h"

static inline int min(const int a, const int b) { return (a < b) ? a : b; }

namespace PeriodicBox
{
    void beckeGradientLauncher(GPUBox* gpu, const int i_thread, const int n_point_total, const GridPoint* points,
                               const int n_atom, const double* atom_xyz, const double* atom_radius, const LatticeVector unit_cell,
                               const double* epsilon_xc, double* gradient)
    {
        const int n_thread = GPUBox::NumGPUstreams();
        const int n_point_job = n_point_total / n_thread + ((i_thread == n_thread - 1) ? (n_point_total % n_thread) : 0 );
        const int i_point_begin = (n_point_total / n_thread) * i_thread;

        double* d_interatomic_quantities;
        float4* d_atoms;
        set_atom_and_interatomic_quantities(gpu, d_interatomic_quantities, d_atoms, n_atom, atom_xyz, atom_radius);
        CUERR;

        double* h_point_xyzw = (double*)gpu->cpuAlloc(sizeof(double) * n_point_job * 4);
        double* h_point_x = h_point_xyzw + n_point_job * 0;
        double* h_point_y = h_point_xyzw + n_point_job * 1;
        double* h_point_z = h_point_xyzw + n_point_job * 2;
        double* h_point_w = h_point_xyzw + n_point_job * 3;

        int* h_i_atom_for_point = (int*)gpu->cpuAlloc(sizeof(int) * n_point_job);
        int* n_point_count_each_atom  = (int*)gpu->cpuAlloc(sizeof(int) * n_atom);
        int* i_point_offset_each_atom = (int*)gpu->cpuAlloc(sizeof(int) * n_atom);
        memset(n_point_count_each_atom,  0, sizeof(int) * n_atom);
        memset(i_point_offset_each_atom, 0, sizeof(int) * n_atom);

        for (int i_point = i_point_begin; i_point < i_point_begin + n_point_job; i_point++) {
            const int i_atom = points[i_point].i_atom;
            n_point_count_each_atom[i_atom]++;
        }

        int i_point_offset = 0;
        for (int i_atom = 0; i_atom < n_atom; i_atom++) {
            i_point_offset_each_atom[i_atom] = i_point_offset;
            i_point_offset += n_point_count_each_atom[i_atom];
        }

        for (int i_point = i_point_begin; i_point < i_point_begin + n_point_job; i_point++) {
            const int i_atom = points[i_point].i_atom;
            const int i_reordered = i_point_offset_each_atom[i_atom]++;
            h_point_x[i_reordered] = points[i_point].x;
            h_point_y[i_reordered] = points[i_point].y;
            h_point_z[i_reordered] = points[i_point].z;
            if (points[i_point].w_total != 0.0)
                h_point_w[i_reordered] = epsilon_xc[i_point] * points[i_point].w_fixed;
            else
                h_point_w[i_reordered] = 0.0;
            h_i_atom_for_point[i_reordered] = i_atom;
        }
        gpu->cpuFree(i_point_offset_each_atom);
        gpu->cpuFree(n_point_count_each_atom);

        double* d_point_xyzw = (double*)gpu->gpuAlloc(sizeof(double) * n_point_job * 4);
        double* d_point_x = d_point_xyzw + n_point_job * 0;
        double* d_point_y = d_point_xyzw + n_point_job * 1;
        double* d_point_z = d_point_xyzw + n_point_job * 2;
        double* d_point_w = d_point_xyzw + n_point_job * 3;
        int* d_i_atom_for_point = (int*)gpu->gpuAlloc(sizeof(int) * n_point_job);
        cudaMemcpy(d_point_xyzw, h_point_xyzw, sizeof(double) * n_point_job * 4, cudaMemcpyHostToDevice); CUERR;
        cudaMemcpy(d_i_atom_for_point, h_i_atom_for_point, sizeof(int) * n_point_job, cudaMemcpyHostToDevice); CUERR;

        // Split n_point_job * n_atom jobs into small pieces
        size_t remaining_gpu_memory_bytes = gpu->MaxAlloc() - 10 * gpu->Align(); // Henry 20240825: I don't under where this 10 comes from, it seems like a very conservative memory strategy 
        remaining_gpu_memory_bytes -= sizeof(float4) * n_atom;  // d_atoms 
        remaining_gpu_memory_bytes -= sizeof(double) * n_atom * n_atom;  // d_interatomic_quantities
        remaining_gpu_memory_bytes -= n_point_job * (4 * sizeof(double) + sizeof(int)); // d_point_xyzw and d_i_atom_for_point

        const int n_bytes_per_point = 4 * n_atom * sizeof(double); // d_gradient_cache_xyz
        const int max_points_for_job = remaining_gpu_memory_bytes / n_bytes_per_point;
        int n_point_per_grid;
        if (n_point_job < max_points_for_job)
            n_point_per_grid = n_point_job;
        else
            n_point_per_grid = (max_points_for_job / weight_gradient_block_x_dimension) * weight_gradient_block_x_dimension;
        n_point_per_grid = min(n_point_per_grid, weight_gradient_max_grid_dimension * weight_gradient_block_x_dimension);

        double* d_gradient_cache_xyz = (double*)gpu->gpuAlloc(sizeof(double) * n_point_per_grid * n_atom * 3);
        double* d_gradient_cache_x = d_gradient_cache_xyz + n_point_per_grid * n_atom * 0;
        double* d_gradient_cache_y = d_gradient_cache_xyz + n_point_per_grid * n_atom * 1;
        double* d_gradient_cache_z = d_gradient_cache_xyz + n_point_per_grid * n_atom * 2;
        cudaMemset(d_gradient_cache_x, 0, sizeof(double) * n_point_per_grid * n_atom); CUERR;
        cudaMemset(d_gradient_cache_y, 0, sizeof(double) * n_point_per_grid * n_atom); CUERR;
        cudaMemset(d_gradient_cache_z, 0, sizeof(double) * n_point_per_grid * n_atom); CUERR;

        dim3 grid_dim;
        grid_dim.y = n_atom / weight_gradient_block_y_dimension;
        if (n_atom % weight_gradient_block_y_dimension != 0)
            grid_dim.y++;
        const dim3 block_dim(weight_gradient_block_x_dimension, weight_gradient_block_y_dimension);

        for (int i_grid = 0; i_grid < n_point_job; i_grid += n_point_per_grid) {
            const int n_point_this_grid = min(n_point_per_grid, n_point_job - i_grid);
            grid_dim.x = n_point_this_grid / weight_gradient_block_x_dimension;
            if (n_point_this_grid % weight_gradient_block_x_dimension != 0)
                grid_dim.x++;

            weight_gradient_compute(grid_dim, block_dim, n_point_this_grid, n_atom, n_point_per_grid,
                                    d_point_x + i_grid, d_point_y + i_grid, d_point_z + i_grid, d_point_w + i_grid, d_atoms,
                                    d_i_atom_for_point + i_grid, d_gradient_cache_x, d_gradient_cache_y, d_gradient_cache_z,
                                    d_interatomic_quantities, default_switch_function_threshold, default_image_cutoff_radius, unit_cell); CUERR;

            // Center atom gradient
            int center_atom_kernel_grid_dim = n_point_this_grid / weight_gradient_center_atom_block_dimension;
            if (n_point_this_grid % weight_gradient_center_atom_block_dimension != 0)
                center_atom_kernel_grid_dim++;
            weight_gradient_center_atom(center_atom_kernel_grid_dim, weight_gradient_center_atom_block_dimension, n_point_this_grid, n_atom,
                                        n_point_per_grid, d_i_atom_for_point + i_grid, d_gradient_cache_x, d_gradient_cache_y, d_gradient_cache_z); CUERR;
        }

        double* h_gradient = (double*)gpu->cpuAlloc(sizeof(double) * n_atom * 3);
        double* d_gradient = (double*)gpu->gpuAlloc(sizeof(double) * n_atom * 3);
        weight_gradient_sum_over_point(n_atom * 3, weight_gradient_sum_over_point_dimension, n_point_per_grid, d_gradient,
                                       d_gradient_cache_x, d_gradient_cache_y, d_gradient_cache_z); CUERR;

        cudaMemcpy(h_gradient, d_gradient, sizeof(double) * n_atom * 3, cudaMemcpyDeviceToHost); CUERR;
        for (int i = 0; i < n_atom * 3; i++)
            gradient[i] += h_gradient[i];

        gpu->gpuFree(d_gradient);
        gpu->cpuFree(h_gradient);
        gpu->gpuFree(d_gradient_cache_xyz);
        gpu->gpuFree(d_i_atom_for_point);
        gpu->gpuFree(d_point_xyzw);
        gpu->cpuFree(h_i_atom_for_point);
        gpu->cpuFree(h_point_xyzw);
        gpu->gpuFree(d_interatomic_quantities);
        gpu->gpuFree(d_atoms);
    }
}
