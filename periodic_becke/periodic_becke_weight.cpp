
#include "periodic_becke_kernel.h"
#include "../periodic_grid.h"
#include "../periodic_types.h" // For LatticeVector

#include <utils.h> // From intbox, for CUERR

#include <gpubox.h>

static inline int min(const int a, const int b) { return (a < b) ? a : b; }

namespace PeriodicBox
{
    const int n_threads_per_block_weight = 256;

    void beckeWeightLauncher(GPUBox* gpu, const int i_thread, const int n_point_total, GridPoint* points,
                             const int n_atom, const double* atom_xyz, const double* atom_radius, const LatticeVector unit_cell)
    {
        const int n_thread = GPUBox::NumGPUstreams();

        double* d_interatomic_quantities;
        float4* d_atoms;
        set_atom_and_interatomic_quantities(gpu, d_interatomic_quantities, d_atoms, n_atom, atom_xyz, atom_radius);
        CUERR;

        const int n_point_job = n_point_total / n_thread + ((i_thread < n_point_total % n_thread) ? 1 : 0);

        double* h_point_x = (double*)gpu->cpuAlloc(sizeof(double) * 4 * n_point_job);
        double* h_point_y = h_point_x + n_point_job;
        double* h_point_z = h_point_y + n_point_job;
        double* h_point_w = h_point_z + n_point_job;
        int* h_point_atom_center_index = (int*)gpu->cpuAlloc(sizeof(int) * n_point_job);

        for (int i_global = i_thread, i_local = 0; i_global < n_point_total; i_global += n_thread, i_local++) {
            h_point_x[i_local] = points[i_global].x;
            h_point_y[i_local] = points[i_global].y;
            h_point_z[i_local] = points[i_global].z;
            h_point_w[i_local] = points[i_global].w_fixed;
            h_point_atom_center_index[i_local]  = points[i_global].i_atom;
        }

        double* d_point_x = (double*)gpu->gpuAlloc(sizeof(double) * 4 * n_point_job);
        double* d_point_y = d_point_x + n_point_job;
        double* d_point_z = d_point_y + n_point_job;
        double* d_point_w = d_point_z + n_point_job;
        int* d_point_atom_center_index = (int*)gpu->gpuAlloc(sizeof(int) * n_point_job);

        cudaMemcpy(d_point_x, h_point_x, 4 * n_point_job * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(d_point_atom_center_index, h_point_atom_center_index, n_point_job * sizeof(int), cudaMemcpyHostToDevice);

        weights((n_point_job + n_threads_per_block_weight - 1) / n_threads_per_block_weight, n_threads_per_block_weight,
                n_point_job, d_point_x, d_point_y, d_point_z, d_atoms, n_atom, d_point_atom_center_index,
                d_point_w, d_interatomic_quantities, default_switch_function_threshold, default_image_cutoff_radius, unit_cell);
        CUERR;

        cudaMemcpy(h_point_w, d_point_w, n_point_job * sizeof(double), cudaMemcpyDeviceToHost);
        CUERR;

        for (int i_global = i_thread, i_local = 0; i_global < n_point_total; i_global += n_thread, i_local++) {
            points[i_global].w_total = h_point_w[i_local];
        }
    }
}
