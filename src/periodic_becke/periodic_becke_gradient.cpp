
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

    //     int bytes_per_pt = 2*(n_atom+1)*sizeof(double) + 3*n_atom*sizeof(double);
    //     int maxPts = remaining_gpu_memory_bytes / bytes_per_pt;
    //     int pts_per_grid;
    //     if (n_point_job < maxPts)
    //         pts_per_grid = n_point_job;
    //     else
    //         pts_per_grid = (maxPts / WGRAD_BLOCK_XDIM_DP) * WGRAD_BLOCK_XDIM_DP;
    //     pts_per_grid = min(pts_per_grid, MAX_WT_GRID_DIM*WGRAD_BLOCK_XDIM_DP);

    //     double*  d_rval = (double*)gpu->gpuAlloc(sizeof(double)*pts_per_grid*(n_atom+1));
    //     double*  d_pval = (double*)gpu->gpuAlloc(sizeof(double)*pts_per_grid*(n_atom+1));
    //     double* d_xGrad = (double*)gpu->gpuAlloc(sizeof(double)*pts_per_grid*n_atom);
    //     double* d_yGrad = (double*)gpu->gpuAlloc(sizeof(double)*pts_per_grid*n_atom);
    //     double* d_zGrad = (double*)gpu->gpuAlloc(sizeof(double)*pts_per_grid*n_atom);
    //     cudaMemset(d_xGrad, 0, sizeof(double)*n_atom*pts_per_grid); CUERR;
    //     cudaMemset(d_yGrad, 0, sizeof(double)*n_atom*pts_per_grid); CUERR;
    //     cudaMemset(d_zGrad, 0, sizeof(double)*n_atom*pts_per_grid); CUERR;

    //     dim3 grid_dim;
    //     grid_dim.y = n_atom / WGRAD_BLOCK_YDIM_DP;
    //     if (n_atom % WGRAD_BLOCK_YDIM_DP)
    //         grid_dim.y++;
    //     dim3 block_dim(WGRAD_BLOCK_XDIM_DP, WGRAD_BLOCK_YDIM_DP);

    //     for (int j=0; j<n_point_job; j+=pts_per_grid) {
    //         int npts = min(pts_per_grid, n_point_job-j);
    //         grid_dim.x = npts / WGRAD_BLOCK_XDIM_DP;
    //         if( npts % WGRAD_BLOCK_XDIM_DP )
    //             grid_dim.x++;
    //         int ptctr_grid_dim = npts / WGRAD_PTCTR_BLOCK;
    //         if( npts % WGRAD_PTCTR_BLOCK )
    //             ptctr_grid_dim++;
    //         cudaThreadSynchronize(); CUERR;

    //         Wgrad_cache_caller(grid_dim, block_dim, npts, n_atom, pts_per_grid,
    //                         d_point_x+j, d_point_y+j, d_point_z+j, d_atoms,
    //                         d_rval, d_pval, d_interatomic_quantities, W_THRE); CUERR;

    //         Wgrad_kernel_caller(grid_dim, block_dim, npts, n_atom, pts_per_grid,
    //                             d_point_x+j, d_point_y+j, d_point_z+j, d_point_w+j, d_atoms,
    //                             d_i_atom_for_point+j, d_xGrad, d_yGrad, d_zGrad,
    //                             d_rval, d_pval, d_interatomic_quantities); CUERR;

    //         Wgrad_ptctr_caller(ptctr_grid_dim, WGRAD_PTCTR_BLOCK, npts, n_atom,
    //                         pts_per_grid, d_i_atom_for_point+j, d_xGrad, d_yGrad, d_zGrad); CUERR;
    //     }

    //     double* h_grads = (double*)gpu->cpuAlloc(sizeof(double)*n_atom*3);
    //     double* h_xGrad = (double*)gpu->cpuAlloc(sizeof(double)*n_atom*pts_per_grid);
    //     double* h_yGrad = (double*)gpu->cpuAlloc(sizeof(double)*n_atom*pts_per_grid);
    //     double* h_zGrad = (double*)gpu->cpuAlloc(sizeof(double)*n_atom*pts_per_grid);
    //     cudaMemcpy(h_xGrad, d_xGrad, sizeof(double)*n_atom*pts_per_grid, cudaMemcpyDeviceToHost); CUERR;
    //     cudaMemcpy(h_yGrad, d_yGrad, sizeof(double)*n_atom*pts_per_grid, cudaMemcpyDeviceToHost); CUERR;
    //     cudaMemcpy(h_zGrad, d_zGrad, sizeof(double)*n_atom*pts_per_grid, cudaMemcpyDeviceToHost); CUERR;

    //     double* d_grads = (double*)gpu->gpuAlloc(sizeof(double)*n_atom*3);
    //     int sums_grid_dim = n_atom*3;
    //     Wgrad_dosums_caller(sums_grid_dim, WGRAD_DOSUMS_BLOCK, pts_per_grid, d_grads,
    //                         d_xGrad, d_yGrad, d_zGrad); CUERR;

    //     cudaMemcpy(h_grads, d_grads, sizeof(double)*n_atom*3, cudaMemcpyDeviceToHost); CUERR;
    //     for (int j=0; j<n_atom*3; ++j)
    //         gradient[j] += (double)h_grads[j];

    //     gpu->gpuFree(d_grads);
    //     gpu->cpuFree(h_zGrad);
    //     gpu->cpuFree(h_yGrad);
    //     gpu->cpuFree(h_xGrad);
    //     gpu->cpuFree(h_grads);
    //     gpu->gpuFree(d_zGrad);
    //     gpu->gpuFree(d_yGrad);
    //     gpu->gpuFree(d_xGrad);
    //     gpu->gpuFree(d_pval);
    //     gpu->gpuFree(d_rval);
        gpu->gpuFree(d_i_atom_for_point);
        gpu->gpuFree(d_point_xyzw);
        gpu->cpuFree(h_i_atom_for_point);
        gpu->cpuFree(h_point_xyzw);
        gpu->gpuFree(d_interatomic_quantities);
        gpu->gpuFree(d_atoms);
    }
}
