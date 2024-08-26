#include "periodic_becke_kernel.h"

#include "../periodic_kernel_data_cu.h"

#include "../gpubox.h"

namespace PeriodicBox
{
    static __device__ double get_r(double x1, double y1, double z1, double x2, double y2, double z2)
    {
        return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1));
    }

    static __device__ double get_one_over_r(double x1, double y1, double z1, double x2, double y2, double z2)
    {
        return rsqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1));
    }

    // This is the form proposed by Becke in the original Becke weight equation
    // J. Chem. Phys. 88, 2547-2553 (1988) https://doi.org/10.1063/1.454033
    static __device__ double switch_function(double m)
    {
        m = 1.5 * m - 0.5 * m*m*m;
        m = 1.5 * m - 0.5 * m*m*m;
        m = 1.5 * m - 0.5 * m*m*m;
        return 0.5 * (1.0 - m);
    }

    static __global__ void weights_kernel(const int n_point, const double* d_point_x, const double* d_point_y, const double* d_point_z,
                                          const float4* d_atoms, const int n_atom,
                                          const int* d_point_atom_center_index, double* d_weight, const double* d_interatomic_quantities,
                                          const double switch_function_threshold, const double image_cutoff_radius,
                                          const PeriodicKernelDataReal<double> periodic_data)
    {
        const int i_point = blockDim.x * blockIdx.x + threadIdx.x;
        if (i_point >= n_point) return;

        const double point[3] { d_point_x[i_point], d_point_y[i_point], d_point_z[i_point]};
        const int i_center_atom = d_point_atom_center_index[i_point];
        const double reference_image[3] { d_atoms[i_center_atom].x, d_atoms[i_center_atom].y, d_atoms[i_center_atom].z };

        double p_sum_b = 0.0;
        double p_a = 0.0;
        for (int i = 0; i < n_atom; i++) {
            const int i_atom_1 = (i_center_atom + i) % n_atom;
            double atom_1[3] = { d_atoms[i_atom_1].x, d_atoms[i_atom_1].y, d_atoms[i_atom_1].z };
            periodic_data.move_to_same_image(reference_image, atom_1);

            const double atom1_to_reference[3] { atom_1[0] - reference_image[0], atom_1[1] - reference_image[1], atom_1[2] - reference_image[2], };
            int image1_positive_bound[3] { 0, 0, 0 };
            int image1_negative_bound[3] { 0, 0, 0 };
            periodic_data.get_cube_bound_real(image1_positive_bound, image1_negative_bound, atom1_to_reference, image_cutoff_radius);

            for (int i_image1_x = -image1_negative_bound[0]; i_image1_x <= image1_positive_bound[0]; i_image1_x++)
                for (int i_image1_y = -image1_negative_bound[1]; i_image1_y <= image1_positive_bound[1]; i_image1_y++)
                    for (int i_image1_z = -image1_negative_bound[2]; i_image1_z <= image1_positive_bound[2]; i_image1_z++) {
                        double lattice_image_1[3];
                        periodic_data.get_absolute_coord_real(lattice_image_1, i_image1_x, i_image1_y, i_image1_z);
                        const double atom_image_1[3] = { atom_1[0] + lattice_image_1[0], atom_1[1] + lattice_image_1[1], atom_1[2] + lattice_image_1[2] };

                        const double ra = get_r(atom_image_1[0], atom_image_1[1], atom_image_1[2], point[0], point[1], point[2]);

                        double p_b = 1.0;
                        for (int i_atom_2 = 0; i_atom_2 < n_atom; i_atom_2++) {
                            const double a_ab = d_interatomic_quantities[i_atom_1 * n_atom + i_atom_2];
                            double atom_2[3] = { d_atoms[i_atom_2].x, d_atoms[i_atom_2].y, d_atoms[i_atom_2].z };
                            periodic_data.move_to_same_image(reference_image, atom_2);

                            const double atom2_to_reference[3] { atom_2[0] - reference_image[0], atom_2[1] - reference_image[1], atom_2[2] - reference_image[2], };
                            int image2_positive_bound[3] { 0, 0, 0 };
                            int image2_negative_bound[3] { 0, 0, 0 };
                            periodic_data.get_cube_bound_real(image2_positive_bound, image2_negative_bound, atom2_to_reference, image_cutoff_radius);

                            for (int i_image2_x = -image2_negative_bound[0]; i_image2_x <= image2_positive_bound[0]; i_image2_x++)
                                for (int i_image2_y = -image2_negative_bound[1]; i_image2_y <= image2_positive_bound[1]; i_image2_y++)
                                    for (int i_image2_z = -image2_negative_bound[2]; i_image2_z <= image2_positive_bound[2]; i_image2_z++) {
                                        double lattice_image_2[3];
                                        periodic_data.get_absolute_coord_real(lattice_image_2, i_image2_x, i_image2_y, i_image2_z);
                                        const double atom_image_2[3] = { atom_2[0] + lattice_image_2[0], atom_2[1] + lattice_image_2[1], atom_2[2] + lattice_image_2[2] };

                                        if ( (i_atom_1 != i_atom_2) ||
                                            !(i_image1_x == i_image2_x && i_image1_y == i_image2_y && i_image1_z == i_image2_z) ) {
                                            const double one_over_rab = get_one_over_r(atom_image_1[0], atom_image_1[1], atom_image_1[2], atom_image_2[0], atom_image_2[1], atom_image_2[2]);
                                            const double rb = get_r(atom_image_2[0], atom_image_2[1], atom_image_2[2], point[0], point[1], point[2]);
                                            const double mu = (ra - rb) * one_over_rab;
                                            // Refer to equation A2 in the original Becke paper for the next equation
                                            const double nu = mu + a_ab * (1.0 - mu * mu);
                                            p_b *= switch_function(nu);
                                            if (p_b < switch_function_threshold) {
                                                p_b = 0.0;
                                                goto jump_out_atom2_image_loop;
                                            }
                                        }
                                    }
                        }

                        jump_out_atom2_image_loop:
                        p_sum_b += p_b;
                        if (i_atom_1 == i_center_atom && i_image1_x == 0 && i_image1_y == 0 && i_image1_z == 0) {
                            p_a = p_b;
                            if (p_a == 0.0)
                                goto jump_out_atom1_image_loop;
                        }
                    }
        }

        jump_out_atom1_image_loop:
        double wt;
        if (p_a == 0.0)
            wt = 0.0;
        else
            wt = d_weight[i_point] * (p_a / p_sum_b);

        d_weight[i_point] = wt;
    }

    void weights(const int n_grid, const int n_block,
                 const int n_point, const double* d_point_x, const double* d_point_y, const double* d_point_z,
                 const float4* d_atoms, const int n_atom, const int* d_point_atom_center_index,
                 double* d_weight, const double* d_interatomic_quantities,
                 const double switch_function_threshold, const double image_cutoff_radius,
                 const LatticeVector unit_cell)
    {
        const PeriodicKernelDataReal<double> periodic_data(NAN, NAN, NAN, unit_cell);
        weights_kernel<<<n_grid, n_block>>>(n_point, d_point_x, d_point_y, d_point_z,
                                            d_atoms, n_atom, d_point_atom_center_index,
                                            d_weight, d_interatomic_quantities,
                                            switch_function_threshold, image_cutoff_radius, periodic_data);
    }

    static __global__ void interatomic_quantity_kernel(const float4* atoms, const int n_atom, double* d_interatomic_quantities)
    {
        const int i_atom_1 = threadIdx.x + blockIdx.x * blockDim.x;
        const int i_atom_2 = threadIdx.y + blockIdx.y * blockDim.y;

        if (i_atom_1 >= n_atom || i_atom_2 >= n_atom) return;

        const float4 atom_1 = atoms[i_atom_1];
        const float4 atom_2 = atoms[i_atom_2];

        double a_ab;
        // Refer to equation A3~A6 in the original Becke paper for the following equations
        const double chi = (double)atom_1.w / (double)atom_2.w;
        const double uab = (chi - 1.0) / (chi + 1.0);
        a_ab = uab / (uab * uab - 1.0);
        if (a_ab >  0.5)  a_ab =  0.5;
        if (a_ab < -0.5)  a_ab = -0.5;

        d_interatomic_quantities[i_atom_1 * n_atom + i_atom_2] = a_ab;
    }

    void set_atom_and_interatomic_quantities(GPUBox* gpu, double*& d_interatomic_quantities, float4*& d_atoms,
                                             const int n_atom, const double* atom_xyz, const double* atom_radius)
    {
        float4* h_atoms = (float4*)gpu->cpuAlloc(sizeof(float4) * n_atom);
        d_atoms = (float4*)gpu->gpuAlloc(sizeof(float4) * n_atom);

        for (int i = 0; i < n_atom; i++) {
            h_atoms[i].x = (float)atom_xyz[3 * i + 0];
            h_atoms[i].y = (float)atom_xyz[3 * i + 1];
            h_atoms[i].z = (float)atom_xyz[3 * i + 2];
            h_atoms[i].w = (float)atom_radius[i];
        }

        cudaMemcpy(d_atoms, h_atoms, n_atom * sizeof(float4), cudaMemcpyHostToDevice);
        gpu->cpuFree(h_atoms);

        const dim3 block_dimension(16, 16);
        const dim3 grid_dimension((n_atom + block_dimension.x - 1) / block_dimension.x, (n_atom + block_dimension.y - 1) / block_dimension.y, 1);

        d_interatomic_quantities = (double*)gpu->gpuAlloc(sizeof(double) * n_atom * n_atom);
        interatomic_quantity_kernel<<<grid_dimension, block_dimension>>>(d_atoms, n_atom, d_interatomic_quantities);
    }

    // For gradient

//     __global__ void Wgrad_cache(const int numPoints, const int numAtoms, const int grad_pitch,
//                                 const double* xpts, const double* ypts, const double* zpts, const float4* d_atoms,
//                                 double* ra_cache, double* p_cache,
//                                 const double2* cache, const double switch_function_threshold)
//     {
//         int ptIndex = threadIdx.x+blockIdx.x*WGRAD_BLOCK_XDIM_DP;
//         if( ptIndex >= numPoints )    return;
//         int a = threadIdx.y + blockIdx.y*WGRAD_BLOCK_YDIM_DP;
//         if( a >= numAtoms )    return;
//         double ptx = xpts[ptIndex];
//         double pty = ypts[ptIndex];
//         double ptz = zpts[ptIndex];
//         float4 temp = d_atoms[a];

//         double ra = get_r(temp.x, temp.y, temp.z, ptx, pty, ptz);
//         ra_cache[ptIndex+grad_pitch*a] = ra;

//         double pVal = 1.0;
//         double2 buff = cache[numAtoms*a];

//         for(int b=0; b<numAtoms; ++b)
//         {
//             double2 reg = buff;
//             int next = b+1;
//             buff = cache[a*numAtoms+next];

//             if( a == b )    continue;

//             temp = d_atoms[b];
//             double rb = get_r(temp.x, temp.y, temp.z, ptx, pty, ptz);
//             double mu = (ra - rb) * reg.x;
//             double nu = mu + reg.y*(1.0 - mu*mu);
//             pVal *= switch_function(nu);
//             if( fabs(pVal) < switch_function_threshold ){
//                 pVal = 0.0;
//                 break;
//             }
//         }
//         p_cache[ptIndex+grad_pitch*a] = pVal;
//     }

//     void Wgrad_cache_caller(const dim3 GRID, const dim3 BLOCK, const int numPoints, const int numAtoms, const int grad_pitch,
//                             const double* xpts, const double* ypts, const double* zpts, const float4* d_atoms,
//                             double* ra_cache, double* p_cache,
//                             const double2* cache, const double switch_function_threshold)
//     {
//         Wgrad_cache<<<GRID, BLOCK>>>(numPoints, numAtoms, grad_pitch, xpts, ypts, zpts, d_atoms, ra_cache, p_cache, cache, switch_function_threshold);
//     }

//     __device__ double3 Dev_grad_nu(int g, double3 rg_vec, double rg, int b, double rb, double d_gb, double a_gb, const float4* d_atoms)
//     {
//         double3 grad;

//         double3 r_gb;
//         r_gb.x = d_atoms[g].x - d_atoms[b].x;
//         r_gb.y = d_atoms[g].y - d_atoms[b].y;
//         r_gb.z = d_atoms[g].z - d_atoms[b].z;

//         double mu = (rg - rb)*d_gb;
//         double coef = 1.0-2.0*a_gb*mu;
//         double r_g_coef;

//         if (rg < 1.0e-14)     r_g_coef = 0.0;
//         else                  r_g_coef = -coef*d_gb/rg;
//         grad.x = r_g_coef * rg_vec.x;
//         grad.y = r_g_coef * rg_vec.y;
//         grad.z = r_g_coef * rg_vec.z;

//         double r_gb_coef = -coef*mu*d_gb*d_gb;
//         grad.x += r_gb_coef * r_gb.x;
//         grad.y += r_gb_coef * r_gb.y;
//         grad.z += r_gb_coef * r_gb.z;
//         return grad;
//     }

    static __device__ double switch_function_derivative(const double m)
    {
        const double f1 = 1.5 *  m - 0.5 * m*m*m;
        const double f2 = 1.5 * f1 - 0.5 * f1*f1*f1;
        const double f3 = 1.5 * f2 - 0.5 * f2*f2*f2;
        const double s = 0.5 * (1.0 - f3);
        if (fabs(s) < 1.0e-14)
            return 0.0;
        else
            return -(27.0/16.0) * (1.0 - f2*f2) * (1.0 - f1*f1) * (1.0 - m*m) / s;
    }
//     __device__ double Dev_comp_t(double mu_ij, double a_ij)
//     {
//         double t;
//         double nu = mu_ij + a_ij*(1.0-mu_ij*mu_ij);
//         double p1 = 1.5*nu - 0.5*nu*nu*nu;
//         double p2 = 1.5*p1 - 0.5*p1*p1*p1;
//         double s_val;
//         double p3 = 1.5*p2 - 0.5*p2*p2*p2;
//         s_val = 0.5 * (1.0 - p3);
//         if (fabs(s_val) < 1.0e-14)
//             return 0.0;
//         t = -(27.0/16.0) * (1.0 - p2*p2) * (1.0 - p1*p1) * (1.0 - nu*nu) / s_val;
//         return t;
//     }

//     __global__ void Wgrad_kernel(const int numPoints, const int numAtoms, const int grad_pitch,
//                                 const double* xpts, const double* ypts, const double* zpts, const double* wpts,
//                                 const float4* d_atoms, const int* ctrs,
//                                 double* xGrad, double* yGrad, double* zGrad,
//                                 const double* ra_cache, const double* p_cache, const double2* cache)
//     {
//         int ptIndex = threadIdx.x+blockIdx.x*WGRAD_BLOCK_XDIM_DP;
//         if( ptIndex >= numPoints )    return;
//         int ctr = ctrs[ptIndex];
//         int grad_atom = threadIdx.y + blockIdx.y*WGRAD_BLOCK_YDIM_DP;
//         if( grad_atom >= ctr )     ++grad_atom;
//         if( grad_atom >= numAtoms )    return;
//         double ptx = xpts[ptIndex];
//         double pty = ypts[ptIndex];
//         double ptz = zpts[ptIndex];
//         double ptw = wpts[ptIndex];
//         float4 tAtom;

//         double p_buff = p_cache[ptIndex];
//         double p_sum = 0.0;
//         double p_a = p_cache[ptIndex + ctr       * grad_pitch];
//         double p_g = p_cache[ptIndex + grad_atom * grad_pitch];
//         double3 grad_c_w_a;    grad_c_w_a.x = grad_c_w_a.y = grad_c_w_a.z = 0.0f;
//         double3 grad_c_p_a;    grad_c_p_a.x = grad_c_p_a.y = grad_c_p_a.z = 0.0f;
//         double3 grad_c_p_g;    grad_c_p_g.x = grad_c_p_g.y = grad_c_p_g.z = 0.0f;
//         tAtom = d_atoms[grad_atom];
//         double3 rg_vec;
//         rg_vec.x = ptx - tAtom.x;
//         rg_vec.y = pty - tAtom.y;
//         rg_vec.z = ptz - tAtom.z;
//         double rg = sqrt(rg_vec.x*rg_vec.x + rg_vec.y*rg_vec.y + rg_vec.z*rg_vec.z);

//         int next = ptIndex;
//         for (int bcenter=0; bcenter<numAtoms; ++bcenter)
//         {
//             double3 grad_c_p_b;
//             double t;

//             next += grad_pitch;
//             double p_b = p_buff;     p_buff = p_cache[next];
//             p_sum += p_b;
//             if( (p_b == 0.0f && p_g == 0.0f) || grad_atom == bcenter )
//                 continue;

//             double rb  = ra_cache[ptIndex+bcenter*grad_pitch];
//             double2 bg = cache[bcenter*numAtoms+grad_atom];
//             double mu_gb = (rg - rb) * bg.x;
//             double3 grad_nu = Dev_grad_nu(grad_atom, rg_vec, rg, bcenter, rb, bg.x, -bg.y, d_atoms);

//             if( p_g != 0.0f )
//             {
//                 t = Dev_comp_t(mu_gb, -bg.y);
//                 grad_c_p_g.x += t*grad_nu.x;
//                 grad_c_p_g.y += t*grad_nu.y;
//                 grad_c_p_g.z += t*grad_nu.z;
//             }

//             if(p_b != 0.0f)
//             {
//                 t = -Dev_comp_t(-mu_gb, bg.y)*p_b;
//                 grad_c_p_b.x = t * grad_nu.x;
//                 grad_c_p_b.y = t * grad_nu.y;
//                 grad_c_p_b.z = t * grad_nu.z;

//                 if(bcenter == ctr)
//                     grad_c_p_a = grad_c_p_b;

//                 grad_c_w_a.x -= grad_c_p_b.x;
//                 grad_c_w_a.y -= grad_c_p_b.y;
//                 grad_c_w_a.z -= grad_c_p_b.z;
//             }
//         }
//         grad_c_w_a.x -= grad_c_p_g.x*p_g;
//         grad_c_w_a.y -= grad_c_p_g.y*p_g;
//         grad_c_w_a.z -= grad_c_p_g.z*p_g;

//         double coef1 = p_a / p_sum;
//         double coef2 = ptw / p_sum;
//         xGrad[ptIndex+grad_atom*grad_pitch] += coef2*(coef1*grad_c_w_a.x + grad_c_p_a.x);
//         yGrad[ptIndex+grad_atom*grad_pitch] += coef2*(coef1*grad_c_w_a.y + grad_c_p_a.y);
//         zGrad[ptIndex+grad_atom*grad_pitch] += coef2*(coef1*grad_c_w_a.z + grad_c_p_a.z);
//     }

//     void Wgrad_kernel_caller(const dim3 GRID, const dim3 BLOCK, const int numPoints, const int numAtoms, const int grad_pitch, 
//                             const double* xpts, const double* ypts, const double* zpts, const double* wpts,
//                             const float4* d_atoms, const int* ctrs,
//                             double* xGrad, double* yGrad, double* zGrad,
//                             const double* ra_cache, const double* p_cache, const double2* cache)
//     {
//         Wgrad_kernel<<<GRID, BLOCK>>>(numPoints, numAtoms, grad_pitch, xpts, ypts, zpts, wpts, d_atoms,
//                                     ctrs, xGrad, yGrad, zGrad, ra_cache, p_cache, cache);
//     }

//     __global__ void Wgrad_ptctr(const int numPoints, const int numAtoms, const int grad_pitch, const int* ctrs,
//                                 double* xGrad, double* yGrad, double* zGrad)
//     {
//         double x=0.0, y=0.0, z=0.0;
//         int ptIndex = threadIdx.x + blockIdx.x*WGRAD_PTCTR_BLOCK;
//         if(ptIndex >= numPoints)    return;
//         int ctr = ctrs[ptIndex];

//         for(int j=0; j<numAtoms; ++j){
//             x -= xGrad[ptIndex+j*grad_pitch];
//             y -= yGrad[ptIndex+j*grad_pitch];
//             z -= zGrad[ptIndex+j*grad_pitch];
//         }

//         xGrad[ptIndex+ctr*grad_pitch] += x;
//         yGrad[ptIndex+ctr*grad_pitch] += y;
//         zGrad[ptIndex+ctr*grad_pitch] += z;
//     }

//     void Wgrad_ptctr_caller(const int GRID, const int BLOCK, const int numPoints, const int numAtoms, const int grad_pitch, const int* ctrs,
//                             double* xGrad, double* yGrad, double* zGrad)
//     {
//         Wgrad_ptctr<<<GRID, BLOCK>>>(numPoints, numAtoms, grad_pitch, ctrs, xGrad, yGrad, zGrad);
//     }

//     __global__ void Wgrad_dosums(const int numPoints, double* final_grad,
//                                 const double* xGrad, const double* yGrad, const double* zGrad)
//     {
//         __shared__ double scr[WGRAD_DOSUMS_BLOCK];
//         int atom = blockIdx.x/3;
//         int dimSelect = blockIdx.x%3;
//         const double *end, *cur, *temp;
//         double sum = 0.0;
//         if( dimSelect == 0 )
//             temp = xGrad;
//         else if( dimSelect == 1 )
//             temp = yGrad;
//         else
//             temp = zGrad;
//         cur = temp + numPoints*atom + threadIdx.x;
//         end = temp + numPoints * (atom+1);

//         while(cur < end)
//         {
//             sum += *cur;
//             cur += WGRAD_DOSUMS_BLOCK;
//         }
//         scr[threadIdx.x] = sum;

//         for(int span=WGRAD_DOSUMS_BLOCK/2; span>1; span /= 2)
//         {
//             __syncthreads();
//             if(threadIdx.x < span)
//                 scr[threadIdx.x] += scr[threadIdx.x+span];
//         }

//         if( !threadIdx.x )
//         {
//             sum = scr[0]+scr[1];
//             final_grad[blockIdx.x] = sum;
//         }
//     }

//     void Wgrad_dosums_caller(const int GRID, const int BLOCK, const int numPoints, double* final_grad,
//                             const double* xGrad, const double* yGrad, const double* zGrad)
//     {
//         Wgrad_dosums<<<GRID, BLOCK>>>(numPoints, final_grad, xGrad, yGrad, zGrad);
//     }

}
