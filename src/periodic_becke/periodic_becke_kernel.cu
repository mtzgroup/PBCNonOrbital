#include "periodic_becke_kernel.h"

#include "../periodic_kernel_data_cu.h"

#include "../gpubox.h"

#define DIVIDE_BY_ZERO_PROTECTION_THRESHOLD 1.0e-14

namespace PeriodicBox
{
    static __device__ double get_r(double x1, double y1, double z1, double x2, double y2, double z2)
    {
        return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1));
    }

    static __device__ double get_one_over_r(double x1, double y1, double z1, double x2, double y2, double z2)
    {
        const double r = get_r(x1, y1, z1, x2, y2, z2);
        return (r > DIVIDE_BY_ZERO_PROTECTION_THRESHOLD) ? (1.0 / r) : 0.0;
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

            const double atom1_to_reference[3] { atom_1[0] - reference_image[0], atom_1[1] - reference_image[1], atom_1[2] - reference_image[2] };
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

                            const double atom2_to_reference[3] { atom_2[0] - reference_image[0], atom_2[1] - reference_image[1], atom_2[2] - reference_image[2] };
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

    __device__ double3 smooth_function_dmudA(const double A[3], const double B[3], const double point[3], const double a_ab, const bool A_equal_B)
    {
        const double Ar[3] { A[0] - point[0], A[1] - point[1], A[2] - point[2] };
        const double Br[3] { B[0] - point[0], B[1] - point[1], B[2] - point[2] };
        const double AB[3] { A[0] - B[0], A[1] - B[1], A[2] - B[2] };
        const double norm_Ar = sqrt(Ar[0] * Ar[0] + Ar[1] * Ar[1] + Ar[2] * Ar[2]);
        const double norm_Br = sqrt(Br[0] * Br[0] + Br[1] * Br[1] + Br[2] * Br[2]);
        const double one_over_Ar = (norm_Ar > DIVIDE_BY_ZERO_PROTECTION_THRESHOLD) ? (1.0 / norm_Ar) : 0.0;
        const double one_over_AB = get_one_over_r(A[0], A[1], A[2], B[0], B[1], B[2]);
        const double normalized_Ar[3] = { Ar[0] * one_over_Ar, Ar[1] * one_over_Ar, Ar[2] * one_over_Ar };
        const double mu = (norm_Ar - norm_Br) * one_over_AB;
        if (!A_equal_B) {
            const double normalized_AB[3] = { AB[0] * one_over_AB, AB[1] * one_over_AB, AB[2] * one_over_AB };
            double3 gradient;
            gradient.x = one_over_AB * (normalized_Ar[0] - mu * normalized_AB[0]);
            gradient.y = one_over_AB * (normalized_Ar[1] - mu * normalized_AB[1]);
            gradient.z = one_over_AB * (normalized_Ar[2] - mu * normalized_AB[2]);
            gradient.x *= (1.0 - 2.0 * a_ab * mu); // TODO: Check if this is correct with different elements
            gradient.y *= (1.0 - 2.0 * a_ab * mu); // TODO: Check if this is correct with different elements
            gradient.z *= (1.0 - 2.0 * a_ab * mu); // TODO: Check if this is correct with different elements
            return gradient;
        } else {
            const double one_over_Br = get_one_over_r(B[0], B[1], B[2], point[0], point[1], point[2]);
            const double normalized_Br[3] = { Br[0] * one_over_Br, Br[1] * one_over_Br, Br[2] * one_over_Br };
            double3 gradient;
            gradient.x = one_over_AB * (normalized_Ar[0] - normalized_Br[0]);
            gradient.y = one_over_AB * (normalized_Ar[1] - normalized_Br[1]);
            gradient.z = one_over_AB * (normalized_Ar[2] - normalized_Br[2]);
            gradient.x *= (1.0 - 2.0 * a_ab * mu); // TODO: Check if this is correct with different elements
            gradient.y *= (1.0 - 2.0 * a_ab * mu); // TODO: Check if this is correct with different elements
            gradient.z *= (1.0 - 2.0 * a_ab * mu); // TODO: Check if this is correct with different elements
            return gradient;
        }
    }

    __device__ double3 smooth_function_dmudB(const double A[3], const double B[3], const double point[3], const double a_ab, const bool A_equal_B)
    {
        return smooth_function_dmudA(B, A, point, a_ab, A_equal_B);
    }

    static __device__ double switch_function_dsdx_over_s(const double m)
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

    static __device__ double switch_function_dsdmu_over_s(const double mu, const double a_ab)
    {
        const double nu = mu + a_ab * (1.0 - mu * mu);
        const double dnudmu = 1.0 - 2.0 * a_ab * mu;
        return switch_function_dsdx_over_s(nu) * dnudmu;
    }

    // Note: This value cannot be cached, because atom_a can be any image of i_atom_a.
    static __device__ double compute_P_A(const double* atom_a, const int i_atom_a, const double* point,
                                         const int n_atom,
                                         const float4* d_atoms,
                                         const double* d_interatomic_quantities, const double switch_function_threshold,
                                         const double image_cutoff_radius, const PeriodicKernelDataReal<double> periodic_data)
    {
        const double ra = get_r(atom_a[0], atom_a[1], atom_a[2], point[0], point[1], point[2]);

        double p_a = 1.0;
        for (int i_atom_b = 0; i_atom_b < n_atom; i_atom_b++) {
            const double a_ab = d_interatomic_quantities[i_atom_a * n_atom + i_atom_b];
            double atom_b[3] = { d_atoms[i_atom_b].x, d_atoms[i_atom_b].y, d_atoms[i_atom_b].z };
            periodic_data.move_to_same_image(atom_a, atom_b);

            const double atom_b_to_reference[3] { atom_b[0] - atom_a[0], atom_b[1] - atom_a[1], atom_b[2] - atom_a[2] };
            int image2_positive_bound[3] { 0, 0, 0 };
            int image2_negative_bound[3] { 0, 0, 0 };
            periodic_data.get_cube_bound_real(image2_positive_bound, image2_negative_bound, atom_b_to_reference, image_cutoff_radius);

            for (int i_image2_x = -image2_negative_bound[0]; i_image2_x <= image2_positive_bound[0]; i_image2_x++)
                for (int i_image2_y = -image2_negative_bound[1]; i_image2_y <= image2_positive_bound[1]; i_image2_y++)
                    for (int i_image2_z = -image2_negative_bound[2]; i_image2_z <= image2_positive_bound[2]; i_image2_z++) {
                        double lattice_image_b[3];
                        periodic_data.get_absolute_coord_real(lattice_image_b, i_image2_x, i_image2_y, i_image2_z);
                        const double atom_image_b[3] = { atom_b[0] + lattice_image_b[0], atom_b[1] + lattice_image_b[1], atom_b[2] + lattice_image_b[2] };

                        if ( (i_atom_a != i_atom_b) ||
                            !(i_image2_x == 0 && i_image2_y == 0 && i_image2_z == 0) ) {
                            const double one_over_rab = get_one_over_r(atom_a[0], atom_a[1], atom_a[2], atom_image_b[0], atom_image_b[1], atom_image_b[2]);
                            const double rb = get_r(atom_image_b[0], atom_image_b[1], atom_image_b[2], point[0], point[1], point[2]);
                            const double mu = (ra - rb) * one_over_rab;
                            // Refer to equation A2 in the original Becke paper for the next equation
                            const double nu = mu + a_ab * (1.0 - mu * mu);
                            p_a *= switch_function(nu);
                            if (p_a < switch_function_threshold) {
                                return 0.0;
                            }
                        }
                    }
        }

        return p_a;
    }

    __global__ void weight_gradient_compute_kernel(const int n_point, const int n_atom, const int n_point_per_grid,
                                                   const double* d_point_x, const double* d_point_y, const double* d_point_z, const double* d_point_w,
                                                   const float4* d_atoms, const int* d_i_atom_for_point,
                                                   double* d_gradient_cache_x, double* d_gradient_cache_y, double* d_gradient_cache_z,
                                                   const double* d_interatomic_quantities, const double switch_function_threshold,
                                                   const double image_cutoff_radius, const PeriodicKernelDataReal<double> periodic_data)
    {
        const int i_point = threadIdx.x + blockIdx.x * blockDim.x;
        if (i_point >= n_point) return;
        const int i_center_atom = d_i_atom_for_point[i_point]; // Called A

        const int i_derivative_atom = threadIdx.y + blockIdx.y * blockDim.y; // Called G
        if (i_derivative_atom == i_center_atom) return;
        if (i_derivative_atom >= n_atom) return;

        const double point[3] { d_point_x[i_point], d_point_y[i_point], d_point_z[i_point]};
        const double atom_a[3] { d_atoms[i_center_atom].x, d_atoms[i_center_atom].y, d_atoms[i_center_atom].z };
        double atom_g[3] { d_atoms[i_derivative_atom].x, d_atoms[i_derivative_atom].y, d_atoms[i_derivative_atom].z };
        periodic_data.move_to_same_image(atom_a, atom_g);

        // $$\frac{\partial P^{PBC}(\vec{A}, \vec{r})}{\partial \vec{G}}$$

        const double rA = get_r(atom_a[0], atom_a[1], atom_a[2], point[0], point[1], point[2]);
        const double a_AG = d_interatomic_quantities[i_center_atom * n_atom + i_derivative_atom];

        const double AG_without_offset[3] { atom_g[0] - atom_a[0], atom_g[1] - atom_a[1], atom_g[2] - atom_a[2] };
        int image1_positive_bound[3] { 0, 0, 0 };
        int image1_negative_bound[3] { 0, 0, 0 };
        periodic_data.get_cube_bound_real(image1_positive_bound, image1_negative_bound, AG_without_offset, image_cutoff_radius);

        const double P_A = compute_P_A(atom_a, i_center_atom, point, n_atom, d_atoms,
                                       d_interatomic_quantities, switch_function_threshold, image_cutoff_radius, periodic_data);
        if (P_A < switch_function_threshold) {
            d_gradient_cache_x[i_point + i_derivative_atom * n_point_per_grid] = 0.0;
            d_gradient_cache_y[i_point + i_derivative_atom * n_point_per_grid] = 0.0;
            d_gradient_cache_z[i_point + i_derivative_atom * n_point_per_grid] = 0.0;
            return;
        }

        double dP_A_dG[3] { 0.0, 0.0, 0.0 };
        for (int i_image1_x = -image1_negative_bound[0]; i_image1_x <= image1_positive_bound[0]; i_image1_x++)
            for (int i_image1_y = -image1_negative_bound[1]; i_image1_y <= image1_positive_bound[1]; i_image1_y++)
                for (int i_image1_z = -image1_negative_bound[2]; i_image1_z <= image1_positive_bound[2]; i_image1_z++) {

                    double lattice_image_g[3];
                    periodic_data.get_absolute_coord_real(lattice_image_g, i_image1_x, i_image1_y, i_image1_z);
                    const double atom_g_image[3] = { atom_g[0] + lattice_image_g[0], atom_g[1] + lattice_image_g[1], atom_g[2] + lattice_image_g[2] };

                    const double one_over_AG = get_one_over_r(atom_a[0], atom_a[1], atom_a[2], atom_g_image[0], atom_g_image[1], atom_g_image[2]);
                    const double rG = get_r(atom_g_image[0], atom_g_image[1], atom_g_image[2], point[0], point[1], point[2]);
                    const double mu = (rA - rG) * one_over_AG;

                    const double dsdmu_over_s = switch_function_dsdmu_over_s(mu, a_AG);
                    const double3 dmudG = smooth_function_dmudB(atom_a, atom_g_image, point, a_AG, false);

                    dP_A_dG[0] += dsdmu_over_s * dmudG.x;
                    dP_A_dG[1] += dsdmu_over_s * dmudG.y;
                    dP_A_dG[2] += dsdmu_over_s * dmudG.z;
                }

        dP_A_dG[0] *= P_A;
        dP_A_dG[1] *= P_A;
        dP_A_dG[2] *= P_A;

        // $$\sum_{B}^{N_{atom}} \sum_{\vec{P}_3 \in Z^3} \frac{\partial P^{PBC}(\vec{B} + \vec{P}_3, \vec{r}) }{\partial \vec{G}}$$

        const double P_G = compute_P_A(atom_g, i_derivative_atom, point, n_atom, d_atoms,
            d_interatomic_quantities, switch_function_threshold, image_cutoff_radius, periodic_data);
        double dP_G_dG[3] { 0.0, 0.0, 0.0 };

        double dP_B_dG_sum[3] { 0.0, 0.0, 0.0 };
        double P_B_sum = 0.0;
        for (int i_atom_b = 0; i_atom_b < n_atom; i_atom_b++) {
            double atom_b[3] { d_atoms[i_atom_b].x, d_atoms[i_atom_b].y, d_atoms[i_atom_b].z };
            periodic_data.move_to_same_image(atom_a, atom_b);

            const double a_BG = d_interatomic_quantities[i_atom_b * n_atom + i_derivative_atom];

            const double BG_without_offset[3] { atom_g[0] - atom_b[0], atom_g[1] - atom_b[1], atom_g[2] - atom_b[2] };
            int image3_positive_bound[3] { 0, 0, 0 };
            int image3_negative_bound[3] { 0, 0, 0 };
            periodic_data.get_cube_bound_real(image3_positive_bound, image3_negative_bound, BG_without_offset, image_cutoff_radius);

            for (int i_image3_x = -image3_negative_bound[0]; i_image3_x <= image3_positive_bound[0]; i_image3_x++)
                for (int i_image3_y = -image3_negative_bound[1]; i_image3_y <= image3_positive_bound[1]; i_image3_y++)
                    for (int i_image3_z = -image3_negative_bound[2]; i_image3_z <= image3_positive_bound[2]; i_image3_z++) {
                        double lattice_image_b[3];
                        periodic_data.get_absolute_coord_real(lattice_image_b, i_image3_x, i_image3_y, i_image3_z);
                        const double atom_b_image[3] = { atom_b[0] + lattice_image_b[0], atom_b[1] + lattice_image_b[1], atom_b[2] + lattice_image_b[2] };

                        const double rB = get_r(atom_b_image[0], atom_b_image[1], atom_b_image[2], point[0], point[1], point[2]);

                        if (P_G >= switch_function_threshold) {
                            const double one_over_BG = get_one_over_r(atom_b_image[0], atom_b_image[1], atom_b_image[2], atom_g[0], atom_g[1], atom_g[2]);
                            const double rG = get_r(atom_g[0], atom_g[1], atom_g[2], point[0], point[1], point[2]);
                            const double mu = (rG - rB) * one_over_BG;

                            const double dsdmu_over_s = switch_function_dsdmu_over_s(mu, a_BG);
                            const double3 dmuGBdG = smooth_function_dmudB(atom_g, atom_b_image, point, a_BG, i_atom_b == i_derivative_atom);
                            if (!(i_atom_b == i_derivative_atom && i_image3_x == 0 && i_image3_y == 0 && i_image3_z == 0)) {
                                dP_G_dG[0] += dsdmu_over_s * dmuGBdG.x;
                                dP_G_dG[1] += dsdmu_over_s * dmuGBdG.y;
                                dP_G_dG[2] += dsdmu_over_s * dmuGBdG.z;
                            }
                        }

                        const double P_B = compute_P_A(atom_b_image, i_atom_b, point, n_atom, d_atoms,
                                                       d_interatomic_quantities, switch_function_threshold, image_cutoff_radius, periodic_data);
                        if (P_B < switch_function_threshold)
                            continue;

                        double dP_B_dG[3] { 0.0, 0.0, 0.0 };
                        for (int i_image2_x = -image3_negative_bound[0]; i_image2_x <= image3_positive_bound[0]; i_image2_x++)
                            for (int i_image2_y = -image3_negative_bound[1]; i_image2_y <= image3_positive_bound[1]; i_image2_y++)
                                for (int i_image2_z = -image3_negative_bound[2]; i_image2_z <= image3_positive_bound[2]; i_image2_z++) {
                                    double lattice_image_g[3];
                                    periodic_data.get_absolute_coord_real(lattice_image_g, i_image2_x, i_image2_y, i_image2_z);
                                    const double atom_g_image[3] = { atom_g[0] + lattice_image_g[0], atom_g[1] + lattice_image_g[1], atom_g[2] + lattice_image_g[2] };

                                    const double one_over_BG = get_one_over_r(atom_b_image[0], atom_b_image[1], atom_b_image[2], atom_g_image[0], atom_g_image[1], atom_g_image[2]);
                                    const double rG = get_r(atom_g_image[0], atom_g_image[1], atom_g_image[2], point[0], point[1], point[2]);
                                    const double mu = (rB - rG) * one_over_BG;

                                    const double dsdmu_over_s = switch_function_dsdmu_over_s(mu, a_BG);
                                    const double3 dmuBGdG = smooth_function_dmudB(atom_b_image, atom_g_image, point, a_BG, false);

                                    if (i_atom_b != i_derivative_atom) {
                                        dP_B_dG[0] += dsdmu_over_s * dmuBGdG.x;
                                        dP_B_dG[1] += dsdmu_over_s * dmuBGdG.y;
                                        dP_B_dG[2] += dsdmu_over_s * dmuBGdG.z;
                                    }
                                }

                        dP_B_dG[0] *= P_B;
                        dP_B_dG[1] *= P_B;
                        dP_B_dG[2] *= P_B;
                        dP_B_dG_sum[0] += dP_B_dG[0];
                        dP_B_dG_sum[1] += dP_B_dG[1];
                        dP_B_dG_sum[2] += dP_B_dG[2];
                        P_B_sum += P_B;
                    }
        }

        dP_B_dG_sum[0] += dP_G_dG[0] * P_G;
        dP_B_dG_sum[1] += dP_G_dG[1] * P_G;
        dP_B_dG_sum[2] += dP_G_dG[2] * P_G;

        // Combine the two pieces

        const double point_weight = d_point_w[i_point];
        // Henry 20250302: I don't understand the overall negative sign.
        d_gradient_cache_x[i_point + i_derivative_atom * n_point_per_grid] += - point_weight / P_B_sum * (dP_A_dG[0] - P_A / P_B_sum * dP_B_dG_sum[0]);
        d_gradient_cache_y[i_point + i_derivative_atom * n_point_per_grid] += - point_weight / P_B_sum * (dP_A_dG[1] - P_A / P_B_sum * dP_B_dG_sum[1]);
        d_gradient_cache_z[i_point + i_derivative_atom * n_point_per_grid] += - point_weight / P_B_sum * (dP_A_dG[2] - P_A / P_B_sum * dP_B_dG_sum[2]);
    }

    void weight_gradient_compute(const dim3 n_grid, const dim3 n_block, const int n_point, const int n_atom, const int n_point_per_grid,
                                 const double* d_point_x, const double* d_point_y, const double* d_point_z, const double* d_point_w,
                                 const float4* d_atoms, const int* d_i_atom_for_point,
                                 double* d_gradient_cache_x, double* d_gradient_cache_y, double* d_gradient_cache_z,
                                 const double* d_interatomic_quantities, const double switch_function_threshold,
                                 const double image_cutoff_radius, const LatticeVector unit_cell)
    {
        const PeriodicKernelDataReal<double> periodic_data(NAN, NAN, NAN, unit_cell);
        weight_gradient_compute_kernel<<<n_grid, n_block>>>(n_point, n_atom, n_point_per_grid, d_point_x, d_point_y, d_point_z, d_point_w, d_atoms, d_i_atom_for_point,
                                                            d_gradient_cache_x, d_gradient_cache_y, d_gradient_cache_z, d_interatomic_quantities, switch_function_threshold, image_cutoff_radius, periodic_data);
    }

    __global__ void weight_gradient_center_atom_kernel(const int n_point, const int n_atom, const int n_point_per_grid, const int* d_i_atom_for_point,
                                                       double* d_gradient_cache_x, double* d_gradient_cache_y, double* d_gradient_cache_z)
    {
        const int i_point = threadIdx.x + blockIdx.x * blockDim.x;
        if (i_point >= n_point) return;
        const int i_center_atom = d_i_atom_for_point[i_point];

        double gx = 0.0;
        double gy = 0.0;
        double gz = 0.0;
        for (int i_atom = 0; i_atom < n_atom; i_atom++) {
            gx -= d_gradient_cache_x[i_point + i_atom * n_point_per_grid];
            gy -= d_gradient_cache_y[i_point + i_atom * n_point_per_grid];
            gz -= d_gradient_cache_z[i_point + i_atom * n_point_per_grid];
        }

        d_gradient_cache_x[i_point + i_center_atom * n_point_per_grid] = gx;
        d_gradient_cache_y[i_point + i_center_atom * n_point_per_grid] = gy;
        d_gradient_cache_z[i_point + i_center_atom * n_point_per_grid] = gz;
    }

    void weight_gradient_center_atom(const int n_grid, const int n_block, const int n_point, const int n_atom, const int n_point_per_grid, const int* d_i_atom_for_point,
                                     double* d_gradient_cache_x, double* d_gradient_cache_y, double* d_gradient_cache_z)
    {
        weight_gradient_center_atom_kernel<<<n_grid, n_block>>>(n_point, n_atom, n_point_per_grid, d_i_atom_for_point, d_gradient_cache_x, d_gradient_cache_y, d_gradient_cache_z);
    }

    // In this kernel, each block maps to (i_atom * 3 + i_xyz), and the whole block goes through all points and perform an internal summation.
    __global__ void weight_gradient_sum_over_point_kernel(const int n_point, double* final_gradient,
                                                          const double* d_gradient_cache_x, const double* d_gradient_cache_y, const double* d_gradient_cache_z)
    {
        __shared__ double shared_partial_output[weight_gradient_sum_over_point_dimension];
        const int i_atom = blockIdx.x / 3;
        const int xyz_dimension_select = blockIdx.x % 3;
        const double* d_gradient_cache = (xyz_dimension_select == 0) ? d_gradient_cache_x : ( (xyz_dimension_select == 1) ? d_gradient_cache_y : d_gradient_cache_z );

        const double* cur = d_gradient_cache + n_point * i_atom + threadIdx.x;
        const double* end = d_gradient_cache + n_point * (i_atom+1);

        double sum = 0.0;
        while (cur < end) {
            sum += *cur;
            cur += weight_gradient_sum_over_point_dimension;
        }
        shared_partial_output[threadIdx.x] = sum;

        for (int span = weight_gradient_sum_over_point_dimension / 2; span > 1; span /= 2) {
            __syncthreads();
            if (threadIdx.x < span)
                shared_partial_output[threadIdx.x] += shared_partial_output[threadIdx.x + span];
        }
        if (threadIdx.x == 0) {
            sum = shared_partial_output[0] + shared_partial_output[1];
            final_gradient[blockIdx.x] = sum;
        }
    }

    void weight_gradient_sum_over_point(const int n_grid, const int n_block, const int n_point, double* final_gradient,
                                        const double* d_gradient_cache_x, const double* d_gradient_cache_y, const double* d_gradient_cache_z)
    {
        weight_gradient_sum_over_point_kernel<<<n_grid, n_block>>>(n_point, final_gradient, d_gradient_cache_x, d_gradient_cache_y, d_gradient_cache_z);
    }

}
