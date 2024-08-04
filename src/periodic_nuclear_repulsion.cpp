#include "periodic_nuclear_repulsion.h"

#include "periodic_lattice_internal.h"
#include "periodic_cutoff.h"

#include <math.h>
#include <stdio.h>

namespace PeriodicBox
{
  double PeriodicNuclearRepulsion(const double full_range_scale,
                                  const double long_range_scale,
                                  const PeriodicParameter periodic_parameter,
                                  const int n_charge,
                                  const double* charge_xyzc)
  {
    const double omega = periodic_parameter.omega;

    const double nuclear_cutoff_real = get_nuclear_cutoff_distance_real(periodic_parameter.thresholds.periodic_charge_cutoff_real, 1.0, omega);
    const double nuclear_cutoff_reciprocal = get_nuclear_cutoff_distance_reciprocal(periodic_parameter.thresholds.periodic_charge_cutoff_reciprocal, 4.0 * M_PI / periodic_parameter.lattice.V_real, omega);
    // printf("nuclear_cutoff_real = %.5f, nuclear_cutoff_reciprocal = %.5f\n", nuclear_cutoff_real, nuclear_cutoff_reciprocal);

    const std::vector<double3> reciprocal_vectors = PeriodicBox::compute_lattice_vectors_within_sphere_reciprocal(periodic_parameter.lattice.unit_cell, nuclear_cutoff_reciprocal);

    double E_real = 0.0;
    double E_reciprocal = 0.0;
    for (int i_charge = 0; i_charge < n_charge; i_charge++) {
      const double q_i = charge_xyzc[i_charge * 4 + 3];
      const double P_i[3] { charge_xyzc[i_charge * 4 + 0], charge_xyzc[i_charge * 4 + 1], charge_xyzc[i_charge * 4 + 2] };
      for (int j_charge = i_charge; j_charge < n_charge; j_charge++) {
        const double q_j = charge_xyzc[j_charge * 4 + 3];
        const double P_j[3] { charge_xyzc[j_charge * 4 + 0], charge_xyzc[j_charge * 4 + 1], charge_xyzc[j_charge * 4 + 2] };

        const double rij[3] = { P_i[0] - P_j[0], P_i[1] - P_j[1], P_i[2] - P_j[2] };

        // Real
        const std::vector<double3> real_vectors = PeriodicBox::compute_lattice_vectors_within_sphere_real(periodic_parameter.lattice.unit_cell, rij, nuclear_cutoff_real, i_charge == j_charge);

        double E_real_ij = 0.0;
        for (int i_R = 0; i_R < real_vectors.size(); i_R++) {
          const double3& image_R = real_vectors[i_R];
          const double r_R[3] = { rij[0] + image_R.x, rij[1] + image_R.y, rij[2] + image_R.z };
          const double r_R_norm = sqrt(r_R[0] * r_R[0] + r_R[1] * r_R[1] + r_R[2] * r_R[2]);

          E_real_ij += erfc(omega * r_R_norm) / r_R_norm;
        }
        if (i_charge == j_charge) {
          E_real_ij /= 2.0;
          E_real_ij -= omega / sqrt(M_PI); // Self
        }
        E_real += E_real_ij * q_i * q_j;

        // Reciporcal
        double E_reciprocal_ij = 0.0;
        for (int i_G = 0; i_G < reciprocal_vectors.size(); i_G++) {
          const double3& image_G = reciprocal_vectors[i_G];
          const double G_2 = image_G.x * image_G.x + image_G.y * image_G.y + image_G.z * image_G.z;
          const double G_dot_delta_R = rij[0] * image_G.x + rij[1] * image_G.y + rij[2] * image_G.z;

          E_reciprocal_ij += exp( -0.25 / (omega * omega) * G_2 ) / G_2 * cos(G_dot_delta_R);
        }
        E_reciprocal_ij -= 0.25 / (omega * omega); // Self, notice that for nucleus only, we don't have $\sum_i q_i = 0$
        if (i_charge == j_charge) {
          E_reciprocal_ij /= 2.0;
        }
        E_reciprocal += E_reciprocal_ij * q_i * q_j;

        // printf("i = %d, j = %d, q_i = %.5f, q_j = %.5f, E_real_ij = %.5f, E_reciprocal_ij = %.5f, E_total_ij = %.5f\n", i_charge, j_charge, q_i, q_j,
        //   E_real_ij, E_reciprocal_ij * 4.0 * M_PI / periodic_parameter.lattice.V_real, E_real_ij + E_reciprocal_ij * 4.0 * M_PI / periodic_parameter.lattice.V_real);
      }
    }

    E_real *= full_range_scale;
    E_reciprocal *= (full_range_scale + long_range_scale) * 4.0 * M_PI / periodic_parameter.lattice.V_real;

    return E_real + E_reciprocal;
  }
}