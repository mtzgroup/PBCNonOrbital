#include "periodic_nuclear_repulsion.h"

#include "periodic_lattice_internal.h"
#include "periodic_cutoff.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

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

  void PeriodicNuclearRepulsionGradient(const double full_range_scale,
                                        const double long_range_scale,
                                        const PeriodicParameter periodic_parameter,
                                        const int n_charge,
                                        const double* charge_xyzc,
                                        double* gradient)
  {
    const double omega = periodic_parameter.omega;

    const double nuclear_cutoff_real = get_nuclear_cutoff_distance_real(periodic_parameter.thresholds.periodic_charge_cutoff_real, 1.0, omega);
    const double nuclear_cutoff_reciprocal = get_nuclear_cutoff_distance_reciprocal(periodic_parameter.thresholds.periodic_charge_cutoff_reciprocal, 4.0 * M_PI / periodic_parameter.lattice.V_real, omega);
    // printf("nuclear_cutoff_real = %.5f, nuclear_cutoff_reciprocal = %.5f\n", nuclear_cutoff_real, nuclear_cutoff_reciprocal);

    const std::vector<double3> reciprocal_vectors = PeriodicBox::compute_lattice_vectors_within_sphere_reciprocal(periodic_parameter.lattice.unit_cell, nuclear_cutoff_reciprocal);

    double* dE_real       = new double[n_charge * 3];
    double* dE_reciprocal = new double[n_charge * 3];
    memset(dE_real,       0, n_charge * 3 * sizeof(double));
    memset(dE_reciprocal, 0, n_charge * 3 * sizeof(double));

    for (int i_charge = 0; i_charge < n_charge; i_charge++) {
      const double q_i = charge_xyzc[i_charge * 4 + 3];
      const double P_i[3] { charge_xyzc[i_charge * 4 + 0], charge_xyzc[i_charge * 4 + 1], charge_xyzc[i_charge * 4 + 2] };
      for (int j_charge = i_charge; j_charge < n_charge; j_charge++) {
        const double q_j = charge_xyzc[j_charge * 4 + 3];
        const double P_j[3] { charge_xyzc[j_charge * 4 + 0], charge_xyzc[j_charge * 4 + 1], charge_xyzc[j_charge * 4 + 2] };

        const double rij[3] = { P_i[0] - P_j[0], P_i[1] - P_j[1], P_i[2] - P_j[2] };

        // Real
        const std::vector<double3> real_vectors = PeriodicBox::compute_lattice_vectors_within_sphere_real(periodic_parameter.lattice.unit_cell, rij, nuclear_cutoff_real, i_charge == j_charge);

        double dE_real_i[3] { 0.0, 0.0, 0.0 };
        double dE_real_j[3] { 0.0, 0.0, 0.0 };
        for (int i_R = 0; i_R < real_vectors.size(); i_R++) {
          const double3& image_R = real_vectors[i_R];
          const double r_R[3] = { rij[0] + image_R.x, rij[1] + image_R.y, rij[2] + image_R.z };
          const double r_2 = r_R[0] * r_R[0] + r_R[1] * r_R[1] + r_R[2] * r_R[2];
          const double r = sqrt(r_2);

          for (int i_xyz = 0; i_xyz < 3; i_xyz++) {
            dE_real_i[i_xyz] += -r_R[i_xyz] / r_2 / r * erfc(omega * r) - 2.0 / sqrt(M_PI) * omega * r_R[i_xyz] / r_2 * exp(-omega * omega * r_2);
            dE_real_j[i_xyz] -= -r_R[i_xyz] / r_2 / r * erfc(omega * r) - 2.0 / sqrt(M_PI) * omega * r_R[i_xyz] / r_2 * exp(-omega * omega * r_2);
          }
        }
        if (i_charge == j_charge) {
          for (int i_xyz = 0; i_xyz < 3; i_xyz++) {
            dE_real_i[i_xyz] /= 2.0;
            dE_real_j[i_xyz] /= 2.0;
          }
        }
        for (int i_xyz = 0; i_xyz < 3; i_xyz++) {
          dE_real[i_charge * 3 + i_xyz] += dE_real_i[i_xyz] * q_i * q_j;
          dE_real[j_charge * 3 + i_xyz] += dE_real_j[i_xyz] * q_i * q_j;
        }

        // Reciporcal
        double dE_reciprocal_i[3] { 0.0, 0.0, 0.0 };
        double dE_reciprocal_j[3] { 0.0, 0.0, 0.0 };
        for (int i_G = 0; i_G < reciprocal_vectors.size(); i_G++) {
          const double3& image_G = reciprocal_vectors[i_G];
          const double G[3] { image_G.x, image_G.y, image_G.z };
          const double G_2 = image_G.x * image_G.x + image_G.y * image_G.y + image_G.z * image_G.z;
          const double G_dot_delta_R = rij[0] * image_G.x + rij[1] * image_G.y + rij[2] * image_G.z;

          for (int i_xyz = 0; i_xyz < 3; i_xyz++) {
            dE_reciprocal_i[i_xyz] += -exp( -0.25 / (omega * omega) * G_2 ) / G_2 * G[i_xyz] * sin(G_dot_delta_R);
            dE_reciprocal_j[i_xyz] -= -exp( -0.25 / (omega * omega) * G_2 ) / G_2 * G[i_xyz] * sin(G_dot_delta_R);
          }
        }
        if (i_charge == j_charge) {
          for (int i_xyz = 0; i_xyz < 3; i_xyz++) {
            dE_reciprocal_i[i_xyz] /= 2.0;
            dE_reciprocal_j[i_xyz] /= 2.0;
          }
        }
        for (int i_xyz = 0; i_xyz < 3; i_xyz++) {
          dE_reciprocal[i_charge * 3 + i_xyz] += dE_reciprocal_i[i_xyz] * q_i * q_j;
          dE_reciprocal[j_charge * 3 + i_xyz] += dE_reciprocal_j[i_xyz] * q_i * q_j;
        }

      }
    }

    for (int i_charge = 0; i_charge < n_charge; i_charge++) {
      for (int i_xyz = 0; i_xyz < 3; i_xyz++) {
        dE_real[i_charge * 3 + i_xyz] *= full_range_scale;
        dE_reciprocal[i_charge * 3 + i_xyz] *= (full_range_scale + long_range_scale) * 4.0 * M_PI / periodic_parameter.lattice.V_real;
      }
    }

    for (int i_charge = 0; i_charge < n_charge; i_charge++) {
      for (int i_xyz = 0; i_xyz < 3; i_xyz++) {
        gradient[i_charge * 3 + i_xyz] +=       dE_real[i_charge * 3 + i_xyz];
        gradient[i_charge * 3 + i_xyz] += dE_reciprocal[i_charge * 3 + i_xyz];
      }
    }

    delete[] dE_real;
    delete[] dE_reciprocal;
  }
}