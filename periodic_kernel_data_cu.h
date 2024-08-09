#ifndef PERIODIC_KERNEL_DATA_CU_H_
#define PERIODIC_KERNEL_DATA_CU_H_

#include "periodic_types.h"

#include <normc.h> // From intbox, for PI

#include <math.h>

#define FABS_MAX(X, Y) ((fabs(X) > fabs(Y)) ? fabs(X) : fabs(Y))
#define CAST(x) (static_cast<FloatType>(x))

namespace PeriodicBox
{
  template <typename FloatType>
  class PeriodicKernelDataGeneral
  {
  protected:
    const FloatType omega_;
    const FloatType one_over_omega_2_;
    const double scale_;
    const FloatType lattice_real_[9];
    const FloatType lattice_reciprocal_[9];

  public:
    PeriodicKernelDataGeneral(const double omega, const double scale, const LatticeVector lattice);

#if defined(__NVCC__)
    __device__ FloatType get_omega() const { return omega_; }
    __device__ FloatType get_one_over_omega_2() const { return one_over_omega_2_; }
    __device__ double get_scale() const { return scale_; }

    __device__
    void get_absolute_coord_real(FloatType absolute[3], const int i_x, const int i_y, const int i_z) const
    {
      absolute[0] = lattice_real_[0] * i_x + lattice_real_[3] * i_y + lattice_real_[6] * i_z;
      absolute[1] = lattice_real_[1] * i_x + lattice_real_[4] * i_y + lattice_real_[7] * i_z;
      absolute[2] = lattice_real_[2] * i_x + lattice_real_[5] * i_y + lattice_real_[8] * i_z;
    }

    __device__
    void get_absolute_coord_reciprocal(FloatType absolute[3], const int i_x, const int i_y, const int i_z) const
    {
      absolute[0] = lattice_reciprocal_[0] * i_x + lattice_reciprocal_[3] * i_y + lattice_reciprocal_[6] * i_z;
      absolute[1] = lattice_reciprocal_[1] * i_x + lattice_reciprocal_[4] * i_y + lattice_reciprocal_[7] * i_z;
      absolute[2] = lattice_reciprocal_[2] * i_x + lattice_reciprocal_[5] * i_y + lattice_reciprocal_[8] * i_z;
    }

    __device__ void get_min_norm_lattice_image_real(FloatType v[3]) const
    {
      const FloatType lattice_index_with_fraction[3] {
        CAST(0.5 / PI) * (lattice_reciprocal_[0] * v[0] + lattice_reciprocal_[1] * v[1] + lattice_reciprocal_[2] * v[2]),
        CAST(0.5 / PI) * (lattice_reciprocal_[3] * v[0] + lattice_reciprocal_[4] * v[1] + lattice_reciprocal_[5] * v[2]),
        CAST(0.5 / PI) * (lattice_reciprocal_[6] * v[0] + lattice_reciprocal_[7] * v[1] + lattice_reciprocal_[8] * v[2]),
      };
      const FloatType lattice_index[3] {
        round(lattice_index_with_fraction[0]),
        round(lattice_index_with_fraction[1]),
        round(lattice_index_with_fraction[2]),
      };

      v[0] -= lattice_real_[0] * lattice_index[0] + lattice_real_[3] * lattice_index[1] + lattice_real_[6] * lattice_index[2];
      v[1] -= lattice_real_[1] * lattice_index[0] + lattice_real_[4] * lattice_index[1] + lattice_real_[7] * lattice_index[2];
      v[2] -= lattice_real_[2] * lattice_index[0] + lattice_real_[5] * lattice_index[1] + lattice_real_[8] * lattice_index[2];
    }

    protected:
    // |P - origin| < radius
    __device__
    void get_cube_bound_general(int positive_bound[3], int negative_bound[3], const FloatType origin_offset_absolute[3], const FloatType radius, const FloatType P_inverse[9]) const
    {
      const FloatType positive_direction_standard_basis[3] {
        radius - origin_offset_absolute[0],
        radius - origin_offset_absolute[1],
        radius - origin_offset_absolute[2],
      };
      const FloatType negative_direction_standard_basis[9] {
        -radius - origin_offset_absolute[0],
        -radius - origin_offset_absolute[1],
        -radius - origin_offset_absolute[2],
      };

      const FloatType positive_direction_lattice_basis[9] {
        P_inverse[0] * positive_direction_standard_basis[0],
        P_inverse[1] * positive_direction_standard_basis[0],
        P_inverse[2] * positive_direction_standard_basis[0],
        P_inverse[3] * positive_direction_standard_basis[1],
        P_inverse[4] * positive_direction_standard_basis[1],
        P_inverse[5] * positive_direction_standard_basis[1],
        P_inverse[6] * positive_direction_standard_basis[2],
        P_inverse[7] * positive_direction_standard_basis[2],
        P_inverse[8] * positive_direction_standard_basis[2],
      };
      const FloatType negative_direction_lattice_basis[9] {
        P_inverse[0] * negative_direction_standard_basis[0],
        P_inverse[1] * negative_direction_standard_basis[0],
        P_inverse[2] * negative_direction_standard_basis[0],
        P_inverse[3] * negative_direction_standard_basis[1],
        P_inverse[4] * negative_direction_standard_basis[1],
        P_inverse[5] * negative_direction_standard_basis[1],
        P_inverse[6] * negative_direction_standard_basis[2],
        P_inverse[7] * negative_direction_standard_basis[2],
        P_inverse[8] * negative_direction_standard_basis[2],
      };

      for (int i = 0; i < 3; i++)
      {
        positive_bound[i] = (int) ceil(FABS_MAX(FABS_MAX(positive_direction_lattice_basis[i + 0],
                                                         positive_direction_lattice_basis[i + 3]),
                                                         positive_direction_lattice_basis[i + 6]));
        negative_bound[i] = (int) ceil(FABS_MAX(FABS_MAX(negative_direction_lattice_basis[i + 0],
                                                         negative_direction_lattice_basis[i + 3]),
                                                         negative_direction_lattice_basis[i + 6]));
      }
    }
#endif
  };

  template <typename FloatType>
  class PeriodicKernelDataReal : public PeriodicKernelDataGeneral<FloatType>
  {
  protected:
    const FloatType log_threshold_P1_;

  public:
    PeriodicKernelDataReal(const double omega, const double real_scale, const double threshold_P1, const LatticeVector lattice);

#if defined(__NVCC__)
    __device__
    void get_cube_bound_real(int positive_bound[3], int negative_bound[3], const FloatType origin_offset_absolute[3], const FloatType radius) const
    {
      const FloatType one_over_two_pi = 0.5 / PI;
      const FloatType P_inverse[9] {
        this->lattice_reciprocal_[0] * one_over_two_pi,
        this->lattice_reciprocal_[3] * one_over_two_pi,
        this->lattice_reciprocal_[6] * one_over_two_pi,
        this->lattice_reciprocal_[1] * one_over_two_pi,
        this->lattice_reciprocal_[4] * one_over_two_pi,
        this->lattice_reciprocal_[7] * one_over_two_pi,
        this->lattice_reciprocal_[2] * one_over_two_pi,
        this->lattice_reciprocal_[5] * one_over_two_pi,
        this->lattice_reciprocal_[8] * one_over_two_pi,
      };

      get_cube_bound_general(positive_bound, negative_bound, origin_offset_absolute, radius, P_inverse);
    }

    __device__
    FloatType get_1e_coulomb_cutoff_distance_real(const FloatType exponent_p, const FloatType prefactor) const
    {
      const FloatType ln_sum = log(fabs(prefactor)) - log_threshold_P1_ + log(CAST(1.0) - CAST(1.0) / sqrt(this->one_over_omega_2_ * exponent_p + CAST(1.0)));
      return (ln_sum <= CAST(0.0)) ? CAST(0.0) : sqrt(ln_sum * (CAST(1.0) / exponent_p + this->one_over_omega_2_));
    }

    __device__
    FloatType get_2e_coulomb_cutoff_distance_real(const FloatType exponent_p, const FloatType exponent_q, const FloatType prefactor) const
    {
      const FloatType exponent_inverse_sum = CAST(1.0) / exponent_p + CAST(1.0) / exponent_q;
      const FloatType ln_sum = log(fabs(prefactor)) - log_threshold_P1_ + log(CAST(1.0) - CAST(1.0) / sqrt(this->one_over_omega_2_ / exponent_inverse_sum + CAST(1.0)));
      return (ln_sum <= CAST(0.0)) ? CAST(0.0) : sqrt(ln_sum * (exponent_inverse_sum + this->one_over_omega_2_));
    }

    __device__
    FloatType get_orbital_cutoff_distance_real(const FloatType exponent_p, const FloatType prefactor) const
    {
      const FloatType ln_sum = log(fabs(prefactor)) - log_threshold_P1_;
      return (ln_sum <= CAST(0.0)) ? CAST(0.0) : sqrt(ln_sum / exponent_p);
    }

    __device__
    void move_to_same_image(const FloatType source[3], FloatType to_move[3]) const
    {
      const FloatType source_index_with_fraction[3] {
        CAST(0.5 / PI) * (this->lattice_reciprocal_[0] * source[0] + this->lattice_reciprocal_[1] * source[1] + this->lattice_reciprocal_[2] * source[2]),
        CAST(0.5 / PI) * (this->lattice_reciprocal_[3] * source[0] + this->lattice_reciprocal_[4] * source[1] + this->lattice_reciprocal_[5] * source[2]),
        CAST(0.5 / PI) * (this->lattice_reciprocal_[6] * source[0] + this->lattice_reciprocal_[7] * source[1] + this->lattice_reciprocal_[8] * source[2]),
      };
      const FloatType source_index[3] {
        round(source_index_with_fraction[0]),
        round(source_index_with_fraction[1]),
        round(source_index_with_fraction[2]),
      };

      const FloatType to_move_index_with_fraction[3] {
        CAST(0.5 / PI) * (this->lattice_reciprocal_[0] * to_move[0] + this->lattice_reciprocal_[1] * to_move[1] + this->lattice_reciprocal_[2] * to_move[2]),
        CAST(0.5 / PI) * (this->lattice_reciprocal_[3] * to_move[0] + this->lattice_reciprocal_[4] * to_move[1] + this->lattice_reciprocal_[5] * to_move[2]),
        CAST(0.5 / PI) * (this->lattice_reciprocal_[6] * to_move[0] + this->lattice_reciprocal_[7] * to_move[1] + this->lattice_reciprocal_[8] * to_move[2]),
      };
      const FloatType to_move_index[3] {
        round(to_move_index_with_fraction[0]),
        round(to_move_index_with_fraction[1]),
        round(to_move_index_with_fraction[2]),
      };

      const FloatType delta_index[3] {
        source_index[0] - to_move_index[0],
        source_index[1] - to_move_index[1],
        source_index[2] - to_move_index[2],
      };
      to_move[0] += this->lattice_real_[0] * delta_index[0] + this->lattice_real_[3] * delta_index[1] + this->lattice_real_[6] * delta_index[2];
      to_move[1] += this->lattice_real_[1] * delta_index[0] + this->lattice_real_[4] * delta_index[1] + this->lattice_real_[7] * delta_index[2];
      to_move[2] += this->lattice_real_[2] * delta_index[0] + this->lattice_real_[5] * delta_index[1] + this->lattice_real_[8] * delta_index[2];
    }
#endif
  };

  template <typename FloatType>
  class PeriodicKernelDataReciprocal : public PeriodicKernelDataGeneral<FloatType>
  {
  protected:
    const FloatType log_threshold_G1_;
    const FloatType one_over_V_real_;
    const FloatType log_V_real_K_min_square_;

  public:
    PeriodicKernelDataReciprocal(const double omega, const double reciprocal_scale, const double threshold_G1,
                                  const double V_real, const double K_min_norm, const LatticeVector lattice);

#if defined(__NVCC__)
    __device__ FloatType get_one_over_V_real() const { return one_over_V_real_; }

    __device__
    void get_cube_bound_reciprocal(int positive_bound[3], int negative_bound[3], const FloatType radius) const
    {
      const FloatType one_over_two_pi = 0.5 / PI;
      const FloatType P_inverse[9] {
        this->lattice_real_[0] * one_over_two_pi,
        this->lattice_real_[3] * one_over_two_pi,
        this->lattice_real_[6] * one_over_two_pi,
        this->lattice_real_[1] * one_over_two_pi,
        this->lattice_real_[4] * one_over_two_pi,
        this->lattice_real_[7] * one_over_two_pi,
        this->lattice_real_[2] * one_over_two_pi,
        this->lattice_real_[5] * one_over_two_pi,
        this->lattice_real_[8] * one_over_two_pi,
      };

      const FloatType origin[3] { 0.0, 0.0, 0.0 };
      get_cube_bound_general(positive_bound, negative_bound, origin, radius, P_inverse);
    }

    __device__
    FloatType get_1e_coulomb_cutoff_distance_reciprocal(const FloatType exponent_p, const FloatType prefactor) const
    {
      const FloatType ln_sum = log(fabs(prefactor)) - log_threshold_G1_ + CAST(2.4102420093340458) // log(2.0 * pow(PI, 1.5))
                        - log_V_real_K_min_square_ - CAST(0.5) * log(exponent_p);
      return (ln_sum <= CAST(0.0)) ? CAST(0.0) : sqrt(ln_sum * CAST(4.0) / (CAST(1.0) / exponent_p + this->one_over_omega_2_));
    }

    __device__
    FloatType get_2e_coulomb_cutoff_distance_reciprocal(const FloatType exponent_p, const FloatType exponent_q, const FloatType prefactor) const
    {
      const FloatType exponent_inverse_sum = CAST(1.0) / exponent_p + CAST(1.0) / exponent_q;
      const FloatType ln_sum = log(fabs(prefactor)) - log_threshold_G1_ + CAST(2.4102420093340458) // log(2.0 * pow(PI, 1.5))
                            - log_V_real_K_min_square_ + CAST(0.5) * log(exponent_inverse_sum);
      return (ln_sum <= CAST(0.0)) ? CAST(0.0) : sqrt(ln_sum * CAST(4.0) / (exponent_inverse_sum + this->one_over_omega_2_));
    }
#endif
  };
}

#endif