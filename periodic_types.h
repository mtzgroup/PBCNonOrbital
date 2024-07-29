#ifndef PERIODIC_TYPES_H_
#define PERIODIC_TYPES_H_

#ifdef __cplusplus
extern "C" {
#endif
  namespace PeriodicBox
  {
    struct LatticeVector
    {
      // in au
      double real_space[9]; // R1.x, R1.y, R1.z, R2.x, R2.y, R2.z, R3.x, R3.y, R3.z
      double reciprocal_space[9]; // K1.x, K1.y, K1.z, K2.x, K2.y, K2.z, K3.x, K3.y, K3.z
    };

    struct LatticeInfo
    {
      double a, b, c; // in Bohr
      double alpha, beta, gamma; // in radius
      LatticeVector unit_cell;
      double V_real; // in Bohr^3
      double K_min_norm; // in Bohr
      int dimension;
    };

    // Cutoff for periodic images (summation over \vec{R} \in R^3), has nothing to do with Schwartz bound
    struct PerioidicThreshold
    {
      double periodic_orbtial_cutoff;
      double periodic_charge_cutoff_real;
      double periodic_charge_cutoff_reciprocal;
    };

    struct PeriodicParameter
    {
      double omega;
      double min_primitive_exponent;
      double max_primitive_exponent;

      LatticeInfo lattice;
      PerioidicThreshold thresholds;
    };
  }

#ifdef __cplusplus
}  // End extern C
#endif

#endif