#include "periodic_kernel_data_cu.h"

// Henry 20240227: If you wonder, why don't we put these simple constructors into the header file?
//                 It has nothing to do with the TeraChem rules or Todd's (and my) opinion against header implementation.
//                 It's because the header file will be included by .cu files.
//                 When compiling the .cu files on an old OS like CentOS 7, we have to use -std=c++98 for the host compiler.
//                 And that's source of evil. The list initializer, as showing up below, is a c++11 feature.
//                 Also, if you wonder why we're using NULL instead of nullptr here and there, it's for the same reason.

namespace PeriodicBox
{
  template <typename FType>
  PeriodicKernelDataGeneral<FType>::PeriodicKernelDataGeneral(const double omega, const double scale, const LatticeVector lattice)
  : omega_(omega),
    one_over_omega_2_(1.0 / omega / omega),
    scale_(scale),
    lattice_real_{ lattice.real_space[0], lattice.real_space[1], lattice.real_space[2],
                    lattice.real_space[3], lattice.real_space[4], lattice.real_space[5],
                    lattice.real_space[6], lattice.real_space[7], lattice.real_space[8], },
    lattice_reciprocal_{ lattice.reciprocal_space[0], lattice.reciprocal_space[1], lattice.reciprocal_space[2],
                          lattice.reciprocal_space[3], lattice.reciprocal_space[4], lattice.reciprocal_space[5],
                          lattice.reciprocal_space[6], lattice.reciprocal_space[7], lattice.reciprocal_space[8], }
  {}
  template PeriodicKernelDataGeneral<float >::PeriodicKernelDataGeneral(const double omega, const double scale, const LatticeVector lattice);
  template PeriodicKernelDataGeneral<double>::PeriodicKernelDataGeneral(const double omega, const double scale, const LatticeVector lattice);

  template <typename FType>
  PeriodicKernelDataReal<FType>::PeriodicKernelDataReal(const double omega, const double real_scale, const double threshold_P1, const LatticeVector lattice)
  : PeriodicKernelDataGeneral<FType>(omega, real_scale, lattice),
    log_threshold_P1_(log(threshold_P1))
  {}
  template PeriodicKernelDataReal<float >::PeriodicKernelDataReal(const double omega, const double real_scale, const double threshold_P1, const LatticeVector lattice);
  template PeriodicKernelDataReal<double>::PeriodicKernelDataReal(const double omega, const double real_scale, const double threshold_P1, const LatticeVector lattice);

  template <typename FType>
  PeriodicKernelDataReciprocal<FType>::PeriodicKernelDataReciprocal(const double omega, const double reciprocal_scale, const double threshold_G1,
                                                                    const double V_real, const double K_min_norm, const LatticeVector lattice)
  : PeriodicKernelDataGeneral<FType>(omega, reciprocal_scale, lattice),
    log_threshold_G1_(log(threshold_G1)),
    one_over_V_real_(1.0 / V_real),
    log_V_real_K_min_square_(log(V_real) + 2.0 * log(K_min_norm))
  {}
  template PeriodicKernelDataReciprocal<float >::PeriodicKernelDataReciprocal(const double omega, const double reciprocal_scale, const double threshold_G1, const double V_real, const double K_min_norm, const LatticeVector lattice);
  template PeriodicKernelDataReciprocal<double>::PeriodicKernelDataReciprocal(const double omega, const double reciprocal_scale, const double threshold_G1, const double V_real, const double K_min_norm, const LatticeVector lattice);
}
