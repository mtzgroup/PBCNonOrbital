#ifndef PERIODIC_CUTOFF_H_
#define PERIODIC_CUTOFF_H_

#include "periodic_types.h"

namespace PeriodicBox
{
  double get_1e_orbital_cutoff_distance_real(const double exponent1, const double exponent2,
                                            const double threshold, const double prefactor);
  double get_2e_orbital_cutoff_distance_real(const double exponent1, const double exponent2,
                                            const double threshold, const double prefactor,
                                            const double min_exponent);
  double get_nuclear_cutoff_distance_real(const double threshold, const double prefactor, const double omega);
  double get_nuclear_cutoff_distance_reciprocal(const double threshold, const double prefactor, const double omega);
}

#endif
