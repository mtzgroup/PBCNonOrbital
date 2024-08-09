#ifndef PERIODIC_NUCLEAR_REPULSION_H_
#define PERIODIC_NUCLEAR_REPULSION_H_

#include "periodic_types.h"

namespace PeriodicBox
{
  double PeriodicNuclearRepulsion(const double full_range_scale,
                                  const double long_range_scale,
                                  const PeriodicParameter periodic_parameter,
                                  const int n_charge,
                                  const double* charge_xyzc);

  void PeriodicNuclearRepulsionGradient(const double full_range_scale,
                                        const double long_range_scale,
                                        const PeriodicParameter periodic_parameter,
                                        const int n_charge,
                                        const double* charge_xyzc,
                                        double* gradient);
}

#endif