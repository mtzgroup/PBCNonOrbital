#include "periodic_cutoff.h"

#include <math.h>

namespace PeriodicBox
{
  double get_1e_orbital_cutoff_distance_real(const double exponent1, const double exponent2,
                                            const double threshold, const double prefactor)
  {
    const double C_P2 = fabs(prefactor);
    const double ln_sum = - log(threshold) + log(C_P2) + log(2.0 * M_PI) - log(exponent1 + exponent2);
    if (ln_sum <= 0.0)
      return 0.0;
    return sqrt(ln_sum * (1.0 / exponent1 + 1.0 / exponent2));
  }

  double get_2e_orbital_cutoff_distance_real(const double exponent1, const double exponent2,
                                            const double threshold, const double prefactor,
                                            const double min_exponent)
  {
    const double C_P2 = fabs(prefactor);
    const double ln_sum = - log(threshold) + log(C_P2) + 0.5 * log(2.0) + 1.25 * log(M_PI)
                          - log(exponent1 + exponent2) - 0.5 * log(exponent1 + exponent2 + min_exponent * 2.0);
    if (ln_sum <= 0.0)
      return 0.0;
    return sqrt(ln_sum * (1.0 / exponent1 + 1.0 / exponent2));
  }

  double get_nuclear_cutoff_distance_real(const double threshold, const double prefactor, const double omega)
  {
    // erfc(x) <= 2/sqrt(pi) * exp(-x^2) / (x + sqrt(x^2 + 4/pi)) < 2/sqrt(pi) * exp(-x^2) / x
    // C * erfc(wr) / r < C * 2/sqrt(pi) * exp(-w^2 * r^2) / r^2 < threshold
    const double C = fabs(prefactor);
    const double log_factors = log(2.0) - 0.5 * log(M_PI) + log(C) - log(threshold);
    // w^2 * r^2 > log_factors + 2 * log(r)
    // Fix point iteration
    double r = sqrt((log_factors <= 1.0) ? 1.0 : log_factors) / omega;
    r = sqrt((log_factors + 2 * log(r) <= 1.0) ? 1.0 : (log_factors + 2 * log(r))) / omega;
    r = sqrt((log_factors + 2 * log(r) <= 1.0) ? 1.0 : (log_factors + 2 * log(r))) / omega;
    r = sqrt((log_factors + 2 * log(r) <= 1.0) ? 1.0 : (log_factors + 2 * log(r))) / omega;
    return r;
  }

  double get_nuclear_cutoff_distance_reciprocal(const double threshold, const double prefactor, const double omega)
  {
    // C * exp(-k^2 / (4 * w^2)) / k^2 < threshold
    const double C = fabs(prefactor);
    const double log_factors = log(C) - log(threshold);
    // k^2 / (4 * w^2) > log_factors + 2 * log(k)
    // Fix point iteration
    double k = sqrt((log_factors <= 1.0) ? 1.0 : log_factors) * 2.0 * omega;
    k = sqrt((log_factors + 2 * log(k) <= 1.0) ? 1.0 : (log_factors + 2 * log(k))) * 2.0 * omega;
    k = sqrt((log_factors + 2 * log(k) <= 1.0) ? 1.0 : (log_factors + 2 * log(k))) * 2.0 * omega;
    k = sqrt((log_factors + 2 * log(k) <= 1.0) ? 1.0 : (log_factors + 2 * log(k))) * 2.0 * omega;
    return k;
  }
}