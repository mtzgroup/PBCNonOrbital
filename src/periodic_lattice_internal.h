#ifndef PERIODIC_LATTICE_INTERNAL_H_
#define PERIODIC_LATTICE_INTERNAL_H_

#include "periodic_types.h"

#include <vector>
#include <vector_types.h> // For double3

namespace PeriodicBox
{
  // |P + origin| < radius
  std::vector<double3> compute_lattice_vectors_within_sphere_real(const LatticeVector unit_cell,
                                                                  const double origin_offset_absolute[3],
                                                                  const double radius,
                                                                  const bool remove_origin = false);

  // |G| < radius, G != 0
  std::vector<double3> compute_lattice_vectors_within_sphere_reciprocal(const LatticeVector unit_cell,
                                                                        const double radius);

  double3 get_image_in_origin_cell(const LatticeVector unit_cell, const double3 point);
}

#endif
