#include "../src/periodic_lattice_export.h"
#include "../src/helper.h"
#include "../src/periodic_grid.h"
#include "../src/periodic_becke/periodic_becke.h"

#include <string.h>

double becke_energy(const int n_point, const PeriodicBox::GridPoint* points, const double* epsilon_xc)
{
    double E_xc = 0.0;
    for (int i_point = 0; i_point < n_point; i_point++) {
        E_xc += points[i_point].w_total * epsilon_xc[i_point];
    }
    return E_xc;
}

int main()
{
    PeriodicBox::LatticeInfo lattice;
    lattice.dimension = 3;
    lattice.a = 2.0 / BohrToAng;
    lattice.b = 10.0 / BohrToAng;
    lattice.c = 10.0 / BohrToAng;
    lattice.alpha = 90.0 * PI / 180.0;
    lattice.beta  = 90.0 * PI / 180.0;
    lattice.gamma = 90.0 * PI / 180.0;
    PeriodicBox::LatticeInfoMethods::calculate_lattice_vector(&lattice);
    PeriodicBox::LatticeInfoMethods::print(lattice);

    PeriodicBox::PeriodicParameter periodic_parameter;
    periodic_parameter.min_primitive_exponent = NAN;
    periodic_parameter.max_primitive_exponent = NAN;
    periodic_parameter.lattice = lattice;
    periodic_parameter.thresholds.periodic_orbtial_cutoff = NAN;
    periodic_parameter.thresholds.periodic_charge_cutoff_real = 1e-14;
    periodic_parameter.thresholds.periodic_charge_cutoff_reciprocal = 1e-14;
    periodic_parameter.omega = 0.2;
    printf("\nomega = %.5e\n", periodic_parameter.omega);

    const int n_atom = 2;
    const double atom_xyz[n_atom * 3] {
        0.0 / BohrToAng, 0.0 / BohrToAng, 0.0 / BohrToAng,
        1.0 / BohrToAng, 0.0 / BohrToAng, 0.0 / BohrToAng,
    };
    const double atom_radius[n_atom] {
        0.35 / BohrToAng,
        0.35 / BohrToAng,
    };

    const int n_grid_point = 4;
    PeriodicBox::GridPoint grid_points[n_grid_point];
    grid_points[0].x = 1.0 / BohrToAng; grid_points[0].y = 0.0 / BohrToAng; grid_points[0].z = 0.0 / BohrToAng;
    grid_points[1].x = 0.4 / BohrToAng; grid_points[1].y = 0.0 / BohrToAng; grid_points[1].z = 0.0 / BohrToAng;
    grid_points[2].x = 1.0 / BohrToAng; grid_points[2].y = 0.0 / BohrToAng; grid_points[2].z = 0.0 / BohrToAng;
    grid_points[3].x = 1.5 / BohrToAng; grid_points[3].y = 0.0 / BohrToAng; grid_points[3].z = 0.0 / BohrToAng;
    grid_points[0].w_fixed = 1.0; grid_points[0].w_total = NAN; grid_points[0].i_atom = 0;
    grid_points[1].w_fixed = 1.0; grid_points[1].w_total = NAN; grid_points[1].i_atom = 0;
    grid_points[2].w_fixed = 1.0; grid_points[2].w_total = NAN; grid_points[2].i_atom = 1;
    grid_points[3].w_fixed = 1.0; grid_points[3].w_total = NAN; grid_points[3].i_atom = 1;

    const double epsilon_xc[n_grid_point] {
        0.0,
        1.0,
        0.0,
        0.0,
    };

    double* analytical_gradient = new double[n_atom * 3];
    double* numerical_gradient  = new double[n_atom * 3];
    memset(analytical_gradient, 0, n_atom * 3 * sizeof(double));
    memset(numerical_gradient,  0, n_atom * 3 * sizeof(double));

    PeriodicBox::BeckeWeightGradient(n_grid_point, grid_points, n_atom, atom_xyz, atom_radius, periodic_parameter.lattice.unit_cell, epsilon_xc, analytical_gradient);

    for (int i_atom = 0; i_atom < n_atom; i_atom++)
        printf("Analytical gradient %2d  %15.10f  %15.10f  %15.10f\n",
            i_atom, analytical_gradient[i_atom * 3 + 0], analytical_gradient[i_atom * 3 + 1], analytical_gradient[i_atom * 3 + 2]);

    double* atom_xyz_copy = new double[n_atom * 3];
    memcpy(atom_xyz_copy, atom_xyz, n_atom * 3 * sizeof(double));

    const double dx = 1e-4;

    for (int i_atom = 0; i_atom < n_atom; i_atom++)
        for (int i_xyz = 0; i_xyz < 3; i_xyz++) {
            atom_xyz_copy[i_atom * 3 + i_xyz] = atom_xyz[i_atom * 3 + i_xyz] + dx;
            for (int i_point = 0; i_point < n_grid_point; i_point++) grid_points[i_point].w_total = NAN;
            BeckeWeights(n_grid_point, grid_points, n_atom, atom_xyz_copy, atom_radius, periodic_parameter.lattice.unit_cell);
            const double E_plus = becke_energy(n_grid_point, grid_points, epsilon_xc);

            atom_xyz_copy[i_atom * 3 + i_xyz] = atom_xyz[i_atom * 3 + i_xyz] - dx;
            for (int i_point = 0; i_point < n_grid_point; i_point++) grid_points[i_point].w_total = NAN;
            BeckeWeights(n_grid_point, grid_points, n_atom, atom_xyz_copy, atom_radius, periodic_parameter.lattice.unit_cell);
            const double E_minus = becke_energy(n_grid_point, grid_points, epsilon_xc);

            numerical_gradient[i_atom * 3 + i_xyz] = (E_plus - E_minus) / dx / 2.0;
            atom_xyz_copy[i_atom * 3 + i_xyz] = atom_xyz[i_atom * 3 + i_xyz];
        }

    for (int i_atom = 0; i_atom < n_atom; i_atom++)
        printf("Numerical gradient %2d  %15.10f  %15.10f  %15.10f\n",
            i_atom, numerical_gradient[i_atom * 3 + 0], numerical_gradient[i_atom * 3 + 1], numerical_gradient[i_atom * 3 + 2]);

    double max_abs_diff = 0.0;
    for (int i_atom = 0; i_atom < n_atom; i_atom++)
        for (int i_xyz = 0; i_xyz < 3; i_xyz++) {
            const double abs_diff = fabs(analytical_gradient[i_atom * 3 + i_xyz] - numerical_gradient[i_atom * 3 + i_xyz]);
            if (abs_diff > max_abs_diff) max_abs_diff = abs_diff;
        }
    printf("\nMax abs diff between analytical and numerical: %.5e\n", max_abs_diff);


    delete[] analytical_gradient;
    delete[] numerical_gradient;
    delete[] atom_xyz_copy;

    return 0;
}
