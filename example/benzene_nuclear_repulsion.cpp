#include "../src/periodic_lattice_export.h"
#include "../src/helper.h"
#include "../src/periodic_nuclear_repulsion.h"

#include <string.h>

int main()
{
    PeriodicBox::LatticeInfo lattice;
    lattice.dimension = 3;
    lattice.a = 5.5220 / BohrToAng;
    lattice.b = 5.4396 / BohrToAng;
    lattice.c = 7.6726 / BohrToAng;
    lattice.alpha =  90.00 * PI / 180.0;
    lattice.beta  = 110.55 * PI / 180.0;
    lattice.gamma =  90.00 * PI / 180.0;
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

    const int n_atom = 12 * 2;
    const double xyzc[n_atom * 4] {
         1.7150183861 / BohrToAng, 5.3416872000 / BohrToAng, 6.3725303821 / BohrToAng, 6.0,
         0.9676735997 / BohrToAng, 5.2813076400 / BohrToAng, 5.8222081416 / BohrToAng, 1.0,
         2.7131641240 / BohrToAng, 4.3516800000 / BohrToAng, 6.3186476562 / BohrToAng, 6.0,
         2.6315693528 / BohrToAng, 3.6227736000 / BohrToAng, 5.7467723254 / BohrToAng, 1.0,
         4.5426155165 / BohrToAng, 0.9519300000 / BohrToAng, 0.0524458532 / BohrToAng, 6.0,
         3.8703881035 / BohrToAng, 1.5932588400 / BohrToAng, 0.0898045432 / BohrToAng, 1.0,
        -0.2329259631 / BohrToAng, 2.6218872000 / BohrToAng, 4.4040147962 / BohrToAng, 6.0,
         1.1137098959 / BohrToAng, 0.0979128000 / BohrToAng, 0.8118330701 / BohrToAng, 6.0,
         3.0616542451 / BohrToAng, 2.8177128000 / BohrToAng, 2.7803486560 / BohrToAng, 6.0,
         0.5144188233 / BohrToAng, 2.5615076400 / BohrToAng, 4.9543370366 / BohrToAng, 1.0,
         1.8610546823 / BohrToAng, 0.1582923600 / BohrToAng, 1.3621553105 / BohrToAng, 1.0,
         2.3143094587 / BohrToAng, 2.8780923600 / BohrToAng, 2.2300264156 / BohrToAng, 1.0,
        -1.2310717010 / BohrToAng, 1.6318800000 / BohrToAng, 4.4578975221 / BohrToAng, 6.0,
         0.1155641580 / BohrToAng, 1.0879200000 / BohrToAng, 0.8657157960 / BohrToAng, 6.0,
         4.0597999830 / BohrToAng, 3.8077200000 / BohrToAng, 2.7264659301 / BohrToAng, 6.0,
        -1.1494769298 / BohrToAng, 0.9029736000 / BohrToAng, 5.0297728529 / BohrToAng, 1.0,
         0.1971589292 / BohrToAng, 1.8168264000 / BohrToAng, 1.4375911268 / BohrToAng, 1.0,
         3.9782052118 / BohrToAng, 4.5366264000 / BohrToAng, 2.1545905993 / BohrToAng, 1.0,
        -0.3672513755 / BohrToAng, 3.6717300000 / BohrToAng, 3.5397358729 / BohrToAng, 6.0,
        -1.7138872345 / BohrToAng, 4.4876700000 / BohrToAng, 7.1319175990 / BohrToAng, 6.0,
         3.1959796574 / BohrToAng, 1.7678700000 / BohrToAng, 3.6446275793 / BohrToAng, 6.0,
         0.3049760375 / BohrToAng, 4.3130588400 / BohrToAng, 3.5023771829 / BohrToAng, 1.0,
        -1.0416598215 / BohrToAng, 3.8463411600 / BohrToAng, 7.0945589090 / BohrToAng, 1.0,
         2.5237522445 / BohrToAng, 1.1265411600 / BohrToAng, 3.6819862692 / BohrToAng, 1.0,
    };

    const double full_range_scale = 1.0;
    const double long_range_scale = 0.0;
    const double E = PeriodicBox::PeriodicNuclearRepulsion(full_range_scale, long_range_scale, periodic_parameter, n_atom, xyzc);
    printf("\nE = %.10f\n", E);

    double* analytical_gradient = new double[n_atom * 3];
    double* numerical_gradient  = new double[n_atom * 3];
    memset(analytical_gradient, 0, n_atom * 3 * sizeof(double));
    memset(numerical_gradient,  0, n_atom * 3 * sizeof(double));

    PeriodicBox::PeriodicNuclearRepulsionGradient(full_range_scale, long_range_scale, periodic_parameter, n_atom, xyzc, analytical_gradient);

    for (int i_atom = 0; i_atom < n_atom; i_atom++)
        printf("Analytical gradient %2d  %15.10f  %15.10f  %15.10f\n",
            i_atom, analytical_gradient[i_atom * 3 + 0], analytical_gradient[i_atom * 3 + 1], analytical_gradient[i_atom * 3 + 2]);

    double* xyzc_copy = new double[n_atom * 4];
    memcpy(xyzc_copy, xyzc, n_atom * 4 * sizeof(double));

    const double dx = 1e-4;

    for (int i_atom = 0; i_atom < n_atom; i_atom++)
        for (int i_xyz = 0; i_xyz < 3; i_xyz++) {
            xyzc_copy[i_atom * 4 + i_xyz] = xyzc[i_atom * 4 + i_xyz] + dx;
            const double E_plus  = PeriodicBox::PeriodicNuclearRepulsion(full_range_scale, long_range_scale, periodic_parameter, n_atom, xyzc_copy);

            xyzc_copy[i_atom * 4 + i_xyz] = xyzc[i_atom * 4 + i_xyz] - dx;
            const double E_minus = PeriodicBox::PeriodicNuclearRepulsion(full_range_scale, long_range_scale, periodic_parameter, n_atom, xyzc_copy);

            numerical_gradient[i_atom * 3 + i_xyz] = (E_plus - E_minus) / dx / 2.0;
            xyzc_copy[i_atom * 4 + i_xyz] = xyzc[i_atom * 4 + i_xyz];
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
    delete[] xyzc_copy;

    return 0;
}
