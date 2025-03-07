#include "../src/periodic_lattice_export.h"
#include "../src/helper.h"
#include "../src/periodic_grid.h"
#include "../src/periodic_becke/periodic_becke.h"

#include <string.h>
#include <algorithm>

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
    lattice.a = 9.4486299386;
    lattice.b = 9.4486299386;
    lattice.c = 9.4486299386;
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
    periodic_parameter.omega = NAN;
    printf("\nomega = %.5e\n", periodic_parameter.omega);

    const int n_atom = 2;
    const double atom_xyz[n_atom * 3] {
        0.0000000000, 0.0000000000, 0.0000000000,
        0.0000000000, 0.0000000000, 4.1951916927,
    };
    const double atom_radius[n_atom] {
        3.0802533600,
        0.7747876550,
    };

    const int n_grid_point = 2114;
    PeriodicBox::GridPoint grid_points[n_grid_point];

    grid_points[   0].x = -16.4716395421; grid_points[   0].y = -16.4716395421; grid_points[   0].z =   0.0000000000; grid_points[   0].w_fixed = 1222.4492551854; grid_points[   0].w_total = 0.0000000000; grid_points[   0].i_atom = 0;
    grid_points[   1].x = -16.4716395421; grid_points[   1].y =   0.0000000000; grid_points[   1].z = -16.4716395421; grid_points[   1].w_fixed = 1222.4492551854; grid_points[   1].w_total = 0.0000000000; grid_points[   1].i_atom = 0;
    grid_points[   2].x = -16.7702682933; grid_points[   2].y =   0.0000000000; grid_points[   2].z =   0.0000000000; grid_points[   2].w_fixed =  238.9289527151; grid_points[   2].w_total = 0.0000000000; grid_points[   2].i_atom = 0;
    grid_points[   3].x =  -7.1135404262; grid_points[   3].y =  -7.1135404262; grid_points[   3].z =  -7.1135404262; grid_points[   3].w_fixed =   69.0587753800; grid_points[   3].w_total = 0.0000000000; grid_points[   3].i_atom = 0;
    grid_points[   4].x =  -6.3448636630; grid_points[   4].y =  -6.3448636630; grid_points[   4].z =   1.9846144277; grid_points[   4].w_fixed =   27.8387029598; grid_points[   4].w_total = 0.0012187121; grid_points[   4].i_atom = 0;
    grid_points[   5].x =  -8.5066870163; grid_points[   5].y =  -8.5066870163; grid_points[   5].z =   2.6608126953; grid_points[   5].w_fixed =   70.1099676542; grid_points[   5].w_total = 0.0000000000; grid_points[   5].i_atom = 0;
    grid_points[   6].x =  -7.0235306998; grid_points[   6].y =  -7.0235306998; grid_points[   6].z =  21.0705920995; grid_points[   6].w_fixed = 1092.3998745831; grid_points[   6].w_total = 0.0000000000; grid_points[   6].i_atom = 0;
    grid_points[   7].x =  -6.3448636630; grid_points[   7].y =   1.9846144277; grid_points[   7].z =  -6.3448636630; grid_points[   7].w_fixed =   27.8387029598; grid_points[   7].w_total = 0.0002691714; grid_points[   7].i_atom = 0;
    grid_points[   8].x =  -8.5066870163; grid_points[   8].y =   2.6608126953; grid_points[   8].z =  -8.5066870163; grid_points[   8].w_fixed =   70.1099676542; grid_points[   8].w_total = 0.0000000000; grid_points[   8].i_atom = 0;
    grid_points[   9].x =  -6.9305700600; grid_points[   9].y =   0.0000000000; grid_points[   9].z =   0.0000000000; grid_points[   9].w_fixed =    4.4485333807; grid_points[   9].w_total = 0.0169320713; grid_points[   9].i_atom = 0;
    grid_points[  10].x =  -6.6888597417; grid_points[  10].y =   1.2829568800; grid_points[  10].z =   1.2829568800; grid_points[  10].w_fixed =    9.5422169009; grid_points[  10].w_total = 0.0783957863; grid_points[  10].i_atom = 0;
    grid_points[  11].x =  -9.1898468013; grid_points[  11].y =   0.0000000000; grid_points[  11].z =   0.0000000000; grid_points[  11].w_fixed =   10.7187035387; grid_points[  11].w_total = 0.0000000000; grid_points[  11].i_atom = 0;
    grid_points[  12].x =  -8.8693420267; grid_points[  12].y =   1.7011843293; grid_points[  12].z =   1.7011843293; grid_points[  12].w_fixed =   22.9918908793; grid_points[  12].w_total = 0.0000018741; grid_points[  12].i_atom = 0;
    grid_points[  13].x =  -7.0235306998; grid_points[  13].y =  21.0705920995; grid_points[  13].z =  -7.0235306998; grid_points[  13].w_fixed = 1092.3998745831; grid_points[  13].w_total = 0.0000000000; grid_points[  13].i_atom = 0;
    grid_points[  14].x =   0.0000000000; grid_points[  14].y = -16.4716395421; grid_points[  14].z = -16.4716395421; grid_points[  14].w_fixed = 1222.4492551854; grid_points[  14].w_total = 0.0000000000; grid_points[  14].i_atom = 0;
    grid_points[  15].x =   0.0000000000; grid_points[  15].y = -16.7702682933; grid_points[  15].z =   0.0000000000; grid_points[  15].w_fixed =  238.9289527151; grid_points[  15].w_total = 0.0000000000; grid_points[  15].i_atom = 0;
    grid_points[  16].x =   1.9846144277; grid_points[  16].y =  -6.3448636630; grid_points[  16].z =  -6.3448636630; grid_points[  16].w_fixed =   27.8387029598; grid_points[  16].w_total = 0.0002691714; grid_points[  16].i_atom = 0;
    grid_points[  17].x =   2.6608126953; grid_points[  17].y =  -8.5066870163; grid_points[  17].z =  -8.5066870163; grid_points[  17].w_fixed =   70.1099676542; grid_points[  17].w_total = 0.0000000000; grid_points[  17].i_atom = 0;
    grid_points[  18].x =   0.0000000000; grid_points[  18].y =  -6.9305700600; grid_points[  18].z =   0.0000000000; grid_points[  18].w_fixed =    4.4485333807; grid_points[  18].w_total = 0.0169320713; grid_points[  18].i_atom = 0;
    grid_points[  19].x =   1.2829568800; grid_points[  19].y =  -6.6888597417; grid_points[  19].z =   1.2829568800; grid_points[  19].w_fixed =    9.5422169009; grid_points[  19].w_total = 0.0783957863; grid_points[  19].i_atom = 0;
    grid_points[  20].x =   0.0000000000; grid_points[  20].y =  -9.1898468013; grid_points[  20].z =   0.0000000000; grid_points[  20].w_fixed =   10.7187035387; grid_points[  20].w_total = 0.0000000000; grid_points[  20].i_atom = 0;
    grid_points[  21].x =   1.7011843293; grid_points[  21].y =  -8.8693420267; grid_points[  21].z =   1.7011843293; grid_points[  21].w_fixed =   22.9918908793; grid_points[  21].w_total = 0.0000018741; grid_points[  21].i_atom = 0;
    grid_points[  22].x =   0.0000000000; grid_points[  22].y =   0.0000000000; grid_points[  22].z = -16.7702682933; grid_points[  22].w_fixed =  238.9289527151; grid_points[  22].w_total = 0.0000000000; grid_points[  22].i_atom = 0;
    grid_points[  23].x =   0.0000000000; grid_points[  23].y =   0.0000000000; grid_points[  23].z =  -6.9305700600; grid_points[  23].w_fixed =    4.4485333807; grid_points[  23].w_total = 0.0000000000; grid_points[  23].i_atom = 0;
    grid_points[  24].x =   1.2829568800; grid_points[  24].y =   1.2829568800; grid_points[  24].z =  -6.6888597417; grid_points[  24].w_fixed =    9.5422169009; grid_points[  24].w_total = 0.0000055082; grid_points[  24].i_atom = 0;
    grid_points[  25].x =   0.0000000000; grid_points[  25].y =   0.0000000000; grid_points[  25].z =  -9.1898468013; grid_points[  25].w_fixed =   10.7187035387; grid_points[  25].w_total = 0.0000000000; grid_points[  25].i_atom = 0;
    grid_points[  26].x =   1.7011843293; grid_points[  26].y =   1.7011843293; grid_points[  26].z =  -8.8693420267; grid_points[  26].w_fixed =   22.9918908793; grid_points[  26].w_total = 0.0000000000; grid_points[  26].i_atom = 0;
    grid_points[  27].x =   0.0036626080; grid_points[  27].y =   0.0000000000; grid_points[  27].z =   0.0000000000; grid_points[  27].w_fixed =    0.0000002129; grid_points[  27].w_total = 0.0000002129; grid_points[  27].i_atom = 0;
    grid_points[  28].x =   0.0157155784; grid_points[  28].y =   0.0000000000; grid_points[  28].z =   0.0000000000; grid_points[  28].w_fixed =    0.0000087099; grid_points[  28].w_total = 0.0000087099; grid_points[  28].i_atom = 0;
    grid_points[  29].x =   0.0380278193; grid_points[  29].y =   0.0000000000; grid_points[  29].z =   0.0000000000; grid_points[  29].w_fixed =    0.0000853157; grid_points[  29].w_total = 0.0000853157; grid_points[  29].i_atom = 0;
    grid_points[  30].x =   0.0729054050; grid_points[  30].y =   0.0000000000; grid_points[  30].z =   0.0000000000; grid_points[  30].w_fixed =    0.0004682262; grid_points[  30].w_total = 0.0004682262; grid_points[  30].i_atom = 0;
    grid_points[  31].x =   0.1232101344; grid_points[  31].y =   0.0000000000; grid_points[  31].z =   0.0000000000; grid_points[  31].w_fixed =    0.0018803520; grid_points[  31].w_total = 0.0018803520; grid_points[  31].i_atom = 0;
    grid_points[  32].x =   0.1925158350; grid_points[  32].y =   0.0000000000; grid_points[  32].z =   0.0000000000; grid_points[  32].w_fixed =    0.0062265394; grid_points[  32].w_total = 0.0062265394; grid_points[  32].i_atom = 0;
    grid_points[  33].x =   0.2853164738; grid_points[  33].y =   0.0000000000; grid_points[  33].z =   0.0000000000; grid_points[  33].w_fixed =    0.0181286111; grid_points[  33].w_total = 0.0181286111; grid_points[  33].i_atom = 0;
    grid_points[  34].x =   0.4073062294; grid_points[  34].y =   0.0000000000; grid_points[  34].z =   0.0000000000; grid_points[  34].w_fixed =    0.0482458934; grid_points[  34].w_total = 0.0482458934; grid_points[  34].i_atom = 0;
    grid_points[  35].x =   0.5657608212; grid_points[  35].y =   0.0000000000; grid_points[  35].z =   0.0000000000; grid_points[  35].w_fixed =    0.0344015986; grid_points[  35].w_total = 0.0344015985; grid_points[  35].i_atom = 0;
    grid_points[  36].x =   0.4000533132; grid_points[  36].y =   0.4000533132; grid_points[  36].z =   0.0000000000; grid_points[  36].w_fixed =    0.0275212788; grid_points[  36].w_total = 0.0275212788; grid_points[  36].i_atom = 0;
    grid_points[  37].x =   0.3266421624; grid_points[  37].y =   0.3266421624; grid_points[  37].z =   0.3266421624; grid_points[  37].w_fixed =    0.0232210790; grid_points[  37].w_total = 0.0232210790; grid_points[  37].i_atom = 0;
    grid_points[  38].x =   0.7700633400; grid_points[  38].y =   0.0000000000; grid_points[  38].z =   0.0000000000; grid_points[  38].w_fixed =    0.0819768394; grid_points[  38].w_total = 0.0819768384; grid_points[  38].i_atom = 0;
    grid_points[  39].x =   0.5445170097; grid_points[  39].y =   0.5445170097; grid_points[  39].z =   0.0000000000; grid_points[  39].w_fixed =    0.0655814715; grid_points[  39].w_total = 0.0655814713; grid_points[  39].i_atom = 0;
    grid_points[  40].x =   0.4445962766; grid_points[  40].y =   0.4445962766; grid_points[  40].z =   0.4445962766; grid_points[  40].w_fixed =    0.0553343666; grid_points[  40].w_total = 0.0553343665; grid_points[  40].i_atom = 0;
    grid_points[  41].x =   1.0324394919; grid_points[  41].y =   0.0000000000; grid_points[  41].z =   0.0000000000; grid_points[  41].w_fixed =    0.1890555184; grid_points[  41].w_total = 0.1890554699; grid_points[  41].i_atom = 0;
    grid_points[  42].x =   0.7300449659; grid_points[  42].y =   0.7300449659; grid_points[  42].z =   0.0000000000; grid_points[  42].w_fixed =    0.1512444148; grid_points[  42].w_total = 0.1512444075; grid_points[  42].i_atom = 0;
    grid_points[  43].x =   0.5960792185; grid_points[  43].y =   0.5960792185; grid_points[  43].z =   0.5960792185; grid_points[  43].w_fixed =    0.1276124749; grid_points[  43].w_total = 0.1276124734; grid_points[  43].i_atom = 0;
    grid_points[  44].x =   1.3690014933; grid_points[  44].y =   0.0000000000; grid_points[  44].z =   0.0000000000; grid_points[  44].w_fixed =    0.4264811437; grid_points[  44].w_total = 0.4264788573; grid_points[  44].i_atom = 0;
    grid_points[  45].x =   0.9680302394; grid_points[  45].y =   0.9680302394; grid_points[  45].z =   0.0000000000; grid_points[  45].w_fixed =    0.3411849149; grid_points[  45].w_total = 0.3411846500; grid_points[  45].i_atom = 0;
    grid_points[  46].x =   0.7903933807; grid_points[  46].y =   0.7903933807; grid_points[  46].z =   0.7903933807; grid_points[  46].w_fixed =    0.2878747720; grid_points[  46].w_total = 0.2878747232; grid_points[  46].i_atom = 0;
    grid_points[  47].x =   1.8012554250; grid_points[  47].y =   0.0000000000; grid_points[  47].z =   0.0000000000; grid_points[  47].w_fixed =    0.2531886212; grid_points[  47].w_total = 0.2531617203; grid_points[  47].i_atom = 0;
    grid_points[  48].x =   1.2736799257; grid_points[  48].y =   1.2736799257; grid_points[  48].z =   0.0000000000; grid_points[  48].w_fixed =    0.4501131044; grid_points[  48].w_total = 0.4501076870; grid_points[  48].i_atom = 0;
    grid_points[  49].x =   1.0399553045; grid_points[  49].y =   1.0399553045; grid_points[  49].z =   1.0399553045; grid_points[  49].w_fixed =    0.4205799265; grid_points[  49].w_total = 0.4205787515; grid_points[  49].i_atom = 0;
    grid_points[  50].x =   0.5430989451; grid_points[  50].y =   0.5430989451; grid_points[  50].z =   1.6292968354; grid_points[  50].w_fixed =    0.4022281471; grid_points[  50].w_total = 0.4021170368; grid_points[  50].i_atom = 0;
    grid_points[  51].x =   0.5430989451; grid_points[  51].y =   1.6292968354; grid_points[  51].z =   0.5430989451; grid_points[  51].w_fixed =    0.4022281471; grid_points[  51].w_total = 0.4022123636; grid_points[  51].i_atom = 0;
    grid_points[  52].x =   1.6292968354; grid_points[  52].y =   0.5430989451; grid_points[  52].z =   0.5430989451; grid_points[  52].w_fixed =    0.4022281471; grid_points[  52].w_total = 0.4022123636; grid_points[  52].i_atom = 0;
    grid_points[  53].x =   2.3583189787; grid_points[  53].y =   0.0000000000; grid_points[  53].z =   0.0000000000; grid_points[  53].w_fixed =    0.5606220456; grid_points[  53].w_total = 0.5595902681; grid_points[  53].i_atom = 0;
    grid_points[  54].x =   1.6675833421; grid_points[  54].y =   1.6675833421; grid_points[  54].z =   0.0000000000; grid_points[  54].w_fixed =    0.9966614143; grid_points[  54].w_total = 0.9964922905; grid_points[  54].i_atom = 0;
    grid_points[  55].x =   1.3615760972; grid_points[  55].y =   1.3615760972; grid_points[  55].z =   1.3615760972; grid_points[  55].w_fixed =    0.9312676753; grid_points[  55].w_total = 0.9312229918; grid_points[  55].i_atom = 0;
    grid_points[  56].x =   0.7110599262; grid_points[  56].y =   0.7110599262; grid_points[  56].z =   2.1331797787; grid_points[  56].w_fixed =    0.8906323100; grid_points[  56].w_total = 0.8798148112; grid_points[  56].i_atom = 0;
    grid_points[  57].x =   0.7110599262; grid_points[  57].y =   2.1331797787; grid_points[  57].z =   0.7110599262; grid_points[  57].w_fixed =    0.8906323100; grid_points[  57].w_total = 0.8900464448; grid_points[  57].i_atom = 0;
    grid_points[  58].x =   2.1331797787; grid_points[  58].y =   0.7110599262; grid_points[  58].z =   0.7110599262; grid_points[  58].w_fixed =    0.8906323100; grid_points[  58].w_total = 0.8900464448; grid_points[  58].i_atom = 0;
    grid_points[  59].x =   1.7783851066; grid_points[  59].y =   1.7783851066; grid_points[  59].z =   1.7783851066; grid_points[  59].w_fixed =    0.9591496581; grid_points[  59].w_total = 0.9586300412; grid_points[  59].i_atom = 0;
    grid_points[  60].x =   2.1266717541; grid_points[  60].y =   2.1266717541; grid_points[  60].z =   0.6652031738; grid_points[  60].w_fixed =    0.9737495508; grid_points[  60].w_total = 0.9722784395; grid_points[  60].i_atom = 0;
    grid_points[  61].x =   2.1266717541; grid_points[  61].y =   0.6652031738; grid_points[  61].z =   2.1266717541; grid_points[  61].w_fixed =    0.9737495508; grid_points[  61].w_total = 0.9690650946; grid_points[  61].i_atom = 0;
    grid_points[  62].x =   0.6652031738; grid_points[  62].y =   2.1266717541; grid_points[  62].z =   2.1266717541; grid_points[  62].w_fixed =    0.9737495508; grid_points[  62].w_total = 0.9690650946; grid_points[  62].i_atom = 0;
    grid_points[  63].x =   1.2188238289; grid_points[  63].y =   1.2188238289; grid_points[  63].z =   2.5528215194; grid_points[  63].w_fixed =    0.9397324606; grid_points[  63].w_total = 0.8796756695; grid_points[  63].i_atom = 0;
    grid_points[  64].x =   1.2188238289; grid_points[  64].y =   2.5528215194; grid_points[  64].z =   1.2188238289; grid_points[  64].w_fixed =    0.9397324606; grid_points[  64].w_total = 0.9361991824; grid_points[  64].i_atom = 0;
    grid_points[  65].x =   2.5528215194; grid_points[  65].y =   1.2188238289; grid_points[  65].z =   1.2188238289; grid_points[  65].w_fixed =    0.9397324606; grid_points[  65].w_total = 0.9361991824; grid_points[  65].i_atom = 0;
    grid_points[  66].x =   1.4734978083; grid_points[  66].y =   2.7049519350; grid_points[  66].z =   0.0000000000; grid_points[  66].w_fixed =    0.9494794437; grid_points[  66].w_total = 0.9424324035; grid_points[  66].i_atom = 0;
    grid_points[  67].x =   2.7049519350; grid_points[  67].y =   1.4734978083; grid_points[  67].z =   0.0000000000; grid_points[  67].w_fixed =    0.9494794437; grid_points[  67].w_total = 0.9424324035; grid_points[  67].i_atom = 0;
    grid_points[  68].x =   2.3227887106; grid_points[  68].y =   2.3227887106; grid_points[  68].z =   2.3227887106; grid_points[  68].w_fixed =    2.1467059397; grid_points[  68].w_total = 2.1405170372; grid_points[  68].i_atom = 0;
    grid_points[  69].x =   0.5702030578; grid_points[  69].y =   0.5702030578; grid_points[  69].z =   2.9728265519; grid_points[  69].w_fixed =    0.8042164697; grid_points[  69].w_total = 0.4464214906; grid_points[  69].i_atom = 0;
    grid_points[  70].x =   2.2808122310; grid_points[  70].y =   2.2808122310; grid_points[  70].z =  11.8913062075; grid_points[  70].w_fixed =   57.9035858166; grid_points[  70].w_total = 0.0000000000; grid_points[  70].i_atom = 0;
    grid_points[  71].x =   0.5702030578; grid_points[  71].y =   2.9728265519; grid_points[  71].z =   0.5702030578; grid_points[  71].w_fixed =    0.8042164697; grid_points[  71].w_total = 0.7903823761; grid_points[  71].i_atom = 0;
    grid_points[  72].x =   0.8688367985; grid_points[  72].y =   2.7776937196; grid_points[  72].z =   2.7776937196; grid_points[  72].w_fixed =    2.1793824633; grid_points[  72].w_total = 2.1192111689; grid_points[  72].i_atom = 0;
    grid_points[  73].x =   2.2808122310; grid_points[  73].y =  11.8913062075; grid_points[  73].z =   2.2808122310; grid_points[  73].w_fixed =   57.9035858166; grid_points[  73].w_total = 0.0000000000; grid_points[  73].i_atom = 0;
    grid_points[  74].x =   3.0802533600; grid_points[  74].y =   0.0000000000; grid_points[  74].z =   0.0000000000; grid_points[  74].w_fixed =    0.3749216611; grid_points[  74].w_total = 0.3657764259; grid_points[  74].i_atom = 0;
    grid_points[  75].x =   2.9728265519; grid_points[  75].y =   0.5702030578; grid_points[  75].z =   0.5702030578; grid_points[  75].w_fixed =    0.8042164697; grid_points[  75].w_total = 0.7903823761; grid_points[  75].i_atom = 0;
    grid_points[  76].x =   2.7776937196; grid_points[  76].y =   0.8688367985; grid_points[  76].z =   2.7776937196; grid_points[  76].w_fixed =    2.1793824633; grid_points[  76].w_total = 2.1192111689; grid_points[  76].i_atom = 0;
    grid_points[  77].x =   2.7776937196; grid_points[  77].y =   2.7776937196; grid_points[  77].z =   0.8688367985; grid_points[  77].w_fixed =    2.1793824633; grid_points[  77].w_total = 2.1482231938; grid_points[  77].i_atom = 0;
    grid_points[  78].x =   3.0411437621; grid_points[  78].y =   3.0411437621; grid_points[  78].z =   3.0411437621; grid_points[  78].w_fixed =    4.8832621034; grid_points[  78].w_total = 4.8025419330; grid_points[  78].i_atom = 0;
    grid_points[  79].x =   9.6823189135; grid_points[  79].y =   9.6823189135; grid_points[  79].z =   9.6823189135; grid_points[  79].w_fixed =  396.8927232114; grid_points[  79].w_total = 0.0000000000; grid_points[  79].i_atom = 0;
    grid_points[  80].x =  10.0265076477; grid_points[  80].y =  10.0265076477; grid_points[  80].z =  30.0795229432; grid_points[  80].w_fixed = 3474.1698452873; grid_points[  80].w_total = 0.0000000000; grid_points[  80].i_atom = 0;
    grid_points[  81].x =  10.0265076477; grid_points[  81].y =  30.0795229432; grid_points[  81].z =  10.0265076477; grid_points[  81].w_fixed = 3474.1698452873; grid_points[  81].w_total = 0.0000000000; grid_points[  81].i_atom = 0;
    grid_points[  82].x =  12.3210134399; grid_points[  82].y =   0.0000000000; grid_points[  82].z =   0.0000000000; grid_points[  82].w_fixed =   26.9943595963; grid_points[  82].w_total = 0.0000000000; grid_points[  82].i_atom = 0;
    grid_points[  83].x =  11.8913062075; grid_points[  83].y =   2.2808122310; grid_points[  83].z =   2.2808122310; grid_points[  83].w_fixed =   57.9035858166; grid_points[  83].w_total = 0.0000000000; grid_points[  83].i_atom = 0;
    grid_points[  84].x =  11.8583704325; grid_points[  84].y =  11.8583704325; grid_points[  84].z =   0.0000000000; grid_points[  84].w_fixed =  424.7625826047; grid_points[  84].w_total = 0.0000000000; grid_points[  84].i_atom = 0;
    grid_points[  85].x =  21.0705920995; grid_points[  85].y =  -7.0235306998; grid_points[  85].z =  -7.0235306998; grid_points[  85].w_fixed = 1092.3998745831; grid_points[  85].w_total = 0.0000000000; grid_points[  85].i_atom = 0;
    grid_points[  86].x =  19.1993004361; grid_points[  86].y =  19.1993004361; grid_points[  86].z =  19.1993004361; grid_points[  86].w_fixed = 3632.6798826243; grid_points[  86].w_total = 0.0000000000; grid_points[  86].i_atom = 0;
    grid_points[  87].x =  30.0795229432; grid_points[  87].y =  10.0265076477; grid_points[  87].z =  10.0265076477; grid_points[  87].w_fixed = 3474.1698452873; grid_points[  87].w_total = 0.0000000000; grid_points[  87].i_atom = 0;
    grid_points[  88].x =  -8.0701458556; grid_points[  88].y =   0.0000000000; grid_points[  88].z =  -4.3961380893; grid_points[  88].w_fixed =   27.1448404564; grid_points[  88].w_total = 0.0000094033; grid_points[  88].i_atom = 0;
    grid_points[  89].x =   0.0000000000; grid_points[  89].y =  -8.0701458556; grid_points[  89].z =  -4.3961380893; grid_points[  89].w_fixed =   27.1448404564; grid_points[  89].w_total = 0.0000094033; grid_points[  89].i_atom = 0;
    grid_points[  90].x =   0.0000000000; grid_points[  90].y =  -8.0701458556; grid_points[  90].z =   4.3961380893; grid_points[  90].w_fixed =   27.1448404564; grid_points[  90].w_total = 0.0000023325; grid_points[  90].i_atom = 0;
    grid_points[  91].x =   0.0000000000; grid_points[  91].y =   0.0000000000; grid_points[  91].z = -33.2541638252; grid_points[  91].w_fixed = 2186.8690181113; grid_points[  91].w_total = 0.0000000000; grid_points[  91].i_atom = 0;
    grid_points[  92].x =   0.0000000000; grid_points[  92].y =   0.0000000000; grid_points[  92].z = -23.2944160349; grid_points[  92].w_fixed =  687.6277060418; grid_points[  92].w_total = 0.0000000000; grid_points[  92].i_atom = 0;
    grid_points[  93].x =   0.0000000000; grid_points[  93].y =   0.0000000000; grid_points[  93].z =  -5.2674155091; grid_points[  93].w_fixed =    1.9088165479; grid_points[  93].w_total = 0.0000000000; grid_points[  93].i_atom = 0;
    grid_points[  94].x =   2.7423536150; grid_points[  94].y =   2.7423536150; grid_points[  94].z =  -5.7438484187; grid_points[  94].w_fixed =   11.1501458946; grid_points[  94].w_total = 0.2040185179; grid_points[  94].i_atom = 0;
    grid_points[  95].x =   0.0000000000; grid_points[  95].y =   0.0000000000; grid_points[  95].z =  -4.0231880620; grid_points[  95].w_fixed =    0.8391251042; grid_points[  95].w_total = 0.1445208656; grid_points[  95].i_atom = 0;
    grid_points[  96].x =   0.7447550142; grid_points[  96].y =   0.7447550142; grid_points[  96].z =  -3.8828754963; grid_points[  96].w_fixed =    1.7999446257; grid_points[  96].w_total = 0.6776187779; grid_points[  96].i_atom = 0;
    grid_points[  97].x =   1.5919331643; grid_points[  97].y =   1.5919331643; grid_points[  97].z =  -3.3342974947; grid_points[  97].w_fixed =    2.1032476402; grid_points[  97].w_total = 1.9796616798; grid_points[  97].i_atom = 0;
    grid_points[  98].x =   1.9245685660; grid_points[  98].y =   0.0000000000; grid_points[  98].z =  -3.5329984457; grid_points[  98].w_fixed =    2.1250626994; grid_points[  98].w_total = 1.7554724438; grid_points[  98].i_atom = 0;
    grid_points[  99].x =   0.0000000000; grid_points[  99].y =   1.9245685660; grid_points[  99].z =  -3.5329984457; grid_points[  99].w_fixed =    2.1250626994; grid_points[  99].w_total = 1.7554724438; grid_points[  99].i_atom = 0;
    grid_points[ 100].x =   0.9750809686; grid_points[ 100].y =   0.9750809686; grid_points[ 100].z =  -5.0837093106; grid_points[ 100].w_fixed =    4.0944598961; grid_points[ 100].w_total = 0.0118697638; grid_points[ 100].i_atom = 0;
    grid_points[ 101].x =   2.0842608672; grid_points[ 101].y =   2.0842608672; grid_points[ 101].z =  -4.3654758527; grid_points[ 101].w_fixed =    4.7844044707; grid_points[ 101].w_total = 2.7579173487; grid_points[ 101].i_atom = 0;
    grid_points[ 102].x =   2.5197684414; grid_points[ 102].y =   0.0000000000; grid_points[ 102].z =  -4.6256278652; grid_points[ 102].w_fixed =    4.8340287112; grid_points[ 102].w_total = 0.9672553950; grid_points[ 102].i_atom = 0;
    grid_points[ 103].x =   0.0000000000; grid_points[ 103].y =   2.5197684414; grid_points[ 103].z =  -4.6256278652; grid_points[ 103].w_fixed =    4.8340287112; grid_points[ 103].w_total = 0.9672553950; grid_points[ 103].i_atom = 0;
    grid_points[ 104].x =   0.7447550142; grid_points[ 104].y =   0.7447550142; grid_points[ 104].z =   3.8828754963; grid_points[ 104].w_fixed =    1.7999446257; grid_points[ 104].w_total = 0.0335908066; grid_points[ 104].i_atom = 0;
    grid_points[ 105].x =   1.5919331643; grid_points[ 105].y =   1.5919331643; grid_points[ 105].z =   3.3342974947; grid_points[ 105].w_fixed =    2.1032476402; grid_points[ 105].w_total = 1.2794689390; grid_points[ 105].i_atom = 0;
    grid_points[ 106].x =   0.9750809686; grid_points[ 106].y =   0.9750809686; grid_points[ 106].z =   5.0837093106; grid_points[ 106].w_fixed =    4.0944598961; grid_points[ 106].w_total = 0.0001490129; grid_points[ 106].i_atom = 0;
    grid_points[ 107].x =   2.0842608672; grid_points[ 107].y =   2.0842608672; grid_points[ 107].z =   4.3654758527; grid_points[ 107].w_fixed =    4.7844044707; grid_points[ 107].w_total = 1.0736736243; grid_points[ 107].i_atom = 0;
    grid_points[ 108].x =   2.7423536150; grid_points[ 108].y =   2.7423536150; grid_points[ 108].z =   5.7438484187; grid_points[ 108].w_fixed =   11.1501458946; grid_points[ 108].w_total = 0.0554063878; grid_points[ 108].i_atom = 0;
    grid_points[ 109].x =  10.8198077400; grid_points[ 109].y =   0.0000000000; grid_points[ 109].z =  -5.8939912332; grid_points[ 109].w_fixed =   68.3625199466; grid_points[ 109].w_total = 0.0000000000; grid_points[ 109].i_atom = 0;
    grid_points[ 110].x =  -7.0235306998; grid_points[ 110].y =  -7.0235306998; grid_points[ 110].z = -21.0705920995; grid_points[ 110].w_fixed = 1092.3998745831; grid_points[ 110].w_total = 0.0000000000; grid_points[ 110].i_atom = 0;
    grid_points[ 111].x =  -6.3448636630; grid_points[ 111].y =  -6.3448636630; grid_points[ 111].z =  -1.9846144277; grid_points[ 111].w_fixed =   27.8387029598; grid_points[ 111].w_total = 0.0013808480; grid_points[ 111].i_atom = 0;
    grid_points[ 112].x =  -8.5066870163; grid_points[ 112].y =  -8.5066870163; grid_points[ 112].z =  -2.6608126953; grid_points[ 112].w_fixed =   70.1099676542; grid_points[ 112].w_total = 0.0000000000; grid_points[ 112].i_atom = 0;
    grid_points[ 113].x =  -7.1135404262; grid_points[ 113].y =  -7.1135404262; grid_points[ 113].z =   7.1135404262; grid_points[ 113].w_fixed =   69.0587753800; grid_points[ 113].w_total = 0.0000000000; grid_points[ 113].i_atom = 0;
    grid_points[ 114].x =  -6.6888597417; grid_points[ 114].y =   1.2829568800; grid_points[ 114].z =  -1.2829568800; grid_points[ 114].w_fixed =    9.5422169009; grid_points[ 114].w_total = 0.0865328484; grid_points[ 114].i_atom = 0;
    grid_points[ 115].x =  -8.8693420267; grid_points[ 115].y =   1.7011843293; grid_points[ 115].z =  -1.7011843293; grid_points[ 115].w_fixed =   22.9918908793; grid_points[ 115].w_total = 0.0000042871; grid_points[ 115].i_atom = 0;
    grid_points[ 116].x =  -6.3448636630; grid_points[ 116].y =   1.9846144277; grid_points[ 116].z =   6.3448636630; grid_points[ 116].w_fixed =   27.8387029598; grid_points[ 116].w_total = 0.0001758536; grid_points[ 116].i_atom = 0;
    grid_points[ 117].x =  -8.5066870163; grid_points[ 117].y =   2.6608126953; grid_points[ 117].z =   8.5066870163; grid_points[ 117].w_fixed =   70.1099676542; grid_points[ 117].w_total = 0.0000000000; grid_points[ 117].i_atom = 0;
    grid_points[ 118].x =  -7.0235306998; grid_points[ 118].y =  21.0705920995; grid_points[ 118].z =   7.0235306998; grid_points[ 118].w_fixed = 1092.3998745831; grid_points[ 118].w_total = 0.0000000000; grid_points[ 118].i_atom = 0;
    grid_points[ 119].x =   0.0000000000; grid_points[ 119].y = -16.4716395421; grid_points[ 119].z =  16.4716395421; grid_points[ 119].w_fixed = 1222.4492551854; grid_points[ 119].w_total = 0.0000000000; grid_points[ 119].i_atom = 0;
    grid_points[ 120].x =   1.2829568800; grid_points[ 120].y =  -6.6888597417; grid_points[ 120].z =  -1.2829568800; grid_points[ 120].w_fixed =    9.5422169009; grid_points[ 120].w_total = 0.0865328484; grid_points[ 120].i_atom = 0;
    grid_points[ 121].x =   1.7011843293; grid_points[ 121].y =  -8.8693420267; grid_points[ 121].z =  -1.7011843293; grid_points[ 121].w_fixed =   22.9918908793; grid_points[ 121].w_total = 0.0000042871; grid_points[ 121].i_atom = 0;
    grid_points[ 122].x =   1.9846144277; grid_points[ 122].y =  -6.3448636630; grid_points[ 122].z =   6.3448636630; grid_points[ 122].w_fixed =   27.8387029598; grid_points[ 122].w_total = 0.0001758536; grid_points[ 122].i_atom = 0;
    grid_points[ 123].x =   2.6608126953; grid_points[ 123].y =  -8.5066870163; grid_points[ 123].z =   8.5066870163; grid_points[ 123].w_fixed =   70.1099676542; grid_points[ 123].w_total = 0.0000000000; grid_points[ 123].i_atom = 0;
    grid_points[ 124].x =   0.0000000000; grid_points[ 124].y =   0.0000000000; grid_points[ 124].z = -12.3210134399; grid_points[ 124].w_fixed =   26.9943595963; grid_points[ 124].w_total = 0.0000000000; grid_points[ 124].i_atom = 0;
    grid_points[ 125].x =   2.2808122310; grid_points[ 125].y =   2.2808122310; grid_points[ 125].z = -11.8913062075; grid_points[ 125].w_fixed =   57.9035858166; grid_points[ 125].w_total = 0.0000000000; grid_points[ 125].i_atom = 0;
    grid_points[ 126].x =   0.0000000000; grid_points[ 126].y =   0.0000000000; grid_points[ 126].z =  -1.3690014933; grid_points[ 126].w_fixed =    0.4264811437; grid_points[ 126].w_total = 0.4264811272; grid_points[ 126].i_atom = 0;
    grid_points[ 127].x =   0.0000000000; grid_points[ 127].y =   0.0000000000; grid_points[ 127].z =  -1.8012554250; grid_points[ 127].w_fixed =    0.2531886212; grid_points[ 127].w_total = 0.2531861451; grid_points[ 127].i_atom = 0;
    grid_points[ 128].x =   0.0000000000; grid_points[ 128].y =   1.2736799257; grid_points[ 128].z =  -1.2736799257; grid_points[ 128].w_fixed =    0.4501131044; grid_points[ 128].w_total = 0.4501105971; grid_points[ 128].i_atom = 0;
    grid_points[ 129].x =   1.2736799257; grid_points[ 129].y =   0.0000000000; grid_points[ 129].z =  -1.2736799257; grid_points[ 129].w_fixed =    0.4501131044; grid_points[ 129].w_total = 0.4501105971; grid_points[ 129].i_atom = 0;
    grid_points[ 130].x =   0.5430989451; grid_points[ 130].y =   0.5430989451; grid_points[ 130].z =  -1.6292968354; grid_points[ 130].w_fixed =    0.4022281471; grid_points[ 130].w_total = 0.4022275009; grid_points[ 130].i_atom = 0;
    grid_points[ 131].x =   0.0000000000; grid_points[ 131].y =   0.0000000000; grid_points[ 131].z =  -2.3583189787; grid_points[ 131].w_fixed =    0.5606220456; grid_points[ 131].w_total = 0.5598142760; grid_points[ 131].i_atom = 0;
    grid_points[ 132].x =   0.0000000000; grid_points[ 132].y =   1.6675833421; grid_points[ 132].z =  -1.6675833421; grid_points[ 132].w_fixed =    0.9966614143; grid_points[ 132].w_total = 0.9965883126; grid_points[ 132].i_atom = 0;
    grid_points[ 133].x =   1.6675833421; grid_points[ 133].y =   0.0000000000; grid_points[ 133].z =  -1.6675833421; grid_points[ 133].w_fixed =    0.9966614143; grid_points[ 133].w_total = 0.9965883126; grid_points[ 133].i_atom = 0;
    grid_points[ 134].x =   1.3615760972; grid_points[ 134].y =   1.3615760972; grid_points[ 134].z =  -1.3615760972; grid_points[ 134].w_fixed =    0.9312676753; grid_points[ 134].w_total = 0.9312379628; grid_points[ 134].i_atom = 0;
    grid_points[ 135].x =   0.7110599262; grid_points[ 135].y =   0.7110599262; grid_points[ 135].z =  -2.1331797787; grid_points[ 135].w_fixed =    0.8906323100; grid_points[ 135].w_total = 0.8904412859; grid_points[ 135].i_atom = 0;
    grid_points[ 136].x =   0.0000000000; grid_points[ 136].y =   0.0000000000; grid_points[ 136].z =  -3.0802533600; grid_points[ 136].w_fixed =    0.3749216611; grid_points[ 136].w_total = 0.3434002227; grid_points[ 136].i_atom = 0;
    grid_points[ 137].x =   1.7783851066; grid_points[ 137].y =   1.7783851066; grid_points[ 137].z =  -1.7783851066; grid_points[ 137].w_fixed =    0.9591496581; grid_points[ 137].w_total = 0.9588589425; grid_points[ 137].i_atom = 0;
    grid_points[ 138].x =   0.5702030578; grid_points[ 138].y =   0.5702030578; grid_points[ 138].z =  -2.9728265519; grid_points[ 138].w_fixed =    0.8042164697; grid_points[ 138].w_total = 0.7676826133; grid_points[ 138].i_atom = 0;
    grid_points[ 139].x =   2.1266717541; grid_points[ 139].y =   0.6652031738; grid_points[ 139].z =  -2.1266717541; grid_points[ 139].w_fixed =    0.9737495508; grid_points[ 139].w_total = 0.9730669536; grid_points[ 139].i_atom = 0;
    grid_points[ 140].x =   0.6652031738; grid_points[ 140].y =   2.1266717541; grid_points[ 140].z =  -2.1266717541; grid_points[ 140].w_fixed =    0.9737495508; grid_points[ 140].w_total = 0.9730669536; grid_points[ 140].i_atom = 0;
    grid_points[ 141].x =   1.2188238289; grid_points[ 141].y =   1.2188238289; grid_points[ 141].z =  -2.5528215194; grid_points[ 141].w_fixed =    0.9397324606; grid_points[ 141].w_total = 0.9373660293; grid_points[ 141].i_atom = 0;
    grid_points[ 142].x =   1.4734978083; grid_points[ 142].y =   0.0000000000; grid_points[ 142].z =  -2.7049519350; grid_points[ 142].w_fixed =    0.9494794437; grid_points[ 142].w_total = 0.9420711494; grid_points[ 142].i_atom = 0;
    grid_points[ 143].x =   2.7049519350; grid_points[ 143].y =   0.0000000000; grid_points[ 143].z =  -1.4734978083; grid_points[ 143].w_fixed =    0.9494794437; grid_points[ 143].w_total = 0.9427515891; grid_points[ 143].i_atom = 0;
    grid_points[ 144].x =   0.0000000000; grid_points[ 144].y =   1.4734978083; grid_points[ 144].z =  -2.7049519350; grid_points[ 144].w_fixed =    0.9494794437; grid_points[ 144].w_total = 0.9420711494; grid_points[ 144].i_atom = 0;
    grid_points[ 145].x =   2.3227887106; grid_points[ 145].y =   2.3227887106; grid_points[ 145].z =  -2.3227887106; grid_points[ 145].w_fixed =    2.1467059397; grid_points[ 145].w_total = 2.1417068903; grid_points[ 145].i_atom = 0;
    grid_points[ 146].x =   0.0000000000; grid_points[ 146].y =   0.0000000000; grid_points[ 146].z =  -0.0036626080; grid_points[ 146].w_fixed =    0.0000002129; grid_points[ 146].w_total = 0.0000002129; grid_points[ 146].i_atom = 0;
    grid_points[ 147].x =   0.0000000000; grid_points[ 147].y =   0.0000000000; grid_points[ 147].z =  -0.0157155784; grid_points[ 147].w_fixed =    0.0000087099; grid_points[ 147].w_total = 0.0000087099; grid_points[ 147].i_atom = 0;
    grid_points[ 148].x =   0.0000000000; grid_points[ 148].y =   0.0000000000; grid_points[ 148].z =  -0.0380278193; grid_points[ 148].w_fixed =    0.0000853157; grid_points[ 148].w_total = 0.0000853157; grid_points[ 148].i_atom = 0;
    grid_points[ 149].x =   0.0000000000; grid_points[ 149].y =   0.0000000000; grid_points[ 149].z =  -0.0729054050; grid_points[ 149].w_fixed =    0.0004682262; grid_points[ 149].w_total = 0.0004682262; grid_points[ 149].i_atom = 0;
    grid_points[ 150].x =   0.0000000000; grid_points[ 150].y =   0.0000000000; grid_points[ 150].z =  -0.1232101344; grid_points[ 150].w_fixed =    0.0018803520; grid_points[ 150].w_total = 0.0018803520; grid_points[ 150].i_atom = 0;
    grid_points[ 151].x =   0.0000000000; grid_points[ 151].y =   0.0000000000; grid_points[ 151].z =  -0.1925158350; grid_points[ 151].w_fixed =    0.0062265394; grid_points[ 151].w_total = 0.0062265394; grid_points[ 151].i_atom = 0;
    grid_points[ 152].x =   0.0000000000; grid_points[ 152].y =   0.0000000000; grid_points[ 152].z =  -0.2853164738; grid_points[ 152].w_fixed =    0.0181286111; grid_points[ 152].w_total = 0.0181286111; grid_points[ 152].i_atom = 0;
    grid_points[ 153].x =   0.0000000000; grid_points[ 153].y =   0.0000000000; grid_points[ 153].z =  -0.4073062294; grid_points[ 153].w_fixed =    0.0482458934; grid_points[ 153].w_total = 0.0482458934; grid_points[ 153].i_atom = 0;
    grid_points[ 154].x =   0.0000000000; grid_points[ 154].y =   0.0000000000; grid_points[ 154].z =  -0.5657608212; grid_points[ 154].w_fixed =    0.0344015986; grid_points[ 154].w_total = 0.0344015986; grid_points[ 154].i_atom = 0;
    grid_points[ 155].x =   0.0000000000; grid_points[ 155].y =   0.4000533132; grid_points[ 155].z =  -0.4000533132; grid_points[ 155].w_fixed =    0.0275212788; grid_points[ 155].w_total = 0.0275212788; grid_points[ 155].i_atom = 0;
    grid_points[ 156].x =   0.4000533132; grid_points[ 156].y =   0.0000000000; grid_points[ 156].z =  -0.4000533132; grid_points[ 156].w_fixed =    0.0275212788; grid_points[ 156].w_total = 0.0275212788; grid_points[ 156].i_atom = 0;
    grid_points[ 157].x =   0.3266421624; grid_points[ 157].y =   0.3266421624; grid_points[ 157].z =  -0.3266421624; grid_points[ 157].w_fixed =    0.0232210790; grid_points[ 157].w_total = 0.0232210790; grid_points[ 157].i_atom = 0;
    grid_points[ 158].x =   0.0000000000; grid_points[ 158].y =   0.0000000000; grid_points[ 158].z =  -0.7700633400; grid_points[ 158].w_fixed =    0.0819768394; grid_points[ 158].w_total = 0.0819768394; grid_points[ 158].i_atom = 0;
    grid_points[ 159].x =   0.0000000000; grid_points[ 159].y =   0.5445170097; grid_points[ 159].z =  -0.5445170097; grid_points[ 159].w_fixed =    0.0655814715; grid_points[ 159].w_total = 0.0655814714; grid_points[ 159].i_atom = 0;
    grid_points[ 160].x =   0.5445170097; grid_points[ 160].y =   0.0000000000; grid_points[ 160].z =  -0.5445170097; grid_points[ 160].w_fixed =    0.0655814715; grid_points[ 160].w_total = 0.0655814714; grid_points[ 160].i_atom = 0;
    grid_points[ 161].x =   0.4445962766; grid_points[ 161].y =   0.4445962766; grid_points[ 161].z =  -0.4445962766; grid_points[ 161].w_fixed =    0.0553343666; grid_points[ 161].w_total = 0.0553343665; grid_points[ 161].i_atom = 0;
    grid_points[ 162].x =   0.0000000000; grid_points[ 162].y =   0.0000000000; grid_points[ 162].z =  -1.0324394919; grid_points[ 162].w_fixed =    0.1890555184; grid_points[ 162].w_total = 0.1890555184; grid_points[ 162].i_atom = 0;
    grid_points[ 163].x =   0.0000000000; grid_points[ 163].y =   0.7300449659; grid_points[ 163].z =  -0.7300449659; grid_points[ 163].w_fixed =    0.1512444148; grid_points[ 163].w_total = 0.1512444110; grid_points[ 163].i_atom = 0;
    grid_points[ 164].x =   0.7300449659; grid_points[ 164].y =   0.0000000000; grid_points[ 164].z =  -0.7300449659; grid_points[ 164].w_fixed =    0.1512444148; grid_points[ 164].w_total = 0.1512444110; grid_points[ 164].i_atom = 0;
    grid_points[ 165].x =   0.5960792185; grid_points[ 165].y =   0.5960792185; grid_points[ 165].z =  -0.5960792185; grid_points[ 165].w_fixed =    0.1276124749; grid_points[ 165].w_total = 0.1276124729; grid_points[ 165].i_atom = 0;
    grid_points[ 166].x =   0.0000000000; grid_points[ 166].y =   0.9680302394; grid_points[ 166].z =  -0.9680302394; grid_points[ 166].w_fixed =    0.3411849149; grid_points[ 166].w_total = 0.3411847831; grid_points[ 166].i_atom = 0;
    grid_points[ 167].x =   0.9680302394; grid_points[ 167].y =   0.0000000000; grid_points[ 167].z =  -0.9680302394; grid_points[ 167].w_fixed =    0.3411849149; grid_points[ 167].w_total = 0.3411847831; grid_points[ 167].i_atom = 0;
    grid_points[ 168].x =   0.7903933807; grid_points[ 168].y =   0.7903933807; grid_points[ 168].z =  -0.7903933807; grid_points[ 168].w_fixed =    0.2878747720; grid_points[ 168].w_total = 0.2878747082; grid_points[ 168].i_atom = 0;
    grid_points[ 169].x =   1.0399553045; grid_points[ 169].y =   1.0399553045; grid_points[ 169].z =  -1.0399553045; grid_points[ 169].w_fixed =    0.4205799265; grid_points[ 169].w_total = 0.4205787392; grid_points[ 169].i_atom = 0;
    grid_points[ 170].x =   0.5430989451; grid_points[ 170].y =   1.6292968354; grid_points[ 170].z =  -0.5430989451; grid_points[ 170].w_fixed =    0.4022281471; grid_points[ 170].w_total = 0.4022099166; grid_points[ 170].i_atom = 0;
    grid_points[ 171].x =   1.6292968354; grid_points[ 171].y =   0.5430989451; grid_points[ 171].z =  -0.5430989451; grid_points[ 171].w_fixed =    0.4022281471; grid_points[ 171].w_total = 0.4022099166; grid_points[ 171].i_atom = 0;
    grid_points[ 172].x =   0.7110599262; grid_points[ 172].y =   2.1331797787; grid_points[ 172].z =  -0.7110599262; grid_points[ 172].w_fixed =    0.8906323100; grid_points[ 172].w_total = 0.8899740741; grid_points[ 172].i_atom = 0;
    grid_points[ 173].x =   2.1331797787; grid_points[ 173].y =   0.7110599262; grid_points[ 173].z =  -0.7110599262; grid_points[ 173].w_fixed =    0.8906323100; grid_points[ 173].w_total = 0.8899740741; grid_points[ 173].i_atom = 0;
    grid_points[ 174].x =   2.1266717541; grid_points[ 174].y =   2.1266717541; grid_points[ 174].z =  -0.6652031738; grid_points[ 174].w_fixed =    0.9737495508; grid_points[ 174].w_total = 0.9721447960; grid_points[ 174].i_atom = 0;
    grid_points[ 175].x =   1.2188238289; grid_points[ 175].y =   2.5528215194; grid_points[ 175].z =  -1.2188238289; grid_points[ 175].w_fixed =    0.9397324606; grid_points[ 175].w_total = 0.9357421258; grid_points[ 175].i_atom = 0;
    grid_points[ 176].x =   2.5528215194; grid_points[ 176].y =   1.2188238289; grid_points[ 176].z =  -1.2188238289; grid_points[ 176].w_fixed =    0.9397324606; grid_points[ 176].w_total = 0.9357421258; grid_points[ 176].i_atom = 0;
    grid_points[ 177].x =   1.2829568800; grid_points[ 177].y =   1.2829568800; grid_points[ 177].z =   6.6888597417; grid_points[ 177].w_fixed =    9.5422169009; grid_points[ 177].w_total = 0.0000003249; grid_points[ 177].i_atom = 0;
    grid_points[ 178].x =   1.7011843293; grid_points[ 178].y =   1.7011843293; grid_points[ 178].z =   8.8693420267; grid_points[ 178].w_fixed =   22.9918908793; grid_points[ 178].w_total = 0.0000000000; grid_points[ 178].i_atom = 0;
    grid_points[ 179].x =   0.8688367985; grid_points[ 179].y =   2.7776937196; grid_points[ 179].z =  -2.7776937196; grid_points[ 179].w_fixed =    2.1793824633; grid_points[ 179].w_total = 2.1646215031; grid_points[ 179].i_atom = 0;
    grid_points[ 180].x =   0.5702030578; grid_points[ 180].y =   2.9728265519; grid_points[ 180].z =  -0.5702030578; grid_points[ 180].w_fixed =    0.8042164697; grid_points[ 180].w_total = 0.7899201355; grid_points[ 180].i_atom = 0;
    grid_points[ 181].x =   0.0000000000; grid_points[ 181].y =  11.8583704325; grid_points[ 181].z = -11.8583704325; grid_points[ 181].w_fixed =  424.7625826047; grid_points[ 181].w_total = 0.0000000000; grid_points[ 181].i_atom = 0;
    grid_points[ 182].x =   2.2808122310; grid_points[ 182].y =  11.8913062075; grid_points[ 182].z =  -2.2808122310; grid_points[ 182].w_fixed =   57.9035858166; grid_points[ 182].w_total = 0.0000000000; grid_points[ 182].i_atom = 0;
    grid_points[ 183].x =   2.7776937196; grid_points[ 183].y =   0.8688367985; grid_points[ 183].z =  -2.7776937196; grid_points[ 183].w_fixed =    2.1793824633; grid_points[ 183].w_total = 2.1646215031; grid_points[ 183].i_atom = 0;
    grid_points[ 184].x =   2.9728265519; grid_points[ 184].y =   0.5702030578; grid_points[ 184].z =  -0.5702030578; grid_points[ 184].w_fixed =    0.8042164697; grid_points[ 184].w_total = 0.7899201355; grid_points[ 184].i_atom = 0;
    grid_points[ 185].x =   3.0411437621; grid_points[ 185].y =   3.0411437621; grid_points[ 185].z =  -3.0411437621; grid_points[ 185].w_fixed =    4.8832621034; grid_points[ 185].w_total = 4.8017896257; grid_points[ 185].i_atom = 0;
    grid_points[ 186].x =   2.7776937196; grid_points[ 186].y =   2.7776937196; grid_points[ 186].z =  -0.8688367985; grid_points[ 186].w_fixed =    2.1793824633; grid_points[ 186].w_total = 2.1466862309; grid_points[ 186].i_atom = 0;
    grid_points[ 187].x =  10.0265076477; grid_points[ 187].y =  10.0265076477; grid_points[ 187].z = -30.0795229432; grid_points[ 187].w_fixed = 3474.1698452873; grid_points[ 187].w_total = 0.0000000000; grid_points[ 187].i_atom = 0;
    grid_points[ 188].x =   9.6823189135; grid_points[ 188].y =   9.6823189135; grid_points[ 188].z =  -9.6823189135; grid_points[ 188].w_fixed =  396.8927232114; grid_points[ 188].w_total = 0.0000000000; grid_points[ 188].i_atom = 0;
    grid_points[ 189].x =  10.0265076477; grid_points[ 189].y =  30.0795229432; grid_points[ 189].z = -10.0265076477; grid_points[ 189].w_fixed = 3474.1698452873; grid_points[ 189].w_total = 0.0000000000; grid_points[ 189].i_atom = 0;
    grid_points[ 190].x =  11.8583704325; grid_points[ 190].y =   0.0000000000; grid_points[ 190].z = -11.8583704325; grid_points[ 190].w_fixed =  424.7625826047; grid_points[ 190].w_total = 0.0000000000; grid_points[ 190].i_atom = 0;
    grid_points[ 191].x =  11.8913062075; grid_points[ 191].y =   2.2808122310; grid_points[ 191].z =  -2.2808122310; grid_points[ 191].w_fixed =   57.9035858166; grid_points[ 191].w_total = 0.0000000000; grid_points[ 191].i_atom = 0;
    grid_points[ 192].x =  21.0705920995; grid_points[ 192].y =  -7.0235306998; grid_points[ 192].z =   7.0235306998; grid_points[ 192].w_fixed = 1092.3998745831; grid_points[ 192].w_total = 0.0000000000; grid_points[ 192].i_atom = 0;
    grid_points[ 193].x =  19.1993004361; grid_points[ 193].y =  19.1993004361; grid_points[ 193].z = -19.1993004361; grid_points[ 193].w_fixed = 3632.6798826243; grid_points[ 193].w_total = 0.0000000000; grid_points[ 193].i_atom = 0;
    grid_points[ 194].x =  30.0795229432; grid_points[ 194].y =  10.0265076477; grid_points[ 194].z = -10.0265076477; grid_points[ 194].w_fixed = 3474.1698452873; grid_points[ 194].w_total = 0.0000000000; grid_points[ 194].i_atom = 0;
    grid_points[ 195].x =  -8.0701458556; grid_points[ 195].y =  -4.3961380893; grid_points[ 195].z =   0.0000000000; grid_points[ 195].w_fixed =   27.1448404564; grid_points[ 195].w_total = 0.0003003544; grid_points[ 195].i_atom = 0;
    grid_points[ 196].x =  -8.0701458556; grid_points[ 196].y =   4.3961380893; grid_points[ 196].z =   0.0000000000; grid_points[ 196].w_fixed =   27.1448404564; grid_points[ 196].w_total = 0.0003003544; grid_points[ 196].i_atom = 0;
    grid_points[ 197].x =   0.0000000000; grid_points[ 197].y = -33.2541638252; grid_points[ 197].z =   0.0000000000; grid_points[ 197].w_fixed = 2186.8690181113; grid_points[ 197].w_total = 0.0000000000; grid_points[ 197].i_atom = 0;
    grid_points[ 198].x =   0.0000000000; grid_points[ 198].y = -23.2944160349; grid_points[ 198].z =   0.0000000000; grid_points[ 198].w_fixed =  687.6277060418; grid_points[ 198].w_total = 0.0000000000; grid_points[ 198].i_atom = 0;
    grid_points[ 199].x =   0.0000000000; grid_points[ 199].y =  -5.2674155091; grid_points[ 199].z =   0.0000000000; grid_points[ 199].w_fixed =    1.9088165479; grid_points[ 199].w_total = 0.4997305992; grid_points[ 199].i_atom = 0;
    grid_points[ 200].x =   2.7423536150; grid_points[ 200].y =  -5.7438484187; grid_points[ 200].z =   2.7423536150; grid_points[ 200].w_fixed =   11.1501458946; grid_points[ 200].w_total = 0.7753792790; grid_points[ 200].i_atom = 0;
    grid_points[ 201].x =   0.0000000000; grid_points[ 201].y =  -4.3961380893; grid_points[ 201].z =  -8.0701458556; grid_points[ 201].w_fixed =   27.1448404564; grid_points[ 201].w_total = 0.0000024398; grid_points[ 201].i_atom = 0;
    grid_points[ 202].x =   0.0000000000; grid_points[ 202].y =  -4.0231880620; grid_points[ 202].z =   0.0000000000; grid_points[ 202].w_fixed =    0.8391251042; grid_points[ 202].w_total = 0.6689355801; grid_points[ 202].i_atom = 0;
    grid_points[ 203].x =   0.7447550142; grid_points[ 203].y =  -3.8828754963; grid_points[ 203].z =   0.7447550142; grid_points[ 203].w_fixed =    1.7999446257; grid_points[ 203].w_total = 1.5204246946; grid_points[ 203].i_atom = 0;
    grid_points[ 204].x =   1.5919331643; grid_points[ 204].y =  -3.3342974947; grid_points[ 204].z =   1.5919331643; grid_points[ 204].w_fixed =    2.1032476402; grid_points[ 204].w_total = 2.0194520110; grid_points[ 204].i_atom = 0;
    grid_points[ 205].x =   1.9245685660; grid_points[ 205].y =  -3.5329984457; grid_points[ 205].z =   0.0000000000; grid_points[ 205].w_fixed =    2.1250626994; grid_points[ 205].w_total = 1.9718171965; grid_points[ 205].i_atom = 0;
    grid_points[ 206].x =   0.0000000000; grid_points[ 206].y =  -3.5329984457; grid_points[ 206].z =   1.9245685660; grid_points[ 206].w_fixed =    2.1250626994; grid_points[ 206].w_total = 1.9844543174; grid_points[ 206].i_atom = 0;
    grid_points[ 207].x =   0.9750809686; grid_points[ 207].y =  -5.0837093106; grid_points[ 207].z =   0.9750809686; grid_points[ 207].w_fixed =    4.0944598961; grid_points[ 207].w_total = 1.3604510062; grid_points[ 207].i_atom = 0;
    grid_points[ 208].x =   2.0842608672; grid_points[ 208].y =  -4.3654758527; grid_points[ 208].z =   2.0842608672; grid_points[ 208].w_fixed =    4.7844044707; grid_points[ 208].w_total = 3.3058550217; grid_points[ 208].i_atom = 0;
    grid_points[ 209].x =   2.5197684414; grid_points[ 209].y =  -4.6256278652; grid_points[ 209].z =   0.0000000000; grid_points[ 209].w_fixed =    4.8340287112; grid_points[ 209].w_total = 2.6697241981; grid_points[ 209].i_atom = 0;
    grid_points[ 210].x =   0.0000000000; grid_points[ 210].y =  -4.6256278652; grid_points[ 210].z =   2.5197684414; grid_points[ 210].w_fixed =    4.8340287112; grid_points[ 210].w_total = 2.6802181871; grid_points[ 210].i_atom = 0;
    grid_points[ 211].x =   0.0000000000; grid_points[ 211].y =   4.3961380893; grid_points[ 211].z =  -8.0701458556; grid_points[ 211].w_fixed =   27.1448404564; grid_points[ 211].w_total = 0.0000024398; grid_points[ 211].i_atom = 0;
    grid_points[ 212].x =   0.7447550142; grid_points[ 212].y =   3.8828754963; grid_points[ 212].z =   0.7447550142; grid_points[ 212].w_fixed =    1.7999446257; grid_points[ 212].w_total = 1.5204246946; grid_points[ 212].i_atom = 0;
    grid_points[ 213].x =   1.5919331643; grid_points[ 213].y =   3.3342974947; grid_points[ 213].z =   1.5919331643; grid_points[ 213].w_fixed =    2.1032476402; grid_points[ 213].w_total = 2.0194520110; grid_points[ 213].i_atom = 0;
    grid_points[ 214].x =   1.9245685660; grid_points[ 214].y =   3.5329984457; grid_points[ 214].z =   0.0000000000; grid_points[ 214].w_fixed =    2.1250626994; grid_points[ 214].w_total = 1.9718171965; grid_points[ 214].i_atom = 0;
    grid_points[ 215].x =   0.9750809686; grid_points[ 215].y =   5.0837093106; grid_points[ 215].z =   0.9750809686; grid_points[ 215].w_fixed =    4.0944598961; grid_points[ 215].w_total = 1.3604510062; grid_points[ 215].i_atom = 0;
    grid_points[ 216].x =   2.0842608672; grid_points[ 216].y =   4.3654758527; grid_points[ 216].z =   2.0842608672; grid_points[ 216].w_fixed =    4.7844044707; grid_points[ 216].w_total = 3.3058550217; grid_points[ 216].i_atom = 0;
    grid_points[ 217].x =   2.5197684414; grid_points[ 217].y =   4.6256278652; grid_points[ 217].z =   0.0000000000; grid_points[ 217].w_fixed =    4.8340287112; grid_points[ 217].w_total = 2.6697241981; grid_points[ 217].i_atom = 0;
    grid_points[ 218].x =   2.7423536150; grid_points[ 218].y =   5.7438484187; grid_points[ 218].z =   2.7423536150; grid_points[ 218].w_fixed =   11.1501458946; grid_points[ 218].w_total = 0.7753792790; grid_points[ 218].i_atom = 0;
    grid_points[ 219].x =  10.8198077400; grid_points[ 219].y =  -5.8939912332; grid_points[ 219].z =   0.0000000000; grid_points[ 219].w_fixed =   68.3625199466; grid_points[ 219].w_total = 0.0000000000; grid_points[ 219].i_atom = 0;
    grid_points[ 220].x =  10.8198077400; grid_points[ 220].y =   5.8939912332; grid_points[ 220].z =   0.0000000000; grid_points[ 220].w_fixed =   68.3625199466; grid_points[ 220].w_total = 0.0000000000; grid_points[ 220].i_atom = 0;
    grid_points[ 221].x =  -7.6162691612; grid_points[ 221].y =  -3.6363256383; grid_points[ 221].z =  -3.6363256383; grid_points[ 221].w_fixed =   26.8661821836; grid_points[ 221].w_total = 0.0010773582; grid_points[ 221].i_atom = 0;
    grid_points[ 222].x =  -7.6162691612; grid_points[ 222].y =  -3.6363256383; grid_points[ 222].z =   3.6363256383; grid_points[ 222].w_fixed =   26.8661821836; grid_points[ 222].w_total = 0.0009330210; grid_points[ 222].i_atom = 0;
    grid_points[ 223].x =  -7.6162691612; grid_points[ 223].y =   3.6363256383; grid_points[ 223].z =  -3.6363256383; grid_points[ 223].w_fixed =   26.8661821836; grid_points[ 223].w_total = 0.0010773582; grid_points[ 223].i_atom = 0;
    grid_points[ 224].x =  -7.6162691612; grid_points[ 224].y =   3.6363256383; grid_points[ 224].z =   3.6363256383; grid_points[ 224].w_fixed =   26.8661821836; grid_points[ 224].w_total = 0.0009330210; grid_points[ 224].i_atom = 0;
    grid_points[ 225].x =   0.0000000000; grid_points[ 225].y = -23.5142447435; grid_points[ 225].z = -23.5142447435; grid_points[ 225].w_fixed = 3887.7671433090; grid_points[ 225].w_total = 0.0000000000; grid_points[ 225].i_atom = 0;
    grid_points[ 226].x =   0.0000000000; grid_points[ 226].y = -23.5142447435; grid_points[ 226].z =  23.5142447435; grid_points[ 226].w_fixed = 3887.7671433090; grid_points[ 226].w_total = 0.0000000000; grid_points[ 226].i_atom = 0;
    grid_points[ 227].x =   0.0000000000; grid_points[ 227].y =  -6.0861418537; grid_points[ 227].z =  -3.3153700687; grid_points[ 227].w_fixed =   11.2657961338; grid_points[ 227].w_total = 0.3549590929; grid_points[ 227].i_atom = 0;
    grid_points[ 228].x =   0.0000000000; grid_points[ 228].y =  -6.0861418537; grid_points[ 228].z =   3.3153700687; grid_points[ 228].w_fixed =   11.2657961338; grid_points[ 228].w_total = 0.3471224190; grid_points[ 228].i_atom = 0;
    grid_points[ 229].x =   0.0000000000; grid_points[ 229].y =  -3.3153700687; grid_points[ 229].z =  -6.0861418537; grid_points[ 229].w_fixed =   11.2657961338; grid_points[ 229].w_total = 0.0242256985; grid_points[ 229].i_atom = 0;
    grid_points[ 230].x =   1.1375367884; grid_points[ 230].y =  -3.6367345381; grid_points[ 230].z =  -3.6367345381; grid_points[ 230].w_fixed =    4.9575936763; grid_points[ 230].w_total = 4.6559094346; grid_points[ 230].i_atom = 0;
    grid_points[ 231].x =   1.4967071411; grid_points[ 231].y =  -4.7850114467; grid_points[ 231].z =  -4.7850114467; grid_points[ 231].w_fixed =   11.5537666422; grid_points[ 231].w_total = 2.9255841962; grid_points[ 231].i_atom = 0;
    grid_points[ 232].x =   1.1375367884; grid_points[ 232].y =  -3.6367345381; grid_points[ 232].z =   3.6367345381; grid_points[ 232].w_fixed =    4.9575936763; grid_points[ 232].w_total = 4.5221022069; grid_points[ 232].i_atom = 0;
    grid_points[ 233].x =   1.4967071411; grid_points[ 233].y =  -4.7850114467; grid_points[ 233].z =   4.7850114467; grid_points[ 233].w_fixed =   11.5537666422; grid_points[ 233].w_total = 1.8853671084; grid_points[ 233].i_atom = 0;
    grid_points[ 234].x =   0.0000000000; grid_points[ 234].y =   3.3153700687; grid_points[ 234].z =  -6.0861418537; grid_points[ 234].w_fixed =   11.2657961338; grid_points[ 234].w_total = 0.0242256985; grid_points[ 234].i_atom = 0;
    grid_points[ 235].x =   1.1375367884; grid_points[ 235].y =   3.6367345381; grid_points[ 235].z =  -3.6367345381; grid_points[ 235].w_fixed =    4.9575936763; grid_points[ 235].w_total = 4.6559094346; grid_points[ 235].i_atom = 0;
    grid_points[ 236].x =   1.4967071411; grid_points[ 236].y =   4.7850114467; grid_points[ 236].z =  -4.7850114467; grid_points[ 236].w_fixed =   11.5537666422; grid_points[ 236].w_total = 2.9255841962; grid_points[ 236].i_atom = 0;
    grid_points[ 237].x =   1.1375367884; grid_points[ 237].y =   3.6367345381; grid_points[ 237].z =   3.6367345381; grid_points[ 237].w_fixed =    4.9575936763; grid_points[ 237].w_total = 4.5221022069; grid_points[ 237].i_atom = 0;
    grid_points[ 238].x =   1.4967071411; grid_points[ 238].y =   4.7850114467; grid_points[ 238].z =   4.7850114467; grid_points[ 238].w_fixed =   11.5537666422; grid_points[ 238].w_total = 1.8853671084; grid_points[ 238].i_atom = 0;
    grid_points[ 239].x =   0.0000000000; grid_points[ 239].y =  23.5142447435; grid_points[ 239].z = -23.5142447435; grid_points[ 239].w_fixed = 3887.7671433090; grid_points[ 239].w_total = 0.0000000000; grid_points[ 239].i_atom = 0;
    grid_points[ 240].x =  10.2112860776; grid_points[ 240].y =  -4.8752953156; grid_points[ 240].z =  -4.8752953156; grid_points[ 240].w_fixed =   67.6607371617; grid_points[ 240].w_total = 0.0000000000; grid_points[ 240].i_atom = 0;
    grid_points[ 241].x =  10.2112860776; grid_points[ 241].y =  -4.8752953156; grid_points[ 241].z =   4.8752953156; grid_points[ 241].w_fixed =   67.6607371617; grid_points[ 241].w_total = 0.0000000000; grid_points[ 241].i_atom = 0;
    grid_points[ 242].x =  10.2112860776; grid_points[ 242].y =   4.8752953156; grid_points[ 242].z =  -4.8752953156; grid_points[ 242].w_fixed =   67.6607371617; grid_points[ 242].w_total = 0.0000000000; grid_points[ 242].i_atom = 0;
    grid_points[ 243].x =  10.2112860776; grid_points[ 243].y =   4.8752953156; grid_points[ 243].z =   4.8752953156; grid_points[ 243].w_fixed =   67.6607371617; grid_points[ 243].w_total = 0.0000000000; grid_points[ 243].i_atom = 0;
    grid_points[ 244].x =   0.0000000000; grid_points[ 244].y =  -5.8939912332; grid_points[ 244].z = -10.8198077400; grid_points[ 244].w_fixed =   68.3625199466; grid_points[ 244].w_total = 0.0000000000; grid_points[ 244].i_atom = 0;
    grid_points[ 245].x =   2.7423536150; grid_points[ 245].y =  -5.7438484187; grid_points[ 245].z =  -2.7423536150; grid_points[ 245].w_fixed =   11.1501458946; grid_points[ 245].w_total = 0.7975816561; grid_points[ 245].i_atom = 0;
    grid_points[ 246].x =   1.5919331643; grid_points[ 246].y =  -3.3342974947; grid_points[ 246].z =  -1.5919331643; grid_points[ 246].w_fixed =    2.1032476402; grid_points[ 246].w_total = 2.0149126278; grid_points[ 246].i_atom = 0;
    grid_points[ 247].x =   0.0000000000; grid_points[ 247].y =  -3.5329984457; grid_points[ 247].z =  -1.9245685660; grid_points[ 247].w_fixed =    2.1250626994; grid_points[ 247].w_total = 1.9769819275; grid_points[ 247].i_atom = 0;
    grid_points[ 248].x =   2.0842608672; grid_points[ 248].y =  -4.3654758527; grid_points[ 248].z =  -2.0842608672; grid_points[ 248].w_fixed =    4.7844044707; grid_points[ 248].w_total = 3.2966811514; grid_points[ 248].i_atom = 0;
    grid_points[ 249].x =   0.0000000000; grid_points[ 249].y =  -4.6256278652; grid_points[ 249].z =  -2.5197684414; grid_points[ 249].w_fixed =    4.8340287112; grid_points[ 249].w_total = 2.6774240028; grid_points[ 249].i_atom = 0;
    grid_points[ 250].x =   0.7447550142; grid_points[ 250].y =  -3.8828754963; grid_points[ 250].z =  -0.7447550142; grid_points[ 250].w_fixed =    1.7999446257; grid_points[ 250].w_total = 1.5182141465; grid_points[ 250].i_atom = 0;
    grid_points[ 251].x =   0.9750809686; grid_points[ 251].y =  -5.0837093106; grid_points[ 251].z =  -0.9750809686; grid_points[ 251].w_fixed =    4.0944598961; grid_points[ 251].w_total = 1.3639957134; grid_points[ 251].i_atom = 0;
    grid_points[ 252].x =   0.0000000000; grid_points[ 252].y =   5.8939912332; grid_points[ 252].z = -10.8198077400; grid_points[ 252].w_fixed =   68.3625199466; grid_points[ 252].w_total = 0.0000000000; grid_points[ 252].i_atom = 0;
    grid_points[ 253].x =   1.5919331643; grid_points[ 253].y =   3.3342974947; grid_points[ 253].z =  -1.5919331643; grid_points[ 253].w_fixed =    2.1032476402; grid_points[ 253].w_total = 2.0149126278; grid_points[ 253].i_atom = 0;
    grid_points[ 254].x =   2.0842608672; grid_points[ 254].y =   4.3654758527; grid_points[ 254].z =  -2.0842608672; grid_points[ 254].w_fixed =    4.7844044707; grid_points[ 254].w_total = 3.2966811514; grid_points[ 254].i_atom = 0;
    grid_points[ 255].x =   2.7423536150; grid_points[ 255].y =   5.7438484187; grid_points[ 255].z =  -2.7423536150; grid_points[ 255].w_fixed =   11.1501458946; grid_points[ 255].w_total = 0.7975816561; grid_points[ 255].i_atom = 0;
    grid_points[ 256].x =   0.7447550142; grid_points[ 256].y =   3.8828754963; grid_points[ 256].z =  -0.7447550142; grid_points[ 256].w_fixed =    1.7999446257; grid_points[ 256].w_total = 1.5182141465; grid_points[ 256].i_atom = 0;
    grid_points[ 257].x =   0.9750809686; grid_points[ 257].y =   5.0837093106; grid_points[ 257].z =  -0.9750809686; grid_points[ 257].w_fixed =    4.0944598961; grid_points[ 257].w_total = 1.3639957134; grid_points[ 257].i_atom = 0;
    grid_points[ 258].x = -16.4716395421; grid_points[ 258].y =  16.4716395421; grid_points[ 258].z =   0.0000000000; grid_points[ 258].w_fixed = 1222.4492551854; grid_points[ 258].w_total = 0.0000000000; grid_points[ 258].i_atom = 0;
    grid_points[ 259].x =  -7.0235306998; grid_points[ 259].y = -21.0705920995; grid_points[ 259].z =  -7.0235306998; grid_points[ 259].w_fixed = 1092.3998745831; grid_points[ 259].w_total = 0.0000000000; grid_points[ 259].i_atom = 0;
    grid_points[ 260].x =  -6.3448636630; grid_points[ 260].y =  -1.9846144277; grid_points[ 260].z =  -6.3448636630; grid_points[ 260].w_fixed =   27.8387029598; grid_points[ 260].w_total = 0.0002691714; grid_points[ 260].i_atom = 0;
    grid_points[ 261].x =  -8.5066870163; grid_points[ 261].y =  -2.6608126953; grid_points[ 261].z =  -8.5066870163; grid_points[ 261].w_fixed =   70.1099676542; grid_points[ 261].w_total = 0.0000000000; grid_points[ 261].i_atom = 0;
    grid_points[ 262].x =  -6.6888597417; grid_points[ 262].y =  -1.2829568800; grid_points[ 262].z =   1.2829568800; grid_points[ 262].w_fixed =    9.5422169009; grid_points[ 262].w_total = 0.0783957863; grid_points[ 262].i_atom = 0;
    grid_points[ 263].x =  -8.8693420267; grid_points[ 263].y =  -1.7011843293; grid_points[ 263].z =   1.7011843293; grid_points[ 263].w_fixed =   22.9918908793; grid_points[ 263].w_total = 0.0000018741; grid_points[ 263].i_atom = 0;
    grid_points[ 264].x =  -6.3448636630; grid_points[ 264].y =   6.3448636630; grid_points[ 264].z =   1.9846144277; grid_points[ 264].w_fixed =   27.8387029598; grid_points[ 264].w_total = 0.0012187121; grid_points[ 264].i_atom = 0;
    grid_points[ 265].x =  -7.1135404262; grid_points[ 265].y =   7.1135404262; grid_points[ 265].z =  -7.1135404262; grid_points[ 265].w_fixed =   69.0587753800; grid_points[ 265].w_total = 0.0000000000; grid_points[ 265].i_atom = 0;
    grid_points[ 266].x =  -8.5066870163; grid_points[ 266].y =   8.5066870163; grid_points[ 266].z =   2.6608126953; grid_points[ 266].w_fixed =   70.1099676542; grid_points[ 266].w_total = 0.0000000000; grid_points[ 266].i_atom = 0;
    grid_points[ 267].x =  -7.0235306998; grid_points[ 267].y =   7.0235306998; grid_points[ 267].z =  21.0705920995; grid_points[ 267].w_fixed = 1092.3998745831; grid_points[ 267].w_total = 0.0000000000; grid_points[ 267].i_atom = 0;
    grid_points[ 268].x =   0.0000000000; grid_points[ 268].y = -12.3210134399; grid_points[ 268].z =   0.0000000000; grid_points[ 268].w_fixed =   26.9943595963; grid_points[ 268].w_total = 0.0000000000; grid_points[ 268].i_atom = 0;
    grid_points[ 269].x =   2.2808122310; grid_points[ 269].y = -11.8913062075; grid_points[ 269].z =   2.2808122310; grid_points[ 269].w_fixed =   57.9035858166; grid_points[ 269].w_total = 0.0000000000; grid_points[ 269].i_atom = 0;
    grid_points[ 270].x =   0.0000000000; grid_points[ 270].y = -11.8583704325; grid_points[ 270].z =  11.8583704325; grid_points[ 270].w_fixed =  424.7625826047; grid_points[ 270].w_total = 0.0000000000; grid_points[ 270].i_atom = 0;
    grid_points[ 271].x =   1.2829568800; grid_points[ 271].y =  -1.2829568800; grid_points[ 271].z =  -6.6888597417; grid_points[ 271].w_fixed =    9.5422169009; grid_points[ 271].w_total = 0.0000055082; grid_points[ 271].i_atom = 0;
    grid_points[ 272].x =   1.7011843293; grid_points[ 272].y =  -1.7011843293; grid_points[ 272].z =  -8.8693420267; grid_points[ 272].w_fixed =   22.9918908793; grid_points[ 272].w_total = 0.0000000000; grid_points[ 272].i_atom = 0;
    grid_points[ 273].x =   0.0000000000; grid_points[ 273].y =  -1.3690014933; grid_points[ 273].z =   0.0000000000; grid_points[ 273].w_fixed =    0.4264811437; grid_points[ 273].w_total = 0.4264788573; grid_points[ 273].i_atom = 0;
    grid_points[ 274].x =   0.0000000000; grid_points[ 274].y =  -1.8012554250; grid_points[ 274].z =   0.0000000000; grid_points[ 274].w_fixed =    0.2531886212; grid_points[ 274].w_total = 0.2531617203; grid_points[ 274].i_atom = 0;
    grid_points[ 275].x =   0.0000000000; grid_points[ 275].y =  -1.2736799257; grid_points[ 275].z =   1.2736799257; grid_points[ 275].w_fixed =    0.4501131044; grid_points[ 275].w_total = 0.4501069125; grid_points[ 275].i_atom = 0;
    grid_points[ 276].x =   1.2736799257; grid_points[ 276].y =  -1.2736799257; grid_points[ 276].z =   0.0000000000; grid_points[ 276].w_fixed =    0.4501131044; grid_points[ 276].w_total = 0.4501076870; grid_points[ 276].i_atom = 0;
    grid_points[ 277].x =   0.5430989451; grid_points[ 277].y =  -1.6292968354; grid_points[ 277].z =   0.5430989451; grid_points[ 277].w_fixed =    0.4022281471; grid_points[ 277].w_total = 0.4022123636; grid_points[ 277].i_atom = 0;
    grid_points[ 278].x =   0.0000000000; grid_points[ 278].y =  -2.3583189787; grid_points[ 278].z =   0.0000000000; grid_points[ 278].w_fixed =    0.5606220456; grid_points[ 278].w_total = 0.5595902681; grid_points[ 278].i_atom = 0;
    grid_points[ 279].x =   0.0000000000; grid_points[ 279].y =  -1.6675833421; grid_points[ 279].z =   1.6675833421; grid_points[ 279].w_fixed =    0.9966614143; grid_points[ 279].w_total = 0.9962426478; grid_points[ 279].i_atom = 0;
    grid_points[ 280].x =   1.6675833421; grid_points[ 280].y =  -1.6675833421; grid_points[ 280].z =   0.0000000000; grid_points[ 280].w_fixed =    0.9966614143; grid_points[ 280].w_total = 0.9964922905; grid_points[ 280].i_atom = 0;
    grid_points[ 281].x =   1.3615760972; grid_points[ 281].y =  -1.3615760972; grid_points[ 281].z =   1.3615760972; grid_points[ 281].w_fixed =    0.9312676753; grid_points[ 281].w_total = 0.9312229918; grid_points[ 281].i_atom = 0;
    grid_points[ 282].x =   0.7110599262; grid_points[ 282].y =  -2.1331797787; grid_points[ 282].z =   0.7110599262; grid_points[ 282].w_fixed =    0.8906323100; grid_points[ 282].w_total = 0.8900464448; grid_points[ 282].i_atom = 0;
    grid_points[ 283].x =   0.0000000000; grid_points[ 283].y =  -3.0802533600; grid_points[ 283].z =   0.0000000000; grid_points[ 283].w_fixed =    0.3749216611; grid_points[ 283].w_total = 0.3657764259; grid_points[ 283].i_atom = 0;
    grid_points[ 284].x =   1.7783851066; grid_points[ 284].y =  -1.7783851066; grid_points[ 284].z =   1.7783851066; grid_points[ 284].w_fixed =    0.9591496581; grid_points[ 284].w_total = 0.9586300412; grid_points[ 284].i_atom = 0;
    grid_points[ 285].x =   0.5702030578; grid_points[ 285].y =  -2.9728265519; grid_points[ 285].z =   0.5702030578; grid_points[ 285].w_fixed =    0.8042164697; grid_points[ 285].w_total = 0.7903823761; grid_points[ 285].i_atom = 0;
    grid_points[ 286].x =   2.1266717541; grid_points[ 286].y =  -2.1266717541; grid_points[ 286].z =   0.6652031738; grid_points[ 286].w_fixed =    0.9737495508; grid_points[ 286].w_total = 0.9722784395; grid_points[ 286].i_atom = 0;
    grid_points[ 287].x =   0.6652031738; grid_points[ 287].y =  -2.1266717541; grid_points[ 287].z =   2.1266717541; grid_points[ 287].w_fixed =    0.9737495508; grid_points[ 287].w_total = 0.9690650946; grid_points[ 287].i_atom = 0;
    grid_points[ 288].x =   1.2188238289; grid_points[ 288].y =  -2.5528215194; grid_points[ 288].z =   1.2188238289; grid_points[ 288].w_fixed =    0.9397324606; grid_points[ 288].w_total = 0.9361991824; grid_points[ 288].i_atom = 0;
    grid_points[ 289].x =   1.4734978083; grid_points[ 289].y =  -2.7049519350; grid_points[ 289].z =   0.0000000000; grid_points[ 289].w_fixed =    0.9494794437; grid_points[ 289].w_total = 0.9424324035; grid_points[ 289].i_atom = 0;
    grid_points[ 290].x =   2.7049519350; grid_points[ 290].y =  -1.4734978083; grid_points[ 290].z =   0.0000000000; grid_points[ 290].w_fixed =    0.9494794437; grid_points[ 290].w_total = 0.9424324035; grid_points[ 290].i_atom = 0;
    grid_points[ 291].x =   0.0000000000; grid_points[ 291].y =  -2.7049519350; grid_points[ 291].z =   1.4734978083; grid_points[ 291].w_fixed =    0.9494794437; grid_points[ 291].w_total = 0.9435891428; grid_points[ 291].i_atom = 0;
    grid_points[ 292].x =   2.3227887106; grid_points[ 292].y =  -2.3227887106; grid_points[ 292].z =   2.3227887106; grid_points[ 292].w_fixed =    2.1467059397; grid_points[ 292].w_total = 2.1405170372; grid_points[ 292].i_atom = 0;
    grid_points[ 293].x =   0.8688367985; grid_points[ 293].y =  -2.7776937196; grid_points[ 293].z =   2.7776937196; grid_points[ 293].w_fixed =    2.1793824633; grid_points[ 293].w_total = 2.1192111689; grid_points[ 293].i_atom = 0;
    grid_points[ 294].x =   2.2808122310; grid_points[ 294].y =  -2.2808122310; grid_points[ 294].z =  11.8913062075; grid_points[ 294].w_fixed =   57.9035858166; grid_points[ 294].w_total = 0.0000000000; grid_points[ 294].i_atom = 0;
    grid_points[ 295].x =   0.0000000000; grid_points[ 295].y =  -0.0036626080; grid_points[ 295].z =   0.0000000000; grid_points[ 295].w_fixed =    0.0000002129; grid_points[ 295].w_total = 0.0000002129; grid_points[ 295].i_atom = 0;
    grid_points[ 296].x =   0.0000000000; grid_points[ 296].y =  -0.0157155784; grid_points[ 296].z =   0.0000000000; grid_points[ 296].w_fixed =    0.0000087099; grid_points[ 296].w_total = 0.0000087099; grid_points[ 296].i_atom = 0;
    grid_points[ 297].x =   0.0000000000; grid_points[ 297].y =  -0.0380278193; grid_points[ 297].z =   0.0000000000; grid_points[ 297].w_fixed =    0.0000853157; grid_points[ 297].w_total = 0.0000853157; grid_points[ 297].i_atom = 0;
    grid_points[ 298].x =   0.0000000000; grid_points[ 298].y =  -0.0729054050; grid_points[ 298].z =   0.0000000000; grid_points[ 298].w_fixed =    0.0004682262; grid_points[ 298].w_total = 0.0004682262; grid_points[ 298].i_atom = 0;
    grid_points[ 299].x =   0.0000000000; grid_points[ 299].y =  -0.1232101344; grid_points[ 299].z =   0.0000000000; grid_points[ 299].w_fixed =    0.0018803520; grid_points[ 299].w_total = 0.0018803520; grid_points[ 299].i_atom = 0;
    grid_points[ 300].x =   0.0000000000; grid_points[ 300].y =  -0.1925158350; grid_points[ 300].z =   0.0000000000; grid_points[ 300].w_fixed =    0.0062265394; grid_points[ 300].w_total = 0.0062265394; grid_points[ 300].i_atom = 0;
    grid_points[ 301].x =   0.0000000000; grid_points[ 301].y =  -0.2853164738; grid_points[ 301].z =   0.0000000000; grid_points[ 301].w_fixed =    0.0181286111; grid_points[ 301].w_total = 0.0181286111; grid_points[ 301].i_atom = 0;
    grid_points[ 302].x =   0.0000000000; grid_points[ 302].y =  -0.4073062294; grid_points[ 302].z =   0.0000000000; grid_points[ 302].w_fixed =    0.0482458934; grid_points[ 302].w_total = 0.0482458934; grid_points[ 302].i_atom = 0;
    grid_points[ 303].x =   0.0000000000; grid_points[ 303].y =  -0.5657608212; grid_points[ 303].z =   0.0000000000; grid_points[ 303].w_fixed =    0.0344015986; grid_points[ 303].w_total = 0.0344015985; grid_points[ 303].i_atom = 0;
    grid_points[ 304].x =   0.0000000000; grid_points[ 304].y =  -0.4000533132; grid_points[ 304].z =   0.4000533132; grid_points[ 304].w_fixed =    0.0275212788; grid_points[ 304].w_total = 0.0275212788; grid_points[ 304].i_atom = 0;
    grid_points[ 305].x =   0.4000533132; grid_points[ 305].y =   0.0000000000; grid_points[ 305].z =   0.4000533132; grid_points[ 305].w_fixed =    0.0275212788; grid_points[ 305].w_total = 0.0275212788; grid_points[ 305].i_atom = 0;
    grid_points[ 306].x =   0.4000533132; grid_points[ 306].y =  -0.4000533132; grid_points[ 306].z =   0.0000000000; grid_points[ 306].w_fixed =    0.0275212788; grid_points[ 306].w_total = 0.0275212788; grid_points[ 306].i_atom = 0;
    grid_points[ 307].x =   0.3266421624; grid_points[ 307].y =  -0.3266421624; grid_points[ 307].z =   0.3266421624; grid_points[ 307].w_fixed =    0.0232210790; grid_points[ 307].w_total = 0.0232210790; grid_points[ 307].i_atom = 0;
    grid_points[ 308].x =   0.0000000000; grid_points[ 308].y =  -0.7700633400; grid_points[ 308].z =   0.0000000000; grid_points[ 308].w_fixed =    0.0819768394; grid_points[ 308].w_total = 0.0819768384; grid_points[ 308].i_atom = 0;
    grid_points[ 309].x =   0.0000000000; grid_points[ 309].y =  -0.5445170097; grid_points[ 309].z =   0.5445170097; grid_points[ 309].w_fixed =    0.0655814715; grid_points[ 309].w_total = 0.0655814714; grid_points[ 309].i_atom = 0;
    grid_points[ 310].x =   0.5445170097; grid_points[ 310].y =   0.0000000000; grid_points[ 310].z =   0.5445170097; grid_points[ 310].w_fixed =    0.0655814715; grid_points[ 310].w_total = 0.0655814714; grid_points[ 310].i_atom = 0;
    grid_points[ 311].x =   0.5445170097; grid_points[ 311].y =  -0.5445170097; grid_points[ 311].z =   0.0000000000; grid_points[ 311].w_fixed =    0.0655814715; grid_points[ 311].w_total = 0.0655814713; grid_points[ 311].i_atom = 0;
    grid_points[ 312].x =   0.4445962766; grid_points[ 312].y =  -0.4445962766; grid_points[ 312].z =   0.4445962766; grid_points[ 312].w_fixed =    0.0553343666; grid_points[ 312].w_total = 0.0553343665; grid_points[ 312].i_atom = 0;
    grid_points[ 313].x =   0.0000000000; grid_points[ 313].y =  -1.0324394919; grid_points[ 313].z =   0.0000000000; grid_points[ 313].w_fixed =    0.1890555184; grid_points[ 313].w_total = 0.1890554699; grid_points[ 313].i_atom = 0;
    grid_points[ 314].x =   0.0000000000; grid_points[ 314].y =  -0.7300449659; grid_points[ 314].z =   0.7300449659; grid_points[ 314].w_fixed =    0.1512444148; grid_points[ 314].w_total = 0.1512444119; grid_points[ 314].i_atom = 0;
    grid_points[ 315].x =   0.7300449659; grid_points[ 315].y =   0.0000000000; grid_points[ 315].z =   0.7300449659; grid_points[ 315].w_fixed =    0.1512444148; grid_points[ 315].w_total = 0.1512444119; grid_points[ 315].i_atom = 0;
    grid_points[ 316].x =   0.7300449659; grid_points[ 316].y =  -0.7300449659; grid_points[ 316].z =   0.0000000000; grid_points[ 316].w_fixed =    0.1512444148; grid_points[ 316].w_total = 0.1512444075; grid_points[ 316].i_atom = 0;
    grid_points[ 317].x =   0.5960792185; grid_points[ 317].y =  -0.5960792185; grid_points[ 317].z =   0.5960792185; grid_points[ 317].w_fixed =    0.1276124749; grid_points[ 317].w_total = 0.1276124734; grid_points[ 317].i_atom = 0;
    grid_points[ 318].x =   0.0000000000; grid_points[ 318].y =  -0.9680302394; grid_points[ 318].z =   0.9680302394; grid_points[ 318].w_fixed =    0.3411849149; grid_points[ 318].w_total = 0.3411847747; grid_points[ 318].i_atom = 0;
    grid_points[ 319].x =   0.9680302394; grid_points[ 319].y =   0.0000000000; grid_points[ 319].z =   0.9680302394; grid_points[ 319].w_fixed =    0.3411849149; grid_points[ 319].w_total = 0.3411847747; grid_points[ 319].i_atom = 0;
    grid_points[ 320].x =   0.9680302394; grid_points[ 320].y =  -0.9680302394; grid_points[ 320].z =   0.0000000000; grid_points[ 320].w_fixed =    0.3411849149; grid_points[ 320].w_total = 0.3411846500; grid_points[ 320].i_atom = 0;
    grid_points[ 321].x =   0.7903933807; grid_points[ 321].y =  -0.7903933807; grid_points[ 321].z =   0.7903933807; grid_points[ 321].w_fixed =    0.2878747720; grid_points[ 321].w_total = 0.2878747232; grid_points[ 321].i_atom = 0;
    grid_points[ 322].x =   1.2736799257; grid_points[ 322].y =   0.0000000000; grid_points[ 322].z =   1.2736799257; grid_points[ 322].w_fixed =    0.4501131044; grid_points[ 322].w_total = 0.4501069125; grid_points[ 322].i_atom = 0;
    grid_points[ 323].x =   1.0399553045; grid_points[ 323].y =  -1.0399553045; grid_points[ 323].z =   1.0399553045; grid_points[ 323].w_fixed =    0.4205799265; grid_points[ 323].w_total = 0.4205787515; grid_points[ 323].i_atom = 0;
    grid_points[ 324].x =   0.5430989451; grid_points[ 324].y =  -0.5430989451; grid_points[ 324].z =   1.6292968354; grid_points[ 324].w_fixed =    0.4022281471; grid_points[ 324].w_total = 0.4021170368; grid_points[ 324].i_atom = 0;
    grid_points[ 325].x =   1.6292968354; grid_points[ 325].y =  -0.5430989451; grid_points[ 325].z =   0.5430989451; grid_points[ 325].w_fixed =    0.4022281471; grid_points[ 325].w_total = 0.4022123636; grid_points[ 325].i_atom = 0;
    grid_points[ 326].x =   1.6675833421; grid_points[ 326].y =   0.0000000000; grid_points[ 326].z =   1.6675833421; grid_points[ 326].w_fixed =    0.9966614143; grid_points[ 326].w_total = 0.9962426478; grid_points[ 326].i_atom = 0;
    grid_points[ 327].x =   0.7110599262; grid_points[ 327].y =  -0.7110599262; grid_points[ 327].z =   2.1331797787; grid_points[ 327].w_fixed =    0.8906323100; grid_points[ 327].w_total = 0.8798148112; grid_points[ 327].i_atom = 0;
    grid_points[ 328].x =   2.1331797787; grid_points[ 328].y =  -0.7110599262; grid_points[ 328].z =   0.7110599262; grid_points[ 328].w_fixed =    0.8906323100; grid_points[ 328].w_total = 0.8900464448; grid_points[ 328].i_atom = 0;
    grid_points[ 329].x =   2.1266717541; grid_points[ 329].y =  -0.6652031738; grid_points[ 329].z =   2.1266717541; grid_points[ 329].w_fixed =    0.9737495508; grid_points[ 329].w_total = 0.9690650946; grid_points[ 329].i_atom = 0;
    grid_points[ 330].x =   1.2188238289; grid_points[ 330].y =  -1.2188238289; grid_points[ 330].z =   2.5528215194; grid_points[ 330].w_fixed =    0.9397324606; grid_points[ 330].w_total = 0.8796756695; grid_points[ 330].i_atom = 0;
    grid_points[ 331].x =   2.5528215194; grid_points[ 331].y =  -1.2188238289; grid_points[ 331].z =   1.2188238289; grid_points[ 331].w_fixed =    0.9397324606; grid_points[ 331].w_total = 0.9361991824; grid_points[ 331].i_atom = 0;
    grid_points[ 332].x =   1.4734978083; grid_points[ 332].y =   0.0000000000; grid_points[ 332].z =   2.7049519350; grid_points[ 332].w_fixed =    0.9494794437; grid_points[ 332].w_total = 0.8169727814; grid_points[ 332].i_atom = 0;
    grid_points[ 333].x =   2.7049519350; grid_points[ 333].y =   0.0000000000; grid_points[ 333].z =   1.4734978083; grid_points[ 333].w_fixed =    0.9494794437; grid_points[ 333].w_total = 0.9435891428; grid_points[ 333].i_atom = 0;
    grid_points[ 334].x =   0.5702030578; grid_points[ 334].y =  -0.5702030578; grid_points[ 334].z =   2.9728265519; grid_points[ 334].w_fixed =    0.8042164697; grid_points[ 334].w_total = 0.4464214906; grid_points[ 334].i_atom = 0;
    grid_points[ 335].x =   1.9846144277; grid_points[ 335].y =   6.3448636630; grid_points[ 335].z =  -6.3448636630; grid_points[ 335].w_fixed =   27.8387029598; grid_points[ 335].w_total = 0.0002691714; grid_points[ 335].i_atom = 0;
    grid_points[ 336].x =   1.2829568800; grid_points[ 336].y =   6.6888597417; grid_points[ 336].z =   1.2829568800; grid_points[ 336].w_fixed =    9.5422169009; grid_points[ 336].w_total = 0.0783957863; grid_points[ 336].i_atom = 0;
    grid_points[ 337].x =   2.6608126953; grid_points[ 337].y =   8.5066870163; grid_points[ 337].z =  -8.5066870163; grid_points[ 337].w_fixed =   70.1099676542; grid_points[ 337].w_total = 0.0000000000; grid_points[ 337].i_atom = 0;
    grid_points[ 338].x =   1.7011843293; grid_points[ 338].y =   8.8693420267; grid_points[ 338].z =   1.7011843293; grid_points[ 338].w_fixed =   22.9918908793; grid_points[ 338].w_total = 0.0000018741; grid_points[ 338].i_atom = 0;
    grid_points[ 339].x =   0.0000000000; grid_points[ 339].y =  16.4716395421; grid_points[ 339].z = -16.4716395421; grid_points[ 339].w_fixed = 1222.4492551854; grid_points[ 339].w_total = 0.0000000000; grid_points[ 339].i_atom = 0;
    grid_points[ 340].x =   2.7776937196; grid_points[ 340].y =  -2.7776937196; grid_points[ 340].z =   0.8688367985; grid_points[ 340].w_fixed =    2.1793824633; grid_points[ 340].w_total = 2.1482231938; grid_points[ 340].i_atom = 0;
    grid_points[ 341].x =   3.0411437621; grid_points[ 341].y =  -3.0411437621; grid_points[ 341].z =   3.0411437621; grid_points[ 341].w_fixed =    4.8832621034; grid_points[ 341].w_total = 4.8025419330; grid_points[ 341].i_atom = 0;
    grid_points[ 342].x =   2.9728265519; grid_points[ 342].y =  -0.5702030578; grid_points[ 342].z =   0.5702030578; grid_points[ 342].w_fixed =    0.8042164697; grid_points[ 342].w_total = 0.7903823761; grid_points[ 342].i_atom = 0;
    grid_points[ 343].x =   2.7776937196; grid_points[ 343].y =  -0.8688367985; grid_points[ 343].z =   2.7776937196; grid_points[ 343].w_fixed =    2.1793824633; grid_points[ 343].w_total = 2.1192111689; grid_points[ 343].i_atom = 0;
    grid_points[ 344].x =  10.0265076477; grid_points[ 344].y = -30.0795229432; grid_points[ 344].z =  10.0265076477; grid_points[ 344].w_fixed = 3474.1698452873; grid_points[ 344].w_total = 0.0000000000; grid_points[ 344].i_atom = 0;
    grid_points[ 345].x =   9.6823189135; grid_points[ 345].y =  -9.6823189135; grid_points[ 345].z =   9.6823189135; grid_points[ 345].w_fixed =  396.8927232114; grid_points[ 345].w_total = 0.0000000000; grid_points[ 345].i_atom = 0;
    grid_points[ 346].x =  10.0265076477; grid_points[ 346].y = -10.0265076477; grid_points[ 346].z =  30.0795229432; grid_points[ 346].w_fixed = 3474.1698452873; grid_points[ 346].w_total = 0.0000000000; grid_points[ 346].i_atom = 0;
    grid_points[ 347].x =  11.8583704325; grid_points[ 347].y = -11.8583704325; grid_points[ 347].z =   0.0000000000; grid_points[ 347].w_fixed =  424.7625826047; grid_points[ 347].w_total = 0.0000000000; grid_points[ 347].i_atom = 0;
    grid_points[ 348].x =  11.8913062075; grid_points[ 348].y =  -2.2808122310; grid_points[ 348].z =   2.2808122310; grid_points[ 348].w_fixed =   57.9035858166; grid_points[ 348].w_total = 0.0000000000; grid_points[ 348].i_atom = 0;
    grid_points[ 349].x =  11.8583704325; grid_points[ 349].y =   0.0000000000; grid_points[ 349].z =  11.8583704325; grid_points[ 349].w_fixed =  424.7625826047; grid_points[ 349].w_total = 0.0000000000; grid_points[ 349].i_atom = 0;
    grid_points[ 350].x =  19.1993004361; grid_points[ 350].y = -19.1993004361; grid_points[ 350].z =  19.1993004361; grid_points[ 350].w_fixed = 3632.6798826243; grid_points[ 350].w_total = 0.0000000000; grid_points[ 350].i_atom = 0;
    grid_points[ 351].x =  21.0705920995; grid_points[ 351].y =   7.0235306998; grid_points[ 351].z =  -7.0235306998; grid_points[ 351].w_fixed = 1092.3998745831; grid_points[ 351].w_total = 0.0000000000; grid_points[ 351].i_atom = 0;
    grid_points[ 352].x =  30.0795229432; grid_points[ 352].y = -10.0265076477; grid_points[ 352].z =  10.0265076477; grid_points[ 352].w_fixed = 3474.1698452873; grid_points[ 352].w_total = 0.0000000000; grid_points[ 352].i_atom = 0;
    grid_points[ 353].x =  -8.0701458556; grid_points[ 353].y =   0.0000000000; grid_points[ 353].z =   4.3961380893; grid_points[ 353].w_fixed =   27.1448404564; grid_points[ 353].w_total = 0.0000023325; grid_points[ 353].i_atom = 0;
    grid_points[ 354].x =   0.0000000000; grid_points[ 354].y = -10.8198077400; grid_points[ 354].z =  -5.8939912332; grid_points[ 354].w_fixed =   68.3625199466; grid_points[ 354].w_total = 0.0000000000; grid_points[ 354].i_atom = 0;
    grid_points[ 355].x =   0.0000000000; grid_points[ 355].y = -10.8198077400; grid_points[ 355].z =   5.8939912332; grid_points[ 355].w_fixed =   68.3625199466; grid_points[ 355].w_total = 0.0000000000; grid_points[ 355].i_atom = 0;
    grid_points[ 356].x =   2.7423536150; grid_points[ 356].y =  -2.7423536150; grid_points[ 356].z =  -5.7438484187; grid_points[ 356].w_fixed =   11.1501458946; grid_points[ 356].w_total = 0.2040185179; grid_points[ 356].i_atom = 0;
    grid_points[ 357].x =   1.5919331643; grid_points[ 357].y =  -1.5919331643; grid_points[ 357].z =  -3.3342974947; grid_points[ 357].w_fixed =    2.1032476402; grid_points[ 357].w_total = 1.9796616798; grid_points[ 357].i_atom = 0;
    grid_points[ 358].x =   0.0000000000; grid_points[ 358].y =  -1.9245685660; grid_points[ 358].z =  -3.5329984457; grid_points[ 358].w_fixed =    2.1250626994; grid_points[ 358].w_total = 1.7554724438; grid_points[ 358].i_atom = 0;
    grid_points[ 359].x =   2.0842608672; grid_points[ 359].y =  -2.0842608672; grid_points[ 359].z =  -4.3654758527; grid_points[ 359].w_fixed =    4.7844044707; grid_points[ 359].w_total = 2.7579173487; grid_points[ 359].i_atom = 0;
    grid_points[ 360].x =   0.0000000000; grid_points[ 360].y =  -2.5197684414; grid_points[ 360].z =  -4.6256278652; grid_points[ 360].w_fixed =    4.8340287112; grid_points[ 360].w_total = 0.9672553950; grid_points[ 360].i_atom = 0;
    grid_points[ 361].x =   1.5919331643; grid_points[ 361].y =  -1.5919331643; grid_points[ 361].z =   3.3342974947; grid_points[ 361].w_fixed =    2.1032476402; grid_points[ 361].w_total = 1.2794689390; grid_points[ 361].i_atom = 0;
    grid_points[ 362].x =   2.0842608672; grid_points[ 362].y =  -2.0842608672; grid_points[ 362].z =   4.3654758527; grid_points[ 362].w_fixed =    4.7844044707; grid_points[ 362].w_total = 1.0736736243; grid_points[ 362].i_atom = 0;
    grid_points[ 363].x =   2.7423536150; grid_points[ 363].y =  -2.7423536150; grid_points[ 363].z =   5.7438484187; grid_points[ 363].w_fixed =   11.1501458946; grid_points[ 363].w_total = 0.0554063878; grid_points[ 363].i_atom = 0;
    grid_points[ 364].x =   0.7447550142; grid_points[ 364].y =  -0.7447550142; grid_points[ 364].z =  -3.8828754963; grid_points[ 364].w_fixed =    1.7999446257; grid_points[ 364].w_total = 0.6776187779; grid_points[ 364].i_atom = 0;
    grid_points[ 365].x =   0.9750809686; grid_points[ 365].y =  -0.9750809686; grid_points[ 365].z =  -5.0837093106; grid_points[ 365].w_fixed =    4.0944598961; grid_points[ 365].w_total = 0.0118697638; grid_points[ 365].i_atom = 0;
    grid_points[ 366].x =   0.7447550142; grid_points[ 366].y =  -0.7447550142; grid_points[ 366].z =   3.8828754963; grid_points[ 366].w_fixed =    1.7999446257; grid_points[ 366].w_total = 0.0335908066; grid_points[ 366].i_atom = 0;
    grid_points[ 367].x =   1.9245685660; grid_points[ 367].y =   0.0000000000; grid_points[ 367].z =   3.5329984457; grid_points[ 367].w_fixed =    2.1250626994; grid_points[ 367].w_total = 0.7185910872; grid_points[ 367].i_atom = 0;
    grid_points[ 368].x =   0.9750809686; grid_points[ 368].y =  -0.9750809686; grid_points[ 368].z =   5.0837093106; grid_points[ 368].w_fixed =    4.0944598961; grid_points[ 368].w_total = 0.0001490129; grid_points[ 368].i_atom = 0;
    grid_points[ 369].x =   2.5197684414; grid_points[ 369].y =   0.0000000000; grid_points[ 369].z =   4.6256278652; grid_points[ 369].w_fixed =    4.8340287112; grid_points[ 369].w_total = 0.2039773991; grid_points[ 369].i_atom = 0;
    grid_points[ 370].x =  10.8198077400; grid_points[ 370].y =   0.0000000000; grid_points[ 370].z =   5.8939912332; grid_points[ 370].w_fixed =   68.3625199466; grid_points[ 370].w_total = 0.0000000000; grid_points[ 370].i_atom = 0;
    grid_points[ 371].x = -16.4716395421; grid_points[ 371].y =   0.0000000000; grid_points[ 371].z =  16.4716395421; grid_points[ 371].w_fixed = 1222.4492551854; grid_points[ 371].w_total = 0.0000000000; grid_points[ 371].i_atom = 0;
    grid_points[ 372].x =  -7.0235306998; grid_points[ 372].y = -21.0705920995; grid_points[ 372].z =   7.0235306998; grid_points[ 372].w_fixed = 1092.3998745831; grid_points[ 372].w_total = 0.0000000000; grid_points[ 372].i_atom = 0;
    grid_points[ 373].x =  -6.6888597417; grid_points[ 373].y =  -1.2829568800; grid_points[ 373].z =  -1.2829568800; grid_points[ 373].w_fixed =    9.5422169009; grid_points[ 373].w_total = 0.0865328484; grid_points[ 373].i_atom = 0;
    grid_points[ 374].x =  -8.8693420267; grid_points[ 374].y =  -1.7011843293; grid_points[ 374].z =  -1.7011843293; grid_points[ 374].w_fixed =   22.9918908793; grid_points[ 374].w_total = 0.0000042871; grid_points[ 374].i_atom = 0;
    grid_points[ 375].x =  -6.3448636630; grid_points[ 375].y =  -1.9846144277; grid_points[ 375].z =   6.3448636630; grid_points[ 375].w_fixed =   27.8387029598; grid_points[ 375].w_total = 0.0001758536; grid_points[ 375].i_atom = 0;
    grid_points[ 376].x =  -8.5066870163; grid_points[ 376].y =  -2.6608126953; grid_points[ 376].z =   8.5066870163; grid_points[ 376].w_fixed =   70.1099676542; grid_points[ 376].w_total = 0.0000000000; grid_points[ 376].i_atom = 0;
    grid_points[ 377].x =  -6.3448636630; grid_points[ 377].y =   6.3448636630; grid_points[ 377].z =  -1.9846144277; grid_points[ 377].w_fixed =   27.8387029598; grid_points[ 377].w_total = 0.0013808480; grid_points[ 377].i_atom = 0;
    grid_points[ 378].x =  -7.0235306998; grid_points[ 378].y =   7.0235306998; grid_points[ 378].z = -21.0705920995; grid_points[ 378].w_fixed = 1092.3998745831; grid_points[ 378].w_total = 0.0000000000; grid_points[ 378].i_atom = 0;
    grid_points[ 379].x =  -8.5066870163; grid_points[ 379].y =   8.5066870163; grid_points[ 379].z =  -2.6608126953; grid_points[ 379].w_fixed =   70.1099676542; grid_points[ 379].w_total = 0.0000000000; grid_points[ 379].i_atom = 0;
    grid_points[ 380].x =  -7.1135404262; grid_points[ 380].y =   7.1135404262; grid_points[ 380].z =   7.1135404262; grid_points[ 380].w_fixed =   69.0587753800; grid_points[ 380].w_total = 0.0000000000; grid_points[ 380].i_atom = 0;
    grid_points[ 381].x =   0.0000000000; grid_points[ 381].y = -11.8583704325; grid_points[ 381].z = -11.8583704325; grid_points[ 381].w_fixed =  424.7625826047; grid_points[ 381].w_total = 0.0000000000; grid_points[ 381].i_atom = 0;
    grid_points[ 382].x =   2.2808122310; grid_points[ 382].y = -11.8913062075; grid_points[ 382].z =  -2.2808122310; grid_points[ 382].w_fixed =   57.9035858166; grid_points[ 382].w_total = 0.0000000000; grid_points[ 382].i_atom = 0;
    grid_points[ 383].x =   2.2808122310; grid_points[ 383].y =  -2.2808122310; grid_points[ 383].z = -11.8913062075; grid_points[ 383].w_fixed =   57.9035858166; grid_points[ 383].w_total = 0.0000000000; grid_points[ 383].i_atom = 0;
    grid_points[ 384].x =   0.0000000000; grid_points[ 384].y =  -1.2736799257; grid_points[ 384].z =  -1.2736799257; grid_points[ 384].w_fixed =    0.4501131044; grid_points[ 384].w_total = 0.4501105971; grid_points[ 384].i_atom = 0;
    grid_points[ 385].x =   0.0000000000; grid_points[ 385].y =  -1.6675833421; grid_points[ 385].z =  -1.6675833421; grid_points[ 385].w_fixed =    0.9966614143; grid_points[ 385].w_total = 0.9965883126; grid_points[ 385].i_atom = 0;
    grid_points[ 386].x =   1.3615760972; grid_points[ 386].y =  -1.3615760972; grid_points[ 386].z =  -1.3615760972; grid_points[ 386].w_fixed =    0.9312676753; grid_points[ 386].w_total = 0.9312379628; grid_points[ 386].i_atom = 0;
    grid_points[ 387].x =   1.7783851066; grid_points[ 387].y =  -1.7783851066; grid_points[ 387].z =  -1.7783851066; grid_points[ 387].w_fixed =    0.9591496581; grid_points[ 387].w_total = 0.9588589425; grid_points[ 387].i_atom = 0;
    grid_points[ 388].x =   0.6652031738; grid_points[ 388].y =  -2.1266717541; grid_points[ 388].z =  -2.1266717541; grid_points[ 388].w_fixed =    0.9737495508; grid_points[ 388].w_total = 0.9730669536; grid_points[ 388].i_atom = 0;
    grid_points[ 389].x =   0.0000000000; grid_points[ 389].y =  -1.4734978083; grid_points[ 389].z =  -2.7049519350; grid_points[ 389].w_fixed =    0.9494794437; grid_points[ 389].w_total = 0.9420711494; grid_points[ 389].i_atom = 0;
    grid_points[ 390].x =   0.0000000000; grid_points[ 390].y =  -2.7049519350; grid_points[ 390].z =  -1.4734978083; grid_points[ 390].w_fixed =    0.9494794437; grid_points[ 390].w_total = 0.9427515891; grid_points[ 390].i_atom = 0;
    grid_points[ 391].x =   2.3227887106; grid_points[ 391].y =  -2.3227887106; grid_points[ 391].z =  -2.3227887106; grid_points[ 391].w_fixed =    2.1467059397; grid_points[ 391].w_total = 2.1417068903; grid_points[ 391].i_atom = 0;
    grid_points[ 392].x =   0.8688367985; grid_points[ 392].y =  -2.7776937196; grid_points[ 392].z =  -2.7776937196; grid_points[ 392].w_fixed =    2.1793824633; grid_points[ 392].w_total = 2.1646215031; grid_points[ 392].i_atom = 0;
    grid_points[ 393].x =   0.5430989451; grid_points[ 393].y =  -1.6292968354; grid_points[ 393].z =  -0.5430989451; grid_points[ 393].w_fixed =    0.4022281471; grid_points[ 393].w_total = 0.4022099166; grid_points[ 393].i_atom = 0;
    grid_points[ 394].x =   0.7110599262; grid_points[ 394].y =  -2.1331797787; grid_points[ 394].z =  -0.7110599262; grid_points[ 394].w_fixed =    0.8906323100; grid_points[ 394].w_total = 0.8899740741; grid_points[ 394].i_atom = 0;
    grid_points[ 395].x =   0.5702030578; grid_points[ 395].y =  -2.9728265519; grid_points[ 395].z =  -0.5702030578; grid_points[ 395].w_fixed =    0.8042164697; grid_points[ 395].w_total = 0.7899201355; grid_points[ 395].i_atom = 0;
    grid_points[ 396].x =   2.1266717541; grid_points[ 396].y =  -2.1266717541; grid_points[ 396].z =  -0.6652031738; grid_points[ 396].w_fixed =    0.9737495508; grid_points[ 396].w_total = 0.9721447960; grid_points[ 396].i_atom = 0;
    grid_points[ 397].x =   1.2188238289; grid_points[ 397].y =  -2.5528215194; grid_points[ 397].z =  -1.2188238289; grid_points[ 397].w_fixed =    0.9397324606; grid_points[ 397].w_total = 0.9357421258; grid_points[ 397].i_atom = 0;
    grid_points[ 398].x =   1.2829568800; grid_points[ 398].y =  -1.2829568800; grid_points[ 398].z =   6.6888597417; grid_points[ 398].w_fixed =    9.5422169009; grid_points[ 398].w_total = 0.0000003249; grid_points[ 398].i_atom = 0;
    grid_points[ 399].x =   1.7011843293; grid_points[ 399].y =  -1.7011843293; grid_points[ 399].z =   8.8693420267; grid_points[ 399].w_fixed =   22.9918908793; grid_points[ 399].w_total = 0.0000000000; grid_points[ 399].i_atom = 0;
    grid_points[ 400].x =   0.5430989451; grid_points[ 400].y =  -0.5430989451; grid_points[ 400].z =  -1.6292968354; grid_points[ 400].w_fixed =    0.4022281471; grid_points[ 400].w_total = 0.4022275009; grid_points[ 400].i_atom = 0;
    grid_points[ 401].x =   0.7110599262; grid_points[ 401].y =  -0.7110599262; grid_points[ 401].z =  -2.1331797787; grid_points[ 401].w_fixed =    0.8906323100; grid_points[ 401].w_total = 0.8904412859; grid_points[ 401].i_atom = 0;
    grid_points[ 402].x =   0.5702030578; grid_points[ 402].y =  -0.5702030578; grid_points[ 402].z =  -2.9728265519; grid_points[ 402].w_fixed =    0.8042164697; grid_points[ 402].w_total = 0.7676826133; grid_points[ 402].i_atom = 0;
    grid_points[ 403].x =   2.1266717541; grid_points[ 403].y =  -0.6652031738; grid_points[ 403].z =  -2.1266717541; grid_points[ 403].w_fixed =    0.9737495508; grid_points[ 403].w_total = 0.9730669536; grid_points[ 403].i_atom = 0;
    grid_points[ 404].x =   1.2188238289; grid_points[ 404].y =  -1.2188238289; grid_points[ 404].z =  -2.5528215194; grid_points[ 404].w_fixed =    0.9397324606; grid_points[ 404].w_total = 0.9373660293; grid_points[ 404].i_atom = 0;
    grid_points[ 405].x =   0.0000000000; grid_points[ 405].y =  -0.4000533132; grid_points[ 405].z =  -0.4000533132; grid_points[ 405].w_fixed =    0.0275212788; grid_points[ 405].w_total = 0.0275212788; grid_points[ 405].i_atom = 0;
    grid_points[ 406].x =   0.3266421624; grid_points[ 406].y =  -0.3266421624; grid_points[ 406].z =  -0.3266421624; grid_points[ 406].w_fixed =    0.0232210790; grid_points[ 406].w_total = 0.0232210790; grid_points[ 406].i_atom = 0;
    grid_points[ 407].x =   0.0000000000; grid_points[ 407].y =  -0.5445170097; grid_points[ 407].z =  -0.5445170097; grid_points[ 407].w_fixed =    0.0655814715; grid_points[ 407].w_total = 0.0655814714; grid_points[ 407].i_atom = 0;
    grid_points[ 408].x =   0.4445962766; grid_points[ 408].y =  -0.4445962766; grid_points[ 408].z =  -0.4445962766; grid_points[ 408].w_fixed =    0.0553343666; grid_points[ 408].w_total = 0.0553343665; grid_points[ 408].i_atom = 0;
    grid_points[ 409].x =   0.0000000000; grid_points[ 409].y =  -0.7300449659; grid_points[ 409].z =  -0.7300449659; grid_points[ 409].w_fixed =    0.1512444148; grid_points[ 409].w_total = 0.1512444110; grid_points[ 409].i_atom = 0;
    grid_points[ 410].x =   0.5960792185; grid_points[ 410].y =  -0.5960792185; grid_points[ 410].z =  -0.5960792185; grid_points[ 410].w_fixed =    0.1276124749; grid_points[ 410].w_total = 0.1276124729; grid_points[ 410].i_atom = 0;
    grid_points[ 411].x =   0.0000000000; grid_points[ 411].y =  -0.9680302394; grid_points[ 411].z =  -0.9680302394; grid_points[ 411].w_fixed =    0.3411849149; grid_points[ 411].w_total = 0.3411847831; grid_points[ 411].i_atom = 0;
    grid_points[ 412].x =   0.7903933807; grid_points[ 412].y =  -0.7903933807; grid_points[ 412].z =  -0.7903933807; grid_points[ 412].w_fixed =    0.2878747720; grid_points[ 412].w_total = 0.2878747082; grid_points[ 412].i_atom = 0;
    grid_points[ 413].x =   1.0399553045; grid_points[ 413].y =  -1.0399553045; grid_points[ 413].z =  -1.0399553045; grid_points[ 413].w_fixed =    0.4205799265; grid_points[ 413].w_total = 0.4205787392; grid_points[ 413].i_atom = 0;
    grid_points[ 414].x =   1.6292968354; grid_points[ 414].y =  -0.5430989451; grid_points[ 414].z =  -0.5430989451; grid_points[ 414].w_fixed =    0.4022281471; grid_points[ 414].w_total = 0.4022099166; grid_points[ 414].i_atom = 0;
    grid_points[ 415].x =   2.1331797787; grid_points[ 415].y =  -0.7110599262; grid_points[ 415].z =  -0.7110599262; grid_points[ 415].w_fixed =    0.8906323100; grid_points[ 415].w_total = 0.8899740741; grid_points[ 415].i_atom = 0;
    grid_points[ 416].x =   2.5528215194; grid_points[ 416].y =  -1.2188238289; grid_points[ 416].z =  -1.2188238289; grid_points[ 416].w_fixed =    0.9397324606; grid_points[ 416].w_total = 0.9357421258; grid_points[ 416].i_atom = 0;
    grid_points[ 417].x =   1.2829568800; grid_points[ 417].y =   6.6888597417; grid_points[ 417].z =  -1.2829568800; grid_points[ 417].w_fixed =    9.5422169009; grid_points[ 417].w_total = 0.0865328484; grid_points[ 417].i_atom = 0;
    grid_points[ 418].x =   1.9846144277; grid_points[ 418].y =   6.3448636630; grid_points[ 418].z =   6.3448636630; grid_points[ 418].w_fixed =   27.8387029598; grid_points[ 418].w_total = 0.0001758536; grid_points[ 418].i_atom = 0;
    grid_points[ 419].x =   1.7011843293; grid_points[ 419].y =   8.8693420267; grid_points[ 419].z =  -1.7011843293; grid_points[ 419].w_fixed =   22.9918908793; grid_points[ 419].w_total = 0.0000042871; grid_points[ 419].i_atom = 0;
    grid_points[ 420].x =   2.6608126953; grid_points[ 420].y =   8.5066870163; grid_points[ 420].z =   8.5066870163; grid_points[ 420].w_fixed =   70.1099676542; grid_points[ 420].w_total = 0.0000000000; grid_points[ 420].i_atom = 0;
    grid_points[ 421].x =   3.0411437621; grid_points[ 421].y =  -3.0411437621; grid_points[ 421].z =  -3.0411437621; grid_points[ 421].w_fixed =    4.8832621034; grid_points[ 421].w_total = 4.8017896257; grid_points[ 421].i_atom = 0;
    grid_points[ 422].x =   2.7776937196; grid_points[ 422].y =  -2.7776937196; grid_points[ 422].z =  -0.8688367985; grid_points[ 422].w_fixed =    2.1793824633; grid_points[ 422].w_total = 2.1466862309; grid_points[ 422].i_atom = 0;
    grid_points[ 423].x =   2.7776937196; grid_points[ 423].y =  -0.8688367985; grid_points[ 423].z =  -2.7776937196; grid_points[ 423].w_fixed =    2.1793824633; grid_points[ 423].w_total = 2.1646215031; grid_points[ 423].i_atom = 0;
    grid_points[ 424].x =   2.9728265519; grid_points[ 424].y =  -0.5702030578; grid_points[ 424].z =  -0.5702030578; grid_points[ 424].w_fixed =    0.8042164697; grid_points[ 424].w_total = 0.7899201355; grid_points[ 424].i_atom = 0;
    grid_points[ 425].x =  10.0265076477; grid_points[ 425].y = -30.0795229432; grid_points[ 425].z = -10.0265076477; grid_points[ 425].w_fixed = 3474.1698452873; grid_points[ 425].w_total = 0.0000000000; grid_points[ 425].i_atom = 0;
    grid_points[ 426].x =  10.0265076477; grid_points[ 426].y = -10.0265076477; grid_points[ 426].z = -30.0795229432; grid_points[ 426].w_fixed = 3474.1698452873; grid_points[ 426].w_total = 0.0000000000; grid_points[ 426].i_atom = 0;
    grid_points[ 427].x =   9.6823189135; grid_points[ 427].y =  -9.6823189135; grid_points[ 427].z =  -9.6823189135; grid_points[ 427].w_fixed =  396.8927232114; grid_points[ 427].w_total = 0.0000000000; grid_points[ 427].i_atom = 0;
    grid_points[ 428].x =  11.8913062075; grid_points[ 428].y =  -2.2808122310; grid_points[ 428].z =  -2.2808122310; grid_points[ 428].w_fixed =   57.9035858166; grid_points[ 428].w_total = 0.0000000000; grid_points[ 428].i_atom = 0;
    grid_points[ 429].x =  19.1993004361; grid_points[ 429].y = -19.1993004361; grid_points[ 429].z = -19.1993004361; grid_points[ 429].w_fixed = 3632.6798826243; grid_points[ 429].w_total = 0.0000000000; grid_points[ 429].i_atom = 0;
    grid_points[ 430].x =  21.0705920995; grid_points[ 430].y =   7.0235306998; grid_points[ 430].z =   7.0235306998; grid_points[ 430].w_fixed = 1092.3998745831; grid_points[ 430].w_total = 0.0000000000; grid_points[ 430].i_atom = 0;
    grid_points[ 431].x =  30.0795229432; grid_points[ 431].y = -10.0265076477; grid_points[ 431].z = -10.0265076477; grid_points[ 431].w_fixed = 3474.1698452873; grid_points[ 431].w_total = 0.0000000000; grid_points[ 431].i_atom = 0;
    grid_points[ 432].x = -33.2541638252; grid_points[ 432].y =   0.0000000000; grid_points[ 432].z =   0.0000000000; grid_points[ 432].w_fixed = 2186.8690181113; grid_points[ 432].w_total = 0.0000000000; grid_points[ 432].i_atom = 0;
    grid_points[ 433].x = -23.2944160349; grid_points[ 433].y =   0.0000000000; grid_points[ 433].z =   0.0000000000; grid_points[ 433].w_fixed =  687.6277060418; grid_points[ 433].w_total = 0.0000000000; grid_points[ 433].i_atom = 0;
    grid_points[ 434].x =  -5.2674155091; grid_points[ 434].y =   0.0000000000; grid_points[ 434].z =   0.0000000000; grid_points[ 434].w_fixed =    1.9088165479; grid_points[ 434].w_total = 0.4997305992; grid_points[ 434].i_atom = 0;
    grid_points[ 435].x =  -5.7438484187; grid_points[ 435].y =   2.7423536150; grid_points[ 435].z =   2.7423536150; grid_points[ 435].w_fixed =   11.1501458946; grid_points[ 435].w_total = 0.7753792790; grid_points[ 435].i_atom = 0;
    grid_points[ 436].x =  -5.8939912332; grid_points[ 436].y =  10.8198077400; grid_points[ 436].z =   0.0000000000; grid_points[ 436].w_fixed =   68.3625199466; grid_points[ 436].w_total = 0.0000000000; grid_points[ 436].i_atom = 0;
    grid_points[ 437].x =  -4.3961380893; grid_points[ 437].y =  -8.0701458556; grid_points[ 437].z =   0.0000000000; grid_points[ 437].w_fixed =   27.1448404564; grid_points[ 437].w_total = 0.0003003544; grid_points[ 437].i_atom = 0;
    grid_points[ 438].x =  -4.3961380893; grid_points[ 438].y =   0.0000000000; grid_points[ 438].z =  -8.0701458556; grid_points[ 438].w_fixed =   27.1448404564; grid_points[ 438].w_total = 0.0000024398; grid_points[ 438].i_atom = 0;
    grid_points[ 439].x =  -4.0231880620; grid_points[ 439].y =   0.0000000000; grid_points[ 439].z =   0.0000000000; grid_points[ 439].w_fixed =    0.8391251042; grid_points[ 439].w_total = 0.6689355801; grid_points[ 439].i_atom = 0;
    grid_points[ 440].x =  -3.8828754963; grid_points[ 440].y =   0.7447550142; grid_points[ 440].z =   0.7447550142; grid_points[ 440].w_fixed =    1.7999446257; grid_points[ 440].w_total = 1.5204246946; grid_points[ 440].i_atom = 0;
    grid_points[ 441].x =  -3.3342974947; grid_points[ 441].y =   1.5919331643; grid_points[ 441].z =   1.5919331643; grid_points[ 441].w_fixed =    2.1032476402; grid_points[ 441].w_total = 2.0194520110; grid_points[ 441].i_atom = 0;
    grid_points[ 442].x =  -3.5329984457; grid_points[ 442].y =   1.9245685660; grid_points[ 442].z =   0.0000000000; grid_points[ 442].w_fixed =    2.1250626994; grid_points[ 442].w_total = 1.9718171965; grid_points[ 442].i_atom = 0;
    grid_points[ 443].x =  -5.0837093106; grid_points[ 443].y =   0.9750809686; grid_points[ 443].z =   0.9750809686; grid_points[ 443].w_fixed =    4.0944598961; grid_points[ 443].w_total = 1.3604510062; grid_points[ 443].i_atom = 0;
    grid_points[ 444].x =  -4.3654758527; grid_points[ 444].y =   2.0842608672; grid_points[ 444].z =   2.0842608672; grid_points[ 444].w_fixed =    4.7844044707; grid_points[ 444].w_total = 3.3058550217; grid_points[ 444].i_atom = 0;
    grid_points[ 445].x =  -4.6256278652; grid_points[ 445].y =   2.5197684414; grid_points[ 445].z =   0.0000000000; grid_points[ 445].w_fixed =    4.8340287112; grid_points[ 445].w_total = 2.6697241981; grid_points[ 445].i_atom = 0;
    grid_points[ 446].x =   4.3961380893; grid_points[ 446].y =  -8.0701458556; grid_points[ 446].z =   0.0000000000; grid_points[ 446].w_fixed =   27.1448404564; grid_points[ 446].w_total = 0.0003003544; grid_points[ 446].i_atom = 0;
    grid_points[ 447].x =   4.3961380893; grid_points[ 447].y =   0.0000000000; grid_points[ 447].z =  -8.0701458556; grid_points[ 447].w_fixed =   27.1448404564; grid_points[ 447].w_total = 0.0000024398; grid_points[ 447].i_atom = 0;
    grid_points[ 448].x =   4.0231880620; grid_points[ 448].y =   0.0000000000; grid_points[ 448].z =   0.0000000000; grid_points[ 448].w_fixed =    0.8391251042; grid_points[ 448].w_total = 0.6689355801; grid_points[ 448].i_atom = 0;
    grid_points[ 449].x =   3.8828754963; grid_points[ 449].y =   0.7447550142; grid_points[ 449].z =   0.7447550142; grid_points[ 449].w_fixed =    1.7999446257; grid_points[ 449].w_total = 1.5204246946; grid_points[ 449].i_atom = 0;
    grid_points[ 450].x =   3.3342974947; grid_points[ 450].y =   1.5919331643; grid_points[ 450].z =   1.5919331643; grid_points[ 450].w_fixed =    2.1032476402; grid_points[ 450].w_total = 2.0194520110; grid_points[ 450].i_atom = 0;
    grid_points[ 451].x =   3.5329984457; grid_points[ 451].y =   1.9245685660; grid_points[ 451].z =   0.0000000000; grid_points[ 451].w_fixed =    2.1250626994; grid_points[ 451].w_total = 1.9718171965; grid_points[ 451].i_atom = 0;
    grid_points[ 452].x =   5.2674155091; grid_points[ 452].y =   0.0000000000; grid_points[ 452].z =   0.0000000000; grid_points[ 452].w_fixed =    1.9088165479; grid_points[ 452].w_total = 0.4997305992; grid_points[ 452].i_atom = 0;
    grid_points[ 453].x =   5.0837093106; grid_points[ 453].y =   0.9750809686; grid_points[ 453].z =   0.9750809686; grid_points[ 453].w_fixed =    4.0944598961; grid_points[ 453].w_total = 1.3604510062; grid_points[ 453].i_atom = 0;
    grid_points[ 454].x =   4.3654758527; grid_points[ 454].y =   2.0842608672; grid_points[ 454].z =   2.0842608672; grid_points[ 454].w_fixed =    4.7844044707; grid_points[ 454].w_total = 3.3058550217; grid_points[ 454].i_atom = 0;
    grid_points[ 455].x =   4.6256278652; grid_points[ 455].y =   2.5197684414; grid_points[ 455].z =   0.0000000000; grid_points[ 455].w_fixed =    4.8340287112; grid_points[ 455].w_total = 2.6697241981; grid_points[ 455].i_atom = 0;
    grid_points[ 456].x =   5.7438484187; grid_points[ 456].y =   2.7423536150; grid_points[ 456].z =   2.7423536150; grid_points[ 456].w_fixed =   11.1501458946; grid_points[ 456].w_total = 0.7753792790; grid_points[ 456].i_atom = 0;
    grid_points[ 457].x =   5.8939912332; grid_points[ 457].y =  10.8198077400; grid_points[ 457].z =   0.0000000000; grid_points[ 457].w_fixed =   68.3625199466; grid_points[ 457].w_total = 0.0000000000; grid_points[ 457].i_atom = 0;
    grid_points[ 458].x =  23.2944160349; grid_points[ 458].y =   0.0000000000; grid_points[ 458].z =   0.0000000000; grid_points[ 458].w_fixed =  687.6277060418; grid_points[ 458].w_total = 0.0000000000; grid_points[ 458].i_atom = 0;
    grid_points[ 459].x =  33.2541638252; grid_points[ 459].y =   0.0000000000; grid_points[ 459].z =   0.0000000000; grid_points[ 459].w_fixed = 2186.8690181113; grid_points[ 459].w_total = 0.0000000000; grid_points[ 459].i_atom = 0;
    grid_points[ 460].x = -23.5142447435; grid_points[ 460].y =   0.0000000000; grid_points[ 460].z = -23.5142447435; grid_points[ 460].w_fixed = 3887.7671433090; grid_points[ 460].w_total = 0.0000000000; grid_points[ 460].i_atom = 0;
    grid_points[ 461].x =  -6.0861418537; grid_points[ 461].y =   0.0000000000; grid_points[ 461].z =  -3.3153700687; grid_points[ 461].w_fixed =   11.2657961338; grid_points[ 461].w_total = 0.3549590929; grid_points[ 461].i_atom = 0;
    grid_points[ 462].x =  -3.6363256383; grid_points[ 462].y =  -7.6162691612; grid_points[ 462].z =  -3.6363256383; grid_points[ 462].w_fixed =   26.8661821836; grid_points[ 462].w_total = 0.0010773582; grid_points[ 462].i_atom = 0;
    grid_points[ 463].x =  -3.6363256383; grid_points[ 463].y =  -7.6162691612; grid_points[ 463].z =   3.6363256383; grid_points[ 463].w_fixed =   26.8661821836; grid_points[ 463].w_total = 0.0009330210; grid_points[ 463].i_atom = 0;
    grid_points[ 464].x =  -3.3153700687; grid_points[ 464].y =   0.0000000000; grid_points[ 464].z =  -6.0861418537; grid_points[ 464].w_fixed =   11.2657961338; grid_points[ 464].w_total = 0.0242256985; grid_points[ 464].i_atom = 0;
    grid_points[ 465].x =  -3.6367345381; grid_points[ 465].y =   1.1375367884; grid_points[ 465].z =  -3.6367345381; grid_points[ 465].w_fixed =    4.9575936763; grid_points[ 465].w_total = 4.6559094346; grid_points[ 465].i_atom = 0;
    grid_points[ 466].x =  -4.7850114467; grid_points[ 466].y =   1.4967071411; grid_points[ 466].z =  -4.7850114467; grid_points[ 466].w_fixed =   11.5537666422; grid_points[ 466].w_total = 2.9255841962; grid_points[ 466].i_atom = 0;
    grid_points[ 467].x =  -3.6367345381; grid_points[ 467].y =   1.1375367884; grid_points[ 467].z =   3.6367345381; grid_points[ 467].w_fixed =    4.9575936763; grid_points[ 467].w_total = 4.5221022069; grid_points[ 467].i_atom = 0;
    grid_points[ 468].x =  -4.7850114467; grid_points[ 468].y =   1.4967071411; grid_points[ 468].z =   4.7850114467; grid_points[ 468].w_fixed =   11.5537666422; grid_points[ 468].w_total = 1.8853671084; grid_points[ 468].i_atom = 0;
    grid_points[ 469].x =  -4.8752953156; grid_points[ 469].y =  10.2112860776; grid_points[ 469].z =  -4.8752953156; grid_points[ 469].w_fixed =   67.6607371617; grid_points[ 469].w_total = 0.0000000000; grid_points[ 469].i_atom = 0;
    grid_points[ 470].x =  -4.8752953156; grid_points[ 470].y =  10.2112860776; grid_points[ 470].z =   4.8752953156; grid_points[ 470].w_fixed =   67.6607371617; grid_points[ 470].w_total = 0.0000000000; grid_points[ 470].i_atom = 0;
    grid_points[ 471].x =   3.6363256383; grid_points[ 471].y =  -7.6162691612; grid_points[ 471].z =  -3.6363256383; grid_points[ 471].w_fixed =   26.8661821836; grid_points[ 471].w_total = 0.0010773582; grid_points[ 471].i_atom = 0;
    grid_points[ 472].x =   3.6363256383; grid_points[ 472].y =  -7.6162691612; grid_points[ 472].z =   3.6363256383; grid_points[ 472].w_fixed =   26.8661821836; grid_points[ 472].w_total = 0.0009330210; grid_points[ 472].i_atom = 0;
    grid_points[ 473].x =   3.3153700687; grid_points[ 473].y =   0.0000000000; grid_points[ 473].z =  -6.0861418537; grid_points[ 473].w_fixed =   11.2657961338; grid_points[ 473].w_total = 0.0242256985; grid_points[ 473].i_atom = 0;
    grid_points[ 474].x =   3.6367345381; grid_points[ 474].y =   1.1375367884; grid_points[ 474].z =  -3.6367345381; grid_points[ 474].w_fixed =    4.9575936763; grid_points[ 474].w_total = 4.6559094346; grid_points[ 474].i_atom = 0;
    grid_points[ 475].x =   4.7850114467; grid_points[ 475].y =   1.4967071411; grid_points[ 475].z =  -4.7850114467; grid_points[ 475].w_fixed =   11.5537666422; grid_points[ 475].w_total = 2.9255841962; grid_points[ 475].i_atom = 0;
    grid_points[ 476].x =   6.0861418537; grid_points[ 476].y =   0.0000000000; grid_points[ 476].z =  -3.3153700687; grid_points[ 476].w_fixed =   11.2657961338; grid_points[ 476].w_total = 0.3549590929; grid_points[ 476].i_atom = 0;
    grid_points[ 477].x =   3.6367345381; grid_points[ 477].y =   1.1375367884; grid_points[ 477].z =   3.6367345381; grid_points[ 477].w_fixed =    4.9575936763; grid_points[ 477].w_total = 4.5221022069; grid_points[ 477].i_atom = 0;
    grid_points[ 478].x =   4.7850114467; grid_points[ 478].y =   1.4967071411; grid_points[ 478].z =   4.7850114467; grid_points[ 478].w_fixed =   11.5537666422; grid_points[ 478].w_total = 1.8853671084; grid_points[ 478].i_atom = 0;
    grid_points[ 479].x =   4.8752953156; grid_points[ 479].y =  10.2112860776; grid_points[ 479].z =  -4.8752953156; grid_points[ 479].w_fixed =   67.6607371617; grid_points[ 479].w_total = 0.0000000000; grid_points[ 479].i_atom = 0;
    grid_points[ 480].x =   4.8752953156; grid_points[ 480].y =  10.2112860776; grid_points[ 480].z =   4.8752953156; grid_points[ 480].w_fixed =   67.6607371617; grid_points[ 480].w_total = 0.0000000000; grid_points[ 480].i_atom = 0;
    grid_points[ 481].x =  23.5142447435; grid_points[ 481].y =   0.0000000000; grid_points[ 481].z = -23.5142447435; grid_points[ 481].w_fixed = 3887.7671433090; grid_points[ 481].w_total = 0.0000000000; grid_points[ 481].i_atom = 0;
    grid_points[ 482].x =  -5.8939912332; grid_points[ 482].y =   0.0000000000; grid_points[ 482].z = -10.8198077400; grid_points[ 482].w_fixed =   68.3625199466; grid_points[ 482].w_total = 0.0000000000; grid_points[ 482].i_atom = 0;
    grid_points[ 483].x =  -5.7438484187; grid_points[ 483].y =   2.7423536150; grid_points[ 483].z =  -2.7423536150; grid_points[ 483].w_fixed =   11.1501458946; grid_points[ 483].w_total = 0.7975816561; grid_points[ 483].i_atom = 0;
    grid_points[ 484].x =  -3.3342974947; grid_points[ 484].y =   1.5919331643; grid_points[ 484].z =  -1.5919331643; grid_points[ 484].w_fixed =    2.1032476402; grid_points[ 484].w_total = 2.0149126278; grid_points[ 484].i_atom = 0;
    grid_points[ 485].x =  -3.5329984457; grid_points[ 485].y =   0.0000000000; grid_points[ 485].z =  -1.9245685660; grid_points[ 485].w_fixed =    2.1250626994; grid_points[ 485].w_total = 1.9769819275; grid_points[ 485].i_atom = 0;
    grid_points[ 486].x =  -4.3654758527; grid_points[ 486].y =   2.0842608672; grid_points[ 486].z =  -2.0842608672; grid_points[ 486].w_fixed =    4.7844044707; grid_points[ 486].w_total = 3.2966811514; grid_points[ 486].i_atom = 0;
    grid_points[ 487].x =  -4.6256278652; grid_points[ 487].y =   0.0000000000; grid_points[ 487].z =  -2.5197684414; grid_points[ 487].w_fixed =    4.8340287112; grid_points[ 487].w_total = 2.6774240028; grid_points[ 487].i_atom = 0;
    grid_points[ 488].x =  -3.8828754963; grid_points[ 488].y =   0.7447550142; grid_points[ 488].z =  -0.7447550142; grid_points[ 488].w_fixed =    1.7999446257; grid_points[ 488].w_total = 1.5182141465; grid_points[ 488].i_atom = 0;
    grid_points[ 489].x =  -5.0837093106; grid_points[ 489].y =   0.9750809686; grid_points[ 489].z =  -0.9750809686; grid_points[ 489].w_fixed =    4.0944598961; grid_points[ 489].w_total = 1.3639957134; grid_points[ 489].i_atom = 0;
    grid_points[ 490].x =   5.8939912332; grid_points[ 490].y =   0.0000000000; grid_points[ 490].z = -10.8198077400; grid_points[ 490].w_fixed =   68.3625199466; grid_points[ 490].w_total = 0.0000000000; grid_points[ 490].i_atom = 0;
    grid_points[ 491].x =   3.3342974947; grid_points[ 491].y =   1.5919331643; grid_points[ 491].z =  -1.5919331643; grid_points[ 491].w_fixed =    2.1032476402; grid_points[ 491].w_total = 2.0149126278; grid_points[ 491].i_atom = 0;
    grid_points[ 492].x =   3.5329984457; grid_points[ 492].y =   0.0000000000; grid_points[ 492].z =  -1.9245685660; grid_points[ 492].w_fixed =    2.1250626994; grid_points[ 492].w_total = 1.9769819275; grid_points[ 492].i_atom = 0;
    grid_points[ 493].x =   4.3654758527; grid_points[ 493].y =   2.0842608672; grid_points[ 493].z =  -2.0842608672; grid_points[ 493].w_fixed =    4.7844044707; grid_points[ 493].w_total = 3.2966811514; grid_points[ 493].i_atom = 0;
    grid_points[ 494].x =   4.6256278652; grid_points[ 494].y =   0.0000000000; grid_points[ 494].z =  -2.5197684414; grid_points[ 494].w_fixed =    4.8340287112; grid_points[ 494].w_total = 2.6774240028; grid_points[ 494].i_atom = 0;
    grid_points[ 495].x =   5.7438484187; grid_points[ 495].y =   2.7423536150; grid_points[ 495].z =  -2.7423536150; grid_points[ 495].w_fixed =   11.1501458946; grid_points[ 495].w_total = 0.7975816561; grid_points[ 495].i_atom = 0;
    grid_points[ 496].x =   3.8828754963; grid_points[ 496].y =   0.7447550142; grid_points[ 496].z =  -0.7447550142; grid_points[ 496].w_fixed =    1.7999446257; grid_points[ 496].w_total = 1.5182141465; grid_points[ 496].i_atom = 0;
    grid_points[ 497].x =   5.0837093106; grid_points[ 497].y =   0.9750809686; grid_points[ 497].z =  -0.9750809686; grid_points[ 497].w_fixed =    4.0944598961; grid_points[ 497].w_total = 1.3639957134; grid_points[ 497].i_atom = 0;
    grid_points[ 498].x = -23.5142447435; grid_points[ 498].y = -23.5142447435; grid_points[ 498].z =   0.0000000000; grid_points[ 498].w_fixed = 3887.7671433090; grid_points[ 498].w_total = 0.0000000000; grid_points[ 498].i_atom = 0;
    grid_points[ 499].x = -23.5142447435; grid_points[ 499].y =  23.5142447435; grid_points[ 499].z =   0.0000000000; grid_points[ 499].w_fixed = 3887.7671433090; grid_points[ 499].w_total = 0.0000000000; grid_points[ 499].i_atom = 0;
    grid_points[ 500].x =  -6.0861418537; grid_points[ 500].y =  -3.3153700687; grid_points[ 500].z =   0.0000000000; grid_points[ 500].w_fixed =   11.2657961338; grid_points[ 500].w_total = 0.4090867032; grid_points[ 500].i_atom = 0;
    grid_points[ 501].x =  -6.0861418537; grid_points[ 501].y =   3.3153700687; grid_points[ 501].z =   0.0000000000; grid_points[ 501].w_fixed =   11.2657961338; grid_points[ 501].w_total = 0.4090867032; grid_points[ 501].i_atom = 0;
    grid_points[ 502].x =  -3.3153700687; grid_points[ 502].y =  -6.0861418537; grid_points[ 502].z =   0.0000000000; grid_points[ 502].w_fixed =   11.2657961338; grid_points[ 502].w_total = 0.4090867032; grid_points[ 502].i_atom = 0;
    grid_points[ 503].x =  -3.6363256383; grid_points[ 503].y =  -3.6363256383; grid_points[ 503].z =  -7.6162691612; grid_points[ 503].w_fixed =   26.8661821836; grid_points[ 503].w_total = 0.0000817454; grid_points[ 503].i_atom = 0;
    grid_points[ 504].x =  -3.6367345381; grid_points[ 504].y =  -3.6367345381; grid_points[ 504].z =   1.1375367884; grid_points[ 504].w_fixed =    4.9575936763; grid_points[ 504].w_total = 4.3797383969; grid_points[ 504].i_atom = 0;
    grid_points[ 505].x =  -4.7850114467; grid_points[ 505].y =  -4.7850114467; grid_points[ 505].z =   1.4967071411; grid_points[ 505].w_fixed =   11.5537666422; grid_points[ 505].w_total = 2.4842047330; grid_points[ 505].i_atom = 0;
    grid_points[ 506].x =  -4.8752953156; grid_points[ 506].y =  -4.8752953156; grid_points[ 506].z =  10.2112860776; grid_points[ 506].w_fixed =   67.6607371617; grid_points[ 506].w_total = 0.0000000000; grid_points[ 506].i_atom = 0;
    grid_points[ 507].x =  -3.6363256383; grid_points[ 507].y =   3.6363256383; grid_points[ 507].z =  -7.6162691612; grid_points[ 507].w_fixed =   26.8661821836; grid_points[ 507].w_total = 0.0000817454; grid_points[ 507].i_atom = 0;
    grid_points[ 508].x =  -3.6367345381; grid_points[ 508].y =   3.6367345381; grid_points[ 508].z =   1.1375367884; grid_points[ 508].w_fixed =    4.9575936763; grid_points[ 508].w_total = 4.3797383969; grid_points[ 508].i_atom = 0;
    grid_points[ 509].x =  -4.7850114467; grid_points[ 509].y =   4.7850114467; grid_points[ 509].z =   1.4967071411; grid_points[ 509].w_fixed =   11.5537666422; grid_points[ 509].w_total = 2.4842047330; grid_points[ 509].i_atom = 0;
    grid_points[ 510].x =  -3.3153700687; grid_points[ 510].y =   6.0861418537; grid_points[ 510].z =   0.0000000000; grid_points[ 510].w_fixed =   11.2657961338; grid_points[ 510].w_total = 0.4090867032; grid_points[ 510].i_atom = 0;
    grid_points[ 511].x =  -4.8752953156; grid_points[ 511].y =   4.8752953156; grid_points[ 511].z =  10.2112860776; grid_points[ 511].w_fixed =   67.6607371617; grid_points[ 511].w_total = 0.0000000000; grid_points[ 511].i_atom = 0;
    grid_points[ 512].x =   3.3153700687; grid_points[ 512].y =  -6.0861418537; grid_points[ 512].z =   0.0000000000; grid_points[ 512].w_fixed =   11.2657961338; grid_points[ 512].w_total = 0.4090867032; grid_points[ 512].i_atom = 0;
    grid_points[ 513].x =   3.6363256383; grid_points[ 513].y =  -3.6363256383; grid_points[ 513].z =  -7.6162691612; grid_points[ 513].w_fixed =   26.8661821836; grid_points[ 513].w_total = 0.0000817454; grid_points[ 513].i_atom = 0;
    grid_points[ 514].x =   3.6367345381; grid_points[ 514].y =  -3.6367345381; grid_points[ 514].z =   1.1375367884; grid_points[ 514].w_fixed =    4.9575936763; grid_points[ 514].w_total = 4.3797383969; grid_points[ 514].i_atom = 0;
    grid_points[ 515].x =   4.7850114467; grid_points[ 515].y =  -4.7850114467; grid_points[ 515].z =   1.4967071411; grid_points[ 515].w_fixed =   11.5537666422; grid_points[ 515].w_total = 2.4842047330; grid_points[ 515].i_atom = 0;
    grid_points[ 516].x =   6.0861418537; grid_points[ 516].y =  -3.3153700687; grid_points[ 516].z =   0.0000000000; grid_points[ 516].w_fixed =   11.2657961338; grid_points[ 516].w_total = 0.4090867032; grid_points[ 516].i_atom = 0;
    grid_points[ 517].x =   4.8752953156; grid_points[ 517].y =  -4.8752953156; grid_points[ 517].z =  10.2112860776; grid_points[ 517].w_fixed =   67.6607371617; grid_points[ 517].w_total = 0.0000000000; grid_points[ 517].i_atom = 0;
    grid_points[ 518].x =   3.6363256383; grid_points[ 518].y =   3.6363256383; grid_points[ 518].z =  -7.6162691612; grid_points[ 518].w_fixed =   26.8661821836; grid_points[ 518].w_total = 0.0000817454; grid_points[ 518].i_atom = 0;
    grid_points[ 519].x =   3.6367345381; grid_points[ 519].y =   3.6367345381; grid_points[ 519].z =   1.1375367884; grid_points[ 519].w_fixed =    4.9575936763; grid_points[ 519].w_total = 4.3797383969; grid_points[ 519].i_atom = 0;
    grid_points[ 520].x =   4.7850114467; grid_points[ 520].y =   4.7850114467; grid_points[ 520].z =   1.4967071411; grid_points[ 520].w_fixed =   11.5537666422; grid_points[ 520].w_total = 2.4842047330; grid_points[ 520].i_atom = 0;
    grid_points[ 521].x =   3.3153700687; grid_points[ 521].y =   6.0861418537; grid_points[ 521].z =   0.0000000000; grid_points[ 521].w_fixed =   11.2657961338; grid_points[ 521].w_total = 0.4090867032; grid_points[ 521].i_atom = 0;
    grid_points[ 522].x =   6.0861418537; grid_points[ 522].y =   3.3153700687; grid_points[ 522].z =   0.0000000000; grid_points[ 522].w_fixed =   11.2657961338; grid_points[ 522].w_total = 0.4090867032; grid_points[ 522].i_atom = 0;
    grid_points[ 523].x =   4.8752953156; grid_points[ 523].y =   4.8752953156; grid_points[ 523].z =  10.2112860776; grid_points[ 523].w_fixed =   67.6607371617; grid_points[ 523].w_total = 0.0000000000; grid_points[ 523].i_atom = 0;
    grid_points[ 524].x =  23.5142447435; grid_points[ 524].y = -23.5142447435; grid_points[ 524].z =   0.0000000000; grid_points[ 524].w_fixed = 3887.7671433090; grid_points[ 524].w_total = 0.0000000000; grid_points[ 524].i_atom = 0;
    grid_points[ 525].x =  23.5142447435; grid_points[ 525].y =  23.5142447435; grid_points[ 525].z =   0.0000000000; grid_points[ 525].w_fixed = 3887.7671433090; grid_points[ 525].w_total = 0.0000000000; grid_points[ 525].i_atom = 0;
    grid_points[ 526].x = -13.4490373684; grid_points[ 526].y = -13.4490373684; grid_points[ 526].z = -13.4490373684; grid_points[ 526].w_fixed = 1142.2409452901; grid_points[ 526].w_total = 0.0000000000; grid_points[ 526].i_atom = 0;
    grid_points[ 527].x = -13.4490373684; grid_points[ 527].y = -13.4490373684; grid_points[ 527].z =  13.4490373684; grid_points[ 527].w_fixed = 1142.2409452901; grid_points[ 527].w_total = 0.0000000000; grid_points[ 527].i_atom = 0;
    grid_points[ 528].x = -15.1692784261; grid_points[ 528].y =  -5.0564261420; grid_points[ 528].z =  -5.0564261420; grid_points[ 528].w_fixed =  379.5745222116; grid_points[ 528].w_total = 0.0000000000; grid_points[ 528].i_atom = 0;
    grid_points[ 529].x = -15.1692784261; grid_points[ 529].y =  -5.0564261420; grid_points[ 529].z =   5.0564261420; grid_points[ 529].w_fixed =  379.5745222116; grid_points[ 529].w_total = 0.0000000000; grid_points[ 529].i_atom = 0;
    grid_points[ 530].x = -15.1692784261; grid_points[ 530].y =   5.0564261420; grid_points[ 530].z =  -5.0564261420; grid_points[ 530].w_fixed =  379.5745222116; grid_points[ 530].w_total = 0.0000000000; grid_points[ 530].i_atom = 0;
    grid_points[ 531].x = -15.1692784261; grid_points[ 531].y =   5.0564261420; grid_points[ 531].z =   5.0564261420; grid_points[ 531].w_fixed =  379.5745222116; grid_points[ 531].w_total = 0.0000000000; grid_points[ 531].i_atom = 0;
    grid_points[ 532].x = -13.4490373684; grid_points[ 532].y =  13.4490373684; grid_points[ 532].z = -13.4490373684; grid_points[ 532].w_fixed = 1142.2409452901; grid_points[ 532].w_total = 0.0000000000; grid_points[ 532].i_atom = 0;
    grid_points[ 533].x = -13.4490373684; grid_points[ 533].y =  13.4490373684; grid_points[ 533].z =  13.4490373684; grid_points[ 533].w_fixed = 1142.2409452901; grid_points[ 533].w_total = 0.0000000000; grid_points[ 533].i_atom = 0;
    grid_points[ 534].x =  -5.3057605245; grid_points[ 534].y =  -5.3057605245; grid_points[ 534].z =  -5.3057605245; grid_points[ 534].w_fixed =   27.4213039728; grid_points[ 534].w_total = 0.0360281880; grid_points[ 534].i_atom = 0;
    grid_points[ 535].x =  -5.3057605245; grid_points[ 535].y =  -5.3057605245; grid_points[ 535].z =   5.3057605245; grid_points[ 535].w_fixed =   27.4213039728; grid_points[ 535].w_total = 0.0257436288; grid_points[ 535].i_atom = 0;
    grid_points[ 536].x =  -5.3057605245; grid_points[ 536].y =   5.3057605245; grid_points[ 536].z =  -5.3057605245; grid_points[ 536].w_fixed =   27.4213039728; grid_points[ 536].w_total = 0.0360281880; grid_points[ 536].i_atom = 0;
    grid_points[ 537].x =  -5.3057605245; grid_points[ 537].y =   5.3057605245; grid_points[ 537].z =   5.3057605245; grid_points[ 537].w_fixed =   27.4213039728; grid_points[ 537].w_total = 0.0257436288; grid_points[ 537].i_atom = 0;
    grid_points[ 538].x =  -5.0564261420; grid_points[ 538].y = -15.1692784261; grid_points[ 538].z =  -5.0564261420; grid_points[ 538].w_fixed =  379.5745222116; grid_points[ 538].w_total = 0.0000000000; grid_points[ 538].i_atom = 0;
    grid_points[ 539].x =  -5.0564261420; grid_points[ 539].y = -15.1692784261; grid_points[ 539].z =   5.0564261420; grid_points[ 539].w_fixed =  379.5745222116; grid_points[ 539].w_total = 0.0000000000; grid_points[ 539].i_atom = 0;
    grid_points[ 540].x =  -5.0564261420; grid_points[ 540].y =  -5.0564261420; grid_points[ 540].z = -15.1692784261; grid_points[ 540].w_fixed =  379.5745222116; grid_points[ 540].w_total = 0.0000000000; grid_points[ 540].i_atom = 0;
    grid_points[ 541].x =  -4.0013664898; grid_points[ 541].y =  -4.0013664898; grid_points[ 541].z =  -4.0013664898; grid_points[ 541].w_fixed =   11.3805354935; grid_points[ 541].w_total = 9.1115392109; grid_points[ 541].i_atom = 0;
    grid_points[ 542].x =  -4.0013664898; grid_points[ 542].y =  -4.0013664898; grid_points[ 542].z =   4.0013664898; grid_points[ 542].w_fixed =   11.3805354935; grid_points[ 542].w_total = 8.8295643245; grid_points[ 542].i_atom = 0;
    grid_points[ 543].x =  -5.0564261420; grid_points[ 543].y =  -5.0564261420; grid_points[ 543].z =  15.1692784261; grid_points[ 543].w_fixed =  379.5745222116; grid_points[ 543].w_total = 0.0000000000; grid_points[ 543].i_atom = 0;
    grid_points[ 544].x =  -5.0564261420; grid_points[ 544].y =   5.0564261420; grid_points[ 544].z = -15.1692784261; grid_points[ 544].w_fixed =  379.5745222116; grid_points[ 544].w_total = 0.0000000000; grid_points[ 544].i_atom = 0;
    grid_points[ 545].x =  -4.0013664898; grid_points[ 545].y =   4.0013664898; grid_points[ 545].z =  -4.0013664898; grid_points[ 545].w_fixed =   11.3805354935; grid_points[ 545].w_total = 9.1115392109; grid_points[ 545].i_atom = 0;
    grid_points[ 546].x =  -4.0013664898; grid_points[ 546].y =   4.0013664898; grid_points[ 546].z =   4.0013664898; grid_points[ 546].w_fixed =   11.3805354935; grid_points[ 546].w_total = 8.8295643245; grid_points[ 546].i_atom = 0;
    grid_points[ 547].x =  -5.0564261420; grid_points[ 547].y =   5.0564261420; grid_points[ 547].z =  15.1692784261; grid_points[ 547].w_fixed =  379.5745222116; grid_points[ 547].w_total = 0.0000000000; grid_points[ 547].i_atom = 0;
    grid_points[ 548].x =  -5.0564261420; grid_points[ 548].y =  15.1692784261; grid_points[ 548].z =  -5.0564261420; grid_points[ 548].w_fixed =  379.5745222116; grid_points[ 548].w_total = 0.0000000000; grid_points[ 548].i_atom = 0;
    grid_points[ 549].x =  -5.0564261420; grid_points[ 549].y =  15.1692784261; grid_points[ 549].z =   5.0564261420; grid_points[ 549].w_fixed =  379.5745222116; grid_points[ 549].w_total = 0.0000000000; grid_points[ 549].i_atom = 0;
    grid_points[ 550].x =   5.0564261420; grid_points[ 550].y = -15.1692784261; grid_points[ 550].z =  -5.0564261420; grid_points[ 550].w_fixed =  379.5745222116; grid_points[ 550].w_total = 0.0000000000; grid_points[ 550].i_atom = 0;
    grid_points[ 551].x =   5.0564261420; grid_points[ 551].y = -15.1692784261; grid_points[ 551].z =   5.0564261420; grid_points[ 551].w_fixed =  379.5745222116; grid_points[ 551].w_total = 0.0000000000; grid_points[ 551].i_atom = 0;
    grid_points[ 552].x =   5.3057605245; grid_points[ 552].y =  -5.3057605245; grid_points[ 552].z =  -5.3057605245; grid_points[ 552].w_fixed =   27.4213039728; grid_points[ 552].w_total = 0.0360281880; grid_points[ 552].i_atom = 0;
    grid_points[ 553].x =   5.3057605245; grid_points[ 553].y =  -5.3057605245; grid_points[ 553].z =   5.3057605245; grid_points[ 553].w_fixed =   27.4213039728; grid_points[ 553].w_total = 0.0257436288; grid_points[ 553].i_atom = 0;
    grid_points[ 554].x =   5.0564261420; grid_points[ 554].y =  -5.0564261420; grid_points[ 554].z = -15.1692784261; grid_points[ 554].w_fixed =  379.5745222116; grid_points[ 554].w_total = 0.0000000000; grid_points[ 554].i_atom = 0;
    grid_points[ 555].x =   4.0013664898; grid_points[ 555].y =  -4.0013664898; grid_points[ 555].z =  -4.0013664898; grid_points[ 555].w_fixed =   11.3805354935; grid_points[ 555].w_total = 9.1115392109; grid_points[ 555].i_atom = 0;
    grid_points[ 556].x =   4.0013664898; grid_points[ 556].y =  -4.0013664898; grid_points[ 556].z =   4.0013664898; grid_points[ 556].w_fixed =   11.3805354935; grid_points[ 556].w_total = 8.8295643245; grid_points[ 556].i_atom = 0;
    grid_points[ 557].x =   5.0564261420; grid_points[ 557].y =  -5.0564261420; grid_points[ 557].z =  15.1692784261; grid_points[ 557].w_fixed =  379.5745222116; grid_points[ 557].w_total = 0.0000000000; grid_points[ 557].i_atom = 0;
    grid_points[ 558].x =   5.0564261420; grid_points[ 558].y =   5.0564261420; grid_points[ 558].z = -15.1692784261; grid_points[ 558].w_fixed =  379.5745222116; grid_points[ 558].w_total = 0.0000000000; grid_points[ 558].i_atom = 0;
    grid_points[ 559].x =   5.3057605245; grid_points[ 559].y =   5.3057605245; grid_points[ 559].z =  -5.3057605245; grid_points[ 559].w_fixed =   27.4213039728; grid_points[ 559].w_total = 0.0360281880; grid_points[ 559].i_atom = 0;
    grid_points[ 560].x =   4.0013664898; grid_points[ 560].y =   4.0013664898; grid_points[ 560].z =  -4.0013664898; grid_points[ 560].w_fixed =   11.3805354935; grid_points[ 560].w_total = 9.1115392109; grid_points[ 560].i_atom = 0;
    grid_points[ 561].x =   4.0013664898; grid_points[ 561].y =   4.0013664898; grid_points[ 561].z =   4.0013664898; grid_points[ 561].w_fixed =   11.3805354935; grid_points[ 561].w_total = 8.8295643245; grid_points[ 561].i_atom = 0;
    grid_points[ 562].x =   5.3057605245; grid_points[ 562].y =   5.3057605245; grid_points[ 562].z =   5.3057605245; grid_points[ 562].w_fixed =   27.4213039728; grid_points[ 562].w_total = 0.0257436288; grid_points[ 562].i_atom = 0;
    grid_points[ 563].x =   5.0564261420; grid_points[ 563].y =   5.0564261420; grid_points[ 563].z =  15.1692784261; grid_points[ 563].w_fixed =  379.5745222116; grid_points[ 563].w_total = 0.0000000000; grid_points[ 563].i_atom = 0;
    grid_points[ 564].x =   5.0564261420; grid_points[ 564].y =  15.1692784261; grid_points[ 564].z =  -5.0564261420; grid_points[ 564].w_fixed =  379.5745222116; grid_points[ 564].w_total = 0.0000000000; grid_points[ 564].i_atom = 0;
    grid_points[ 565].x =   5.0564261420; grid_points[ 565].y =  15.1692784261; grid_points[ 565].z =   5.0564261420; grid_points[ 565].w_fixed =  379.5745222116; grid_points[ 565].w_total = 0.0000000000; grid_points[ 565].i_atom = 0;
    grid_points[ 566].x =  13.4490373684; grid_points[ 566].y = -13.4490373684; grid_points[ 566].z = -13.4490373684; grid_points[ 566].w_fixed = 1142.2409452901; grid_points[ 566].w_total = 0.0000000000; grid_points[ 566].i_atom = 0;
    grid_points[ 567].x =  13.4490373684; grid_points[ 567].y = -13.4490373684; grid_points[ 567].z =  13.4490373684; grid_points[ 567].w_fixed = 1142.2409452901; grid_points[ 567].w_total = 0.0000000000; grid_points[ 567].i_atom = 0;
    grid_points[ 568].x =  13.4490373684; grid_points[ 568].y =  13.4490373684; grid_points[ 568].z = -13.4490373684; grid_points[ 568].w_fixed = 1142.2409452901; grid_points[ 568].w_total = 0.0000000000; grid_points[ 568].i_atom = 0;
    grid_points[ 569].x =  13.4490373684; grid_points[ 569].y =  13.4490373684; grid_points[ 569].z =  13.4490373684; grid_points[ 569].w_fixed = 1142.2409452901; grid_points[ 569].w_total = 0.0000000000; grid_points[ 569].i_atom = 0;
    grid_points[ 570].x =  15.1692784261; grid_points[ 570].y =  -5.0564261420; grid_points[ 570].z =  -5.0564261420; grid_points[ 570].w_fixed =  379.5745222116; grid_points[ 570].w_total = 0.0000000000; grid_points[ 570].i_atom = 0;
    grid_points[ 571].x =  15.1692784261; grid_points[ 571].y =  -5.0564261420; grid_points[ 571].z =   5.0564261420; grid_points[ 571].w_fixed =  379.5745222116; grid_points[ 571].w_total = 0.0000000000; grid_points[ 571].i_atom = 0;
    grid_points[ 572].x =  15.1692784261; grid_points[ 572].y =   5.0564261420; grid_points[ 572].z =  -5.0564261420; grid_points[ 572].w_fixed =  379.5745222116; grid_points[ 572].w_total = 0.0000000000; grid_points[ 572].i_atom = 0;
    grid_points[ 573].x =  15.1692784261; grid_points[ 573].y =   5.0564261420; grid_points[ 573].z =   5.0564261420; grid_points[ 573].w_fixed =  379.5745222116; grid_points[ 573].w_total = 0.0000000000; grid_points[ 573].i_atom = 0;
    grid_points[ 574].x =  -4.8752953156; grid_points[ 574].y =  -4.8752953156; grid_points[ 574].z = -10.2112860776; grid_points[ 574].w_fixed =   67.6607371617; grid_points[ 574].w_total = 0.0000000000; grid_points[ 574].i_atom = 0;
    grid_points[ 575].x =  -4.7850114467; grid_points[ 575].y =  -4.7850114467; grid_points[ 575].z =  -1.4967071411; grid_points[ 575].w_fixed =   11.5537666422; grid_points[ 575].w_total = 2.4869662759; grid_points[ 575].i_atom = 0;
    grid_points[ 576].x =  -3.6367345381; grid_points[ 576].y =  -3.6367345381; grid_points[ 576].z =  -1.1375367884; grid_points[ 576].w_fixed =    4.9575936763; grid_points[ 576].w_total = 4.3703668413; grid_points[ 576].i_atom = 0;
    grid_points[ 577].x =  -3.6363256383; grid_points[ 577].y =  -3.6363256383; grid_points[ 577].z =   7.6162691612; grid_points[ 577].w_fixed =   26.8661821836; grid_points[ 577].w_total = 0.0000343844; grid_points[ 577].i_atom = 0;
    grid_points[ 578].x =  -4.8752953156; grid_points[ 578].y =   4.8752953156; grid_points[ 578].z = -10.2112860776; grid_points[ 578].w_fixed =   67.6607371617; grid_points[ 578].w_total = 0.0000000000; grid_points[ 578].i_atom = 0;
    grid_points[ 579].x =  -4.7850114467; grid_points[ 579].y =   4.7850114467; grid_points[ 579].z =  -1.4967071411; grid_points[ 579].w_fixed =   11.5537666422; grid_points[ 579].w_total = 2.4869662759; grid_points[ 579].i_atom = 0;
    grid_points[ 580].x =  -3.6367345381; grid_points[ 580].y =   3.6367345381; grid_points[ 580].z =  -1.1375367884; grid_points[ 580].w_fixed =    4.9575936763; grid_points[ 580].w_total = 4.3703668413; grid_points[ 580].i_atom = 0;
    grid_points[ 581].x =  -3.6363256383; grid_points[ 581].y =   3.6363256383; grid_points[ 581].z =   7.6162691612; grid_points[ 581].w_fixed =   26.8661821836; grid_points[ 581].w_total = 0.0000343844; grid_points[ 581].i_atom = 0;
    grid_points[ 582].x =   4.8752953156; grid_points[ 582].y =  -4.8752953156; grid_points[ 582].z = -10.2112860776; grid_points[ 582].w_fixed =   67.6607371617; grid_points[ 582].w_total = 0.0000000000; grid_points[ 582].i_atom = 0;
    grid_points[ 583].x =   4.7850114467; grid_points[ 583].y =  -4.7850114467; grid_points[ 583].z =  -1.4967071411; grid_points[ 583].w_fixed =   11.5537666422; grid_points[ 583].w_total = 2.4869662759; grid_points[ 583].i_atom = 0;
    grid_points[ 584].x =   3.6367345381; grid_points[ 584].y =  -3.6367345381; grid_points[ 584].z =  -1.1375367884; grid_points[ 584].w_fixed =    4.9575936763; grid_points[ 584].w_total = 4.3703668413; grid_points[ 584].i_atom = 0;
    grid_points[ 585].x =   3.6363256383; grid_points[ 585].y =  -3.6363256383; grid_points[ 585].z =   7.6162691612; grid_points[ 585].w_fixed =   26.8661821836; grid_points[ 585].w_total = 0.0000343844; grid_points[ 585].i_atom = 0;
    grid_points[ 586].x =   4.8752953156; grid_points[ 586].y =   4.8752953156; grid_points[ 586].z = -10.2112860776; grid_points[ 586].w_fixed =   67.6607371617; grid_points[ 586].w_total = 0.0000000000; grid_points[ 586].i_atom = 0;
    grid_points[ 587].x =   4.7850114467; grid_points[ 587].y =   4.7850114467; grid_points[ 587].z =  -1.4967071411; grid_points[ 587].w_fixed =   11.5537666422; grid_points[ 587].w_total = 2.4869662759; grid_points[ 587].i_atom = 0;
    grid_points[ 588].x =   3.6367345381; grid_points[ 588].y =   3.6367345381; grid_points[ 588].z =  -1.1375367884; grid_points[ 588].w_fixed =    4.9575936763; grid_points[ 588].w_total = 4.3703668413; grid_points[ 588].i_atom = 0;
    grid_points[ 589].x =   3.6363256383; grid_points[ 589].y =   3.6363256383; grid_points[ 589].z =   7.6162691612; grid_points[ 589].w_fixed =   26.8661821836; grid_points[ 589].w_total = 0.0000343844; grid_points[ 589].i_atom = 0;
    grid_points[ 590].x =  -5.8939912332; grid_points[ 590].y = -10.8198077400; grid_points[ 590].z =   0.0000000000; grid_points[ 590].w_fixed =   68.3625199466; grid_points[ 590].w_total = 0.0000000000; grid_points[ 590].i_atom = 0;
    grid_points[ 591].x =  -5.7438484187; grid_points[ 591].y =  -2.7423536150; grid_points[ 591].z =   2.7423536150; grid_points[ 591].w_fixed =   11.1501458946; grid_points[ 591].w_total = 0.7753792790; grid_points[ 591].i_atom = 0;
    grid_points[ 592].x =  -5.8939912332; grid_points[ 592].y =   0.0000000000; grid_points[ 592].z =  10.8198077400; grid_points[ 592].w_fixed =   68.3625199466; grid_points[ 592].w_total = 0.0000000000; grid_points[ 592].i_atom = 0;
    grid_points[ 593].x =  -3.3342974947; grid_points[ 593].y =  -1.5919331643; grid_points[ 593].z =   1.5919331643; grid_points[ 593].w_fixed =    2.1032476402; grid_points[ 593].w_total = 2.0194520110; grid_points[ 593].i_atom = 0;
    grid_points[ 594].x =  -3.5329984457; grid_points[ 594].y =  -1.9245685660; grid_points[ 594].z =   0.0000000000; grid_points[ 594].w_fixed =    2.1250626994; grid_points[ 594].w_total = 1.9718171965; grid_points[ 594].i_atom = 0;
    grid_points[ 595].x =  -4.3654758527; grid_points[ 595].y =  -2.0842608672; grid_points[ 595].z =   2.0842608672; grid_points[ 595].w_fixed =    4.7844044707; grid_points[ 595].w_total = 3.3058550217; grid_points[ 595].i_atom = 0;
    grid_points[ 596].x =  -4.6256278652; grid_points[ 596].y =  -2.5197684414; grid_points[ 596].z =   0.0000000000; grid_points[ 596].w_fixed =    4.8340287112; grid_points[ 596].w_total = 2.6697241981; grid_points[ 596].i_atom = 0;
    grid_points[ 597].x =  -3.8828754963; grid_points[ 597].y =  -0.7447550142; grid_points[ 597].z =   0.7447550142; grid_points[ 597].w_fixed =    1.7999446257; grid_points[ 597].w_total = 1.5204246946; grid_points[ 597].i_atom = 0;
    grid_points[ 598].x =  -3.5329984457; grid_points[ 598].y =   0.0000000000; grid_points[ 598].z =   1.9245685660; grid_points[ 598].w_fixed =    2.1250626994; grid_points[ 598].w_total = 1.9844543174; grid_points[ 598].i_atom = 0;
    grid_points[ 599].x =  -5.0837093106; grid_points[ 599].y =  -0.9750809686; grid_points[ 599].z =   0.9750809686; grid_points[ 599].w_fixed =    4.0944598961; grid_points[ 599].w_total = 1.3604510062; grid_points[ 599].i_atom = 0;
    grid_points[ 600].x =  -4.6256278652; grid_points[ 600].y =   0.0000000000; grid_points[ 600].z =   2.5197684414; grid_points[ 600].w_fixed =    4.8340287112; grid_points[ 600].w_total = 2.6802181871; grid_points[ 600].i_atom = 0;
    grid_points[ 601].x =  -4.3961380893; grid_points[ 601].y =   8.0701458556; grid_points[ 601].z =   0.0000000000; grid_points[ 601].w_fixed =   27.1448404564; grid_points[ 601].w_total = 0.0003003544; grid_points[ 601].i_atom = 0;
    grid_points[ 602].x =   5.8939912332; grid_points[ 602].y = -10.8198077400; grid_points[ 602].z =   0.0000000000; grid_points[ 602].w_fixed =   68.3625199466; grid_points[ 602].w_total = 0.0000000000; grid_points[ 602].i_atom = 0;
    grid_points[ 603].x =   3.3342974947; grid_points[ 603].y =  -1.5919331643; grid_points[ 603].z =   1.5919331643; grid_points[ 603].w_fixed =    2.1032476402; grid_points[ 603].w_total = 2.0194520110; grid_points[ 603].i_atom = 0;
    grid_points[ 604].x =   3.5329984457; grid_points[ 604].y =  -1.9245685660; grid_points[ 604].z =   0.0000000000; grid_points[ 604].w_fixed =    2.1250626994; grid_points[ 604].w_total = 1.9718171965; grid_points[ 604].i_atom = 0;
    grid_points[ 605].x =   4.3654758527; grid_points[ 605].y =  -2.0842608672; grid_points[ 605].z =   2.0842608672; grid_points[ 605].w_fixed =    4.7844044707; grid_points[ 605].w_total = 3.3058550217; grid_points[ 605].i_atom = 0;
    grid_points[ 606].x =   4.6256278652; grid_points[ 606].y =  -2.5197684414; grid_points[ 606].z =   0.0000000000; grid_points[ 606].w_fixed =    4.8340287112; grid_points[ 606].w_total = 2.6697241981; grid_points[ 606].i_atom = 0;
    grid_points[ 607].x =   5.7438484187; grid_points[ 607].y =  -2.7423536150; grid_points[ 607].z =   2.7423536150; grid_points[ 607].w_fixed =   11.1501458946; grid_points[ 607].w_total = 0.7753792790; grid_points[ 607].i_atom = 0;
    grid_points[ 608].x =   3.8828754963; grid_points[ 608].y =  -0.7447550142; grid_points[ 608].z =   0.7447550142; grid_points[ 608].w_fixed =    1.7999446257; grid_points[ 608].w_total = 1.5204246946; grid_points[ 608].i_atom = 0;
    grid_points[ 609].x =   3.5329984457; grid_points[ 609].y =   0.0000000000; grid_points[ 609].z =   1.9245685660; grid_points[ 609].w_fixed =    2.1250626994; grid_points[ 609].w_total = 1.9844543174; grid_points[ 609].i_atom = 0;
    grid_points[ 610].x =   5.0837093106; grid_points[ 610].y =  -0.9750809686; grid_points[ 610].z =   0.9750809686; grid_points[ 610].w_fixed =    4.0944598961; grid_points[ 610].w_total = 1.3604510062; grid_points[ 610].i_atom = 0;
    grid_points[ 611].x =   4.6256278652; grid_points[ 611].y =   0.0000000000; grid_points[ 611].z =   2.5197684414; grid_points[ 611].w_fixed =    4.8340287112; grid_points[ 611].w_total = 2.6802181871; grid_points[ 611].i_atom = 0;
    grid_points[ 612].x =   5.8939912332; grid_points[ 612].y =   0.0000000000; grid_points[ 612].z =  10.8198077400; grid_points[ 612].w_fixed =   68.3625199466; grid_points[ 612].w_total = 0.0000000000; grid_points[ 612].i_atom = 0;
    grid_points[ 613].x =   4.3961380893; grid_points[ 613].y =   8.0701458556; grid_points[ 613].z =   0.0000000000; grid_points[ 613].w_fixed =   27.1448404564; grid_points[ 613].w_total = 0.0003003544; grid_points[ 613].i_atom = 0;
    grid_points[ 614].x = -23.5142447435; grid_points[ 614].y =   0.0000000000; grid_points[ 614].z =  23.5142447435; grid_points[ 614].w_fixed = 3887.7671433090; grid_points[ 614].w_total = 0.0000000000; grid_points[ 614].i_atom = 0;
    grid_points[ 615].x =  -6.0861418537; grid_points[ 615].y =   0.0000000000; grid_points[ 615].z =   3.3153700687; grid_points[ 615].w_fixed =   11.2657961338; grid_points[ 615].w_total = 0.3471224190; grid_points[ 615].i_atom = 0;
    grid_points[ 616].x =  -4.8752953156; grid_points[ 616].y = -10.2112860776; grid_points[ 616].z =  -4.8752953156; grid_points[ 616].w_fixed =   67.6607371617; grid_points[ 616].w_total = 0.0000000000; grid_points[ 616].i_atom = 0;
    grid_points[ 617].x =  -4.8752953156; grid_points[ 617].y = -10.2112860776; grid_points[ 617].z =   4.8752953156; grid_points[ 617].w_fixed =   67.6607371617; grid_points[ 617].w_total = 0.0000000000; grid_points[ 617].i_atom = 0;
    grid_points[ 618].x =  -4.7850114467; grid_points[ 618].y =  -1.4967071411; grid_points[ 618].z =  -4.7850114467; grid_points[ 618].w_fixed =   11.5537666422; grid_points[ 618].w_total = 2.9255841962; grid_points[ 618].i_atom = 0;
    grid_points[ 619].x =  -4.7850114467; grid_points[ 619].y =  -1.4967071411; grid_points[ 619].z =   4.7850114467; grid_points[ 619].w_fixed =   11.5537666422; grid_points[ 619].w_total = 1.8853671084; grid_points[ 619].i_atom = 0;
    grid_points[ 620].x =  -3.6367345381; grid_points[ 620].y =  -1.1375367884; grid_points[ 620].z =  -3.6367345381; grid_points[ 620].w_fixed =    4.9575936763; grid_points[ 620].w_total = 4.6559094346; grid_points[ 620].i_atom = 0;
    grid_points[ 621].x =  -3.6367345381; grid_points[ 621].y =  -1.1375367884; grid_points[ 621].z =   3.6367345381; grid_points[ 621].w_fixed =    4.9575936763; grid_points[ 621].w_total = 4.5221022069; grid_points[ 621].i_atom = 0;
    grid_points[ 622].x =  -3.3153700687; grid_points[ 622].y =   0.0000000000; grid_points[ 622].z =   6.0861418537; grid_points[ 622].w_fixed =   11.2657961338; grid_points[ 622].w_total = 0.0045977423; grid_points[ 622].i_atom = 0;
    grid_points[ 623].x =  -3.6363256383; grid_points[ 623].y =   7.6162691612; grid_points[ 623].z =  -3.6363256383; grid_points[ 623].w_fixed =   26.8661821836; grid_points[ 623].w_total = 0.0010773582; grid_points[ 623].i_atom = 0;
    grid_points[ 624].x =  -3.6363256383; grid_points[ 624].y =   7.6162691612; grid_points[ 624].z =   3.6363256383; grid_points[ 624].w_fixed =   26.8661821836; grid_points[ 624].w_total = 0.0009330210; grid_points[ 624].i_atom = 0;
    grid_points[ 625].x =   4.8752953156; grid_points[ 625].y = -10.2112860776; grid_points[ 625].z =  -4.8752953156; grid_points[ 625].w_fixed =   67.6607371617; grid_points[ 625].w_total = 0.0000000000; grid_points[ 625].i_atom = 0;
    grid_points[ 626].x =   4.8752953156; grid_points[ 626].y = -10.2112860776; grid_points[ 626].z =   4.8752953156; grid_points[ 626].w_fixed =   67.6607371617; grid_points[ 626].w_total = 0.0000000000; grid_points[ 626].i_atom = 0;
    grid_points[ 627].x =   4.7850114467; grid_points[ 627].y =  -1.4967071411; grid_points[ 627].z =  -4.7850114467; grid_points[ 627].w_fixed =   11.5537666422; grid_points[ 627].w_total = 2.9255841962; grid_points[ 627].i_atom = 0;
    grid_points[ 628].x =   4.7850114467; grid_points[ 628].y =  -1.4967071411; grid_points[ 628].z =   4.7850114467; grid_points[ 628].w_fixed =   11.5537666422; grid_points[ 628].w_total = 1.8853671084; grid_points[ 628].i_atom = 0;
    grid_points[ 629].x =   3.6367345381; grid_points[ 629].y =  -1.1375367884; grid_points[ 629].z =  -3.6367345381; grid_points[ 629].w_fixed =    4.9575936763; grid_points[ 629].w_total = 4.6559094346; grid_points[ 629].i_atom = 0;
    grid_points[ 630].x =   3.6367345381; grid_points[ 630].y =  -1.1375367884; grid_points[ 630].z =   3.6367345381; grid_points[ 630].w_fixed =    4.9575936763; grid_points[ 630].w_total = 4.5221022069; grid_points[ 630].i_atom = 0;
    grid_points[ 631].x =   3.3153700687; grid_points[ 631].y =   0.0000000000; grid_points[ 631].z =   6.0861418537; grid_points[ 631].w_fixed =   11.2657961338; grid_points[ 631].w_total = 0.0045977423; grid_points[ 631].i_atom = 0;
    grid_points[ 632].x =   6.0861418537; grid_points[ 632].y =   0.0000000000; grid_points[ 632].z =   3.3153700687; grid_points[ 632].w_fixed =   11.2657961338; grid_points[ 632].w_total = 0.3471224190; grid_points[ 632].i_atom = 0;
    grid_points[ 633].x =   3.6363256383; grid_points[ 633].y =   7.6162691612; grid_points[ 633].z =  -3.6363256383; grid_points[ 633].w_fixed =   26.8661821836; grid_points[ 633].w_total = 0.0010773582; grid_points[ 633].i_atom = 0;
    grid_points[ 634].x =   3.6363256383; grid_points[ 634].y =   7.6162691612; grid_points[ 634].z =   3.6363256383; grid_points[ 634].w_fixed =   26.8661821836; grid_points[ 634].w_total = 0.0009330210; grid_points[ 634].i_atom = 0;
    grid_points[ 635].x =  23.5142447435; grid_points[ 635].y =   0.0000000000; grid_points[ 635].z =  23.5142447435; grid_points[ 635].w_fixed = 3887.7671433090; grid_points[ 635].w_total = 0.0000000000; grid_points[ 635].i_atom = 0;
    grid_points[ 636].x =  -5.7438484187; grid_points[ 636].y =  -2.7423536150; grid_points[ 636].z =  -2.7423536150; grid_points[ 636].w_fixed =   11.1501458946; grid_points[ 636].w_total = 0.7975816561; grid_points[ 636].i_atom = 0;
    grid_points[ 637].x =  -3.3342974947; grid_points[ 637].y =  -1.5919331643; grid_points[ 637].z =  -1.5919331643; grid_points[ 637].w_fixed =    2.1032476402; grid_points[ 637].w_total = 2.0149126278; grid_points[ 637].i_atom = 0;
    grid_points[ 638].x =  -4.3654758527; grid_points[ 638].y =  -2.0842608672; grid_points[ 638].z =  -2.0842608672; grid_points[ 638].w_fixed =    4.7844044707; grid_points[ 638].w_total = 3.2966811514; grid_points[ 638].i_atom = 0;
    grid_points[ 639].x =  -3.8828754963; grid_points[ 639].y =  -0.7447550142; grid_points[ 639].z =  -0.7447550142; grid_points[ 639].w_fixed =    1.7999446257; grid_points[ 639].w_total = 1.5182141465; grid_points[ 639].i_atom = 0;
    grid_points[ 640].x =  -5.0837093106; grid_points[ 640].y =  -0.9750809686; grid_points[ 640].z =  -0.9750809686; grid_points[ 640].w_fixed =    4.0944598961; grid_points[ 640].w_total = 1.3639957134; grid_points[ 640].i_atom = 0;
    grid_points[ 641].x =  -4.3961380893; grid_points[ 641].y =   0.0000000000; grid_points[ 641].z =   8.0701458556; grid_points[ 641].w_fixed =   27.1448404564; grid_points[ 641].w_total = 0.0000008316; grid_points[ 641].i_atom = 0;
    grid_points[ 642].x =   3.3342974947; grid_points[ 642].y =  -1.5919331643; grid_points[ 642].z =  -1.5919331643; grid_points[ 642].w_fixed =    2.1032476402; grid_points[ 642].w_total = 2.0149126278; grid_points[ 642].i_atom = 0;
    grid_points[ 643].x =   4.3654758527; grid_points[ 643].y =  -2.0842608672; grid_points[ 643].z =  -2.0842608672; grid_points[ 643].w_fixed =    4.7844044707; grid_points[ 643].w_total = 3.2966811514; grid_points[ 643].i_atom = 0;
    grid_points[ 644].x =   5.7438484187; grid_points[ 644].y =  -2.7423536150; grid_points[ 644].z =  -2.7423536150; grid_points[ 644].w_fixed =   11.1501458946; grid_points[ 644].w_total = 0.7975816561; grid_points[ 644].i_atom = 0;
    grid_points[ 645].x =   3.8828754963; grid_points[ 645].y =  -0.7447550142; grid_points[ 645].z =  -0.7447550142; grid_points[ 645].w_fixed =    1.7999446257; grid_points[ 645].w_total = 1.5182141465; grid_points[ 645].i_atom = 0;
    grid_points[ 646].x =   5.0837093106; grid_points[ 646].y =  -0.9750809686; grid_points[ 646].z =  -0.9750809686; grid_points[ 646].w_fixed =    4.0944598961; grid_points[ 646].w_total = 1.3639957134; grid_points[ 646].i_atom = 0;
    grid_points[ 647].x =   4.3961380893; grid_points[ 647].y =   0.0000000000; grid_points[ 647].z =   8.0701458556; grid_points[ 647].w_fixed =   27.1448404564; grid_points[ 647].w_total = 0.0000008316; grid_points[ 647].i_atom = 0;
    grid_points[ 648].x = -30.0795229432; grid_points[ 648].y =  10.0265076477; grid_points[ 648].z =  10.0265076477; grid_points[ 648].w_fixed = 3474.1698452873; grid_points[ 648].w_total = 0.0000000000; grid_points[ 648].i_atom = 0;
    grid_points[ 649].x = -21.0705920995; grid_points[ 649].y =  -7.0235306998; grid_points[ 649].z =  -7.0235306998; grid_points[ 649].w_fixed = 1092.3998745831; grid_points[ 649].w_total = 0.0000000000; grid_points[ 649].i_atom = 0;
    grid_points[ 650].x = -19.1993004361; grid_points[ 650].y =  19.1993004361; grid_points[ 650].z =  19.1993004361; grid_points[ 650].w_fixed = 3632.6798826243; grid_points[ 650].w_total = 0.0000000000; grid_points[ 650].i_atom = 0;
    grid_points[ 651].x = -12.3210134399; grid_points[ 651].y =   0.0000000000; grid_points[ 651].z =   0.0000000000; grid_points[ 651].w_fixed =   26.9943595963; grid_points[ 651].w_total = 0.0000000000; grid_points[ 651].i_atom = 0;
    grid_points[ 652].x = -11.8913062075; grid_points[ 652].y =   2.2808122310; grid_points[ 652].z =   2.2808122310; grid_points[ 652].w_fixed =   57.9035858166; grid_points[ 652].w_total = 0.0000000000; grid_points[ 652].i_atom = 0;
    grid_points[ 653].x =  -9.6823189135; grid_points[ 653].y =   9.6823189135; grid_points[ 653].z =   9.6823189135; grid_points[ 653].w_fixed =  396.8927232114; grid_points[ 653].w_total = 0.0000000000; grid_points[ 653].i_atom = 0;
    grid_points[ 654].x = -10.0265076477; grid_points[ 654].y =  10.0265076477; grid_points[ 654].z =  30.0795229432; grid_points[ 654].w_fixed = 3474.1698452873; grid_points[ 654].w_total = 0.0000000000; grid_points[ 654].i_atom = 0;
    grid_points[ 655].x = -11.8583704325; grid_points[ 655].y =  11.8583704325; grid_points[ 655].z =   0.0000000000; grid_points[ 655].w_fixed =  424.7625826047; grid_points[ 655].w_total = 0.0000000000; grid_points[ 655].i_atom = 0;
    grid_points[ 656].x = -10.0265076477; grid_points[ 656].y =  30.0795229432; grid_points[ 656].z =  10.0265076477; grid_points[ 656].w_fixed = 3474.1698452873; grid_points[ 656].w_total = 0.0000000000; grid_points[ 656].i_atom = 0;
    grid_points[ 657].x =  -1.9846144277; grid_points[ 657].y =  -6.3448636630; grid_points[ 657].z =  -6.3448636630; grid_points[ 657].w_fixed =   27.8387029598; grid_points[ 657].w_total = 0.0002691714; grid_points[ 657].i_atom = 0;
    grid_points[ 658].x =  -2.6608126953; grid_points[ 658].y =  -8.5066870163; grid_points[ 658].z =  -8.5066870163; grid_points[ 658].w_fixed =   70.1099676542; grid_points[ 658].w_total = 0.0000000000; grid_points[ 658].i_atom = 0;
    grid_points[ 659].x =  -1.2829568800; grid_points[ 659].y =  -6.6888597417; grid_points[ 659].z =   1.2829568800; grid_points[ 659].w_fixed =    9.5422169009; grid_points[ 659].w_total = 0.0783957863; grid_points[ 659].i_atom = 0;
    grid_points[ 660].x =  -1.7011843293; grid_points[ 660].y =  -8.8693420267; grid_points[ 660].z =   1.7011843293; grid_points[ 660].w_fixed =   22.9918908793; grid_points[ 660].w_total = 0.0000018741; grid_points[ 660].i_atom = 0;
    grid_points[ 661].x =  -1.2829568800; grid_points[ 661].y =   1.2829568800; grid_points[ 661].z =  -6.6888597417; grid_points[ 661].w_fixed =    9.5422169009; grid_points[ 661].w_total = 0.0000055082; grid_points[ 661].i_atom = 0;
    grid_points[ 662].x =  -1.7011843293; grid_points[ 662].y =   1.7011843293; grid_points[ 662].z =  -8.8693420267; grid_points[ 662].w_fixed =   22.9918908793; grid_points[ 662].w_total = 0.0000000000; grid_points[ 662].i_atom = 0;
    grid_points[ 663].x =  -1.3690014933; grid_points[ 663].y =   0.0000000000; grid_points[ 663].z =   0.0000000000; grid_points[ 663].w_fixed =    0.4264811437; grid_points[ 663].w_total = 0.4264788573; grid_points[ 663].i_atom = 0;
    grid_points[ 664].x =  -1.8012554250; grid_points[ 664].y =   0.0000000000; grid_points[ 664].z =   0.0000000000; grid_points[ 664].w_fixed =    0.2531886212; grid_points[ 664].w_total = 0.2531617203; grid_points[ 664].i_atom = 0;
    grid_points[ 665].x =  -1.2736799257; grid_points[ 665].y =   1.2736799257; grid_points[ 665].z =   0.0000000000; grid_points[ 665].w_fixed =    0.4501131044; grid_points[ 665].w_total = 0.4501076870; grid_points[ 665].i_atom = 0;
    grid_points[ 666].x =  -1.6292968354; grid_points[ 666].y =   0.5430989451; grid_points[ 666].z =   0.5430989451; grid_points[ 666].w_fixed =    0.4022281471; grid_points[ 666].w_total = 0.4022123636; grid_points[ 666].i_atom = 0;
    grid_points[ 667].x =  -2.3583189787; grid_points[ 667].y =   0.0000000000; grid_points[ 667].z =   0.0000000000; grid_points[ 667].w_fixed =    0.5606220456; grid_points[ 667].w_total = 0.5595902681; grid_points[ 667].i_atom = 0;
    grid_points[ 668].x =  -1.6675833421; grid_points[ 668].y =   1.6675833421; grid_points[ 668].z =   0.0000000000; grid_points[ 668].w_fixed =    0.9966614143; grid_points[ 668].w_total = 0.9964922905; grid_points[ 668].i_atom = 0;
    grid_points[ 669].x =  -1.3615760972; grid_points[ 669].y =   1.3615760972; grid_points[ 669].z =   1.3615760972; grid_points[ 669].w_fixed =    0.9312676753; grid_points[ 669].w_total = 0.9312229918; grid_points[ 669].i_atom = 0;
    grid_points[ 670].x =  -2.1331797787; grid_points[ 670].y =   0.7110599262; grid_points[ 670].z =   0.7110599262; grid_points[ 670].w_fixed =    0.8906323100; grid_points[ 670].w_total = 0.8900464448; grid_points[ 670].i_atom = 0;
    grid_points[ 671].x =  -3.0802533600; grid_points[ 671].y =   0.0000000000; grid_points[ 671].z =   0.0000000000; grid_points[ 671].w_fixed =    0.3749216611; grid_points[ 671].w_total = 0.3657764259; grid_points[ 671].i_atom = 0;
    grid_points[ 672].x =  -1.7783851066; grid_points[ 672].y =   1.7783851066; grid_points[ 672].z =   1.7783851066; grid_points[ 672].w_fixed =    0.9591496581; grid_points[ 672].w_total = 0.9586300412; grid_points[ 672].i_atom = 0;
    grid_points[ 673].x =  -2.9728265519; grid_points[ 673].y =   0.5702030578; grid_points[ 673].z =   0.5702030578; grid_points[ 673].w_fixed =    0.8042164697; grid_points[ 673].w_total = 0.7903823761; grid_points[ 673].i_atom = 0;
    grid_points[ 674].x =  -2.1266717541; grid_points[ 674].y =   2.1266717541; grid_points[ 674].z =   0.6652031738; grid_points[ 674].w_fixed =    0.9737495508; grid_points[ 674].w_total = 0.9722784395; grid_points[ 674].i_atom = 0;
    grid_points[ 675].x =  -2.1266717541; grid_points[ 675].y =   0.6652031738; grid_points[ 675].z =   2.1266717541; grid_points[ 675].w_fixed =    0.9737495508; grid_points[ 675].w_total = 0.9690650946; grid_points[ 675].i_atom = 0;
    grid_points[ 676].x =  -2.5528215194; grid_points[ 676].y =   1.2188238289; grid_points[ 676].z =   1.2188238289; grid_points[ 676].w_fixed =    0.9397324606; grid_points[ 676].w_total = 0.9361991824; grid_points[ 676].i_atom = 0;
    grid_points[ 677].x =  -1.4734978083; grid_points[ 677].y =   2.7049519350; grid_points[ 677].z =   0.0000000000; grid_points[ 677].w_fixed =    0.9494794437; grid_points[ 677].w_total = 0.9424324035; grid_points[ 677].i_atom = 0;
    grid_points[ 678].x =  -2.7049519350; grid_points[ 678].y =   1.4734978083; grid_points[ 678].z =   0.0000000000; grid_points[ 678].w_fixed =    0.9494794437; grid_points[ 678].w_total = 0.9424324035; grid_points[ 678].i_atom = 0;
    grid_points[ 679].x =  -2.3227887106; grid_points[ 679].y =   2.3227887106; grid_points[ 679].z =   2.3227887106; grid_points[ 679].w_fixed =    2.1467059397; grid_points[ 679].w_total = 2.1405170372; grid_points[ 679].i_atom = 0;
    grid_points[ 680].x =  -2.7776937196; grid_points[ 680].y =   0.8688367985; grid_points[ 680].z =   2.7776937196; grid_points[ 680].w_fixed =    2.1793824633; grid_points[ 680].w_total = 2.1192111689; grid_points[ 680].i_atom = 0;
    grid_points[ 681].x =  -2.2808122310; grid_points[ 681].y =   2.2808122310; grid_points[ 681].z =  11.8913062075; grid_points[ 681].w_fixed =   57.9035858166; grid_points[ 681].w_total = 0.0000000000; grid_points[ 681].i_atom = 0;
    grid_points[ 682].x =  -2.7776937196; grid_points[ 682].y =   2.7776937196; grid_points[ 682].z =   0.8688367985; grid_points[ 682].w_fixed =    2.1793824633; grid_points[ 682].w_total = 2.1482231938; grid_points[ 682].i_atom = 0;
    grid_points[ 683].x =  -3.0411437621; grid_points[ 683].y =   3.0411437621; grid_points[ 683].z =   3.0411437621; grid_points[ 683].w_fixed =    4.8832621034; grid_points[ 683].w_total = 4.8025419330; grid_points[ 683].i_atom = 0;
    grid_points[ 684].x =  -2.2808122310; grid_points[ 684].y =  11.8913062075; grid_points[ 684].z =   2.2808122310; grid_points[ 684].w_fixed =   57.9035858166; grid_points[ 684].w_total = 0.0000000000; grid_points[ 684].i_atom = 0;
    grid_points[ 685].x =  -0.0036626080; grid_points[ 685].y =   0.0000000000; grid_points[ 685].z =   0.0000000000; grid_points[ 685].w_fixed =    0.0000002129; grid_points[ 685].w_total = 0.0000002129; grid_points[ 685].i_atom = 0;
    grid_points[ 686].x =   0.0000000000; grid_points[ 686].y =   0.0036626080; grid_points[ 686].z =   0.0000000000; grid_points[ 686].w_fixed =    0.0000002129; grid_points[ 686].w_total = 0.0000002129; grid_points[ 686].i_atom = 0;
    grid_points[ 687].x =  -0.0157155784; grid_points[ 687].y =   0.0000000000; grid_points[ 687].z =   0.0000000000; grid_points[ 687].w_fixed =    0.0000087099; grid_points[ 687].w_total = 0.0000087099; grid_points[ 687].i_atom = 0;
    grid_points[ 688].x =   0.0000000000; grid_points[ 688].y =   0.0157155784; grid_points[ 688].z =   0.0000000000; grid_points[ 688].w_fixed =    0.0000087099; grid_points[ 688].w_total = 0.0000087099; grid_points[ 688].i_atom = 0;
    grid_points[ 689].x =  -0.0380278193; grid_points[ 689].y =   0.0000000000; grid_points[ 689].z =   0.0000000000; grid_points[ 689].w_fixed =    0.0000853157; grid_points[ 689].w_total = 0.0000853157; grid_points[ 689].i_atom = 0;
    grid_points[ 690].x =   0.0000000000; grid_points[ 690].y =   0.0380278193; grid_points[ 690].z =   0.0000000000; grid_points[ 690].w_fixed =    0.0000853157; grid_points[ 690].w_total = 0.0000853157; grid_points[ 690].i_atom = 0;
    grid_points[ 691].x =  -0.0729054050; grid_points[ 691].y =   0.0000000000; grid_points[ 691].z =   0.0000000000; grid_points[ 691].w_fixed =    0.0004682262; grid_points[ 691].w_total = 0.0004682262; grid_points[ 691].i_atom = 0;
    grid_points[ 692].x =   0.0000000000; grid_points[ 692].y =   0.0729054050; grid_points[ 692].z =   0.0000000000; grid_points[ 692].w_fixed =    0.0004682262; grid_points[ 692].w_total = 0.0004682262; grid_points[ 692].i_atom = 0;
    grid_points[ 693].x =  -0.1232101344; grid_points[ 693].y =   0.0000000000; grid_points[ 693].z =   0.0000000000; grid_points[ 693].w_fixed =    0.0018803520; grid_points[ 693].w_total = 0.0018803520; grid_points[ 693].i_atom = 0;
    grid_points[ 694].x =   0.0000000000; grid_points[ 694].y =   0.1232101344; grid_points[ 694].z =   0.0000000000; grid_points[ 694].w_fixed =    0.0018803520; grid_points[ 694].w_total = 0.0018803520; grid_points[ 694].i_atom = 0;
    grid_points[ 695].x =  -0.1925158350; grid_points[ 695].y =   0.0000000000; grid_points[ 695].z =   0.0000000000; grid_points[ 695].w_fixed =    0.0062265394; grid_points[ 695].w_total = 0.0062265394; grid_points[ 695].i_atom = 0;
    grid_points[ 696].x =   0.0000000000; grid_points[ 696].y =   0.1925158350; grid_points[ 696].z =   0.0000000000; grid_points[ 696].w_fixed =    0.0062265394; grid_points[ 696].w_total = 0.0062265394; grid_points[ 696].i_atom = 0;
    grid_points[ 697].x =  -0.2853164738; grid_points[ 697].y =   0.0000000000; grid_points[ 697].z =   0.0000000000; grid_points[ 697].w_fixed =    0.0181286111; grid_points[ 697].w_total = 0.0181286111; grid_points[ 697].i_atom = 0;
    grid_points[ 698].x =   0.0000000000; grid_points[ 698].y =   0.2853164738; grid_points[ 698].z =   0.0000000000; grid_points[ 698].w_fixed =    0.0181286111; grid_points[ 698].w_total = 0.0181286111; grid_points[ 698].i_atom = 0;
    grid_points[ 699].x =  -0.4073062294; grid_points[ 699].y =   0.0000000000; grid_points[ 699].z =   0.0000000000; grid_points[ 699].w_fixed =    0.0482458934; grid_points[ 699].w_total = 0.0482458934; grid_points[ 699].i_atom = 0;
    grid_points[ 700].x =   0.0000000000; grid_points[ 700].y =   0.4073062294; grid_points[ 700].z =   0.0000000000; grid_points[ 700].w_fixed =    0.0482458934; grid_points[ 700].w_total = 0.0482458934; grid_points[ 700].i_atom = 0;
    grid_points[ 701].x =  -0.5657608212; grid_points[ 701].y =   0.0000000000; grid_points[ 701].z =   0.0000000000; grid_points[ 701].w_fixed =    0.0344015986; grid_points[ 701].w_total = 0.0344015985; grid_points[ 701].i_atom = 0;
    grid_points[ 702].x =   0.0000000000; grid_points[ 702].y =   0.5657608212; grid_points[ 702].z =   0.0000000000; grid_points[ 702].w_fixed =    0.0344015986; grid_points[ 702].w_total = 0.0344015985; grid_points[ 702].i_atom = 0;
    grid_points[ 703].x =   0.0000000000; grid_points[ 703].y =   0.4000533132; grid_points[ 703].z =   0.4000533132; grid_points[ 703].w_fixed =    0.0275212788; grid_points[ 703].w_total = 0.0275212788; grid_points[ 703].i_atom = 0;
    grid_points[ 704].x =  -0.4000533132; grid_points[ 704].y =   0.4000533132; grid_points[ 704].z =   0.0000000000; grid_points[ 704].w_fixed =    0.0275212788; grid_points[ 704].w_total = 0.0275212788; grid_points[ 704].i_atom = 0;
    grid_points[ 705].x =  -0.3266421624; grid_points[ 705].y =   0.3266421624; grid_points[ 705].z =   0.3266421624; grid_points[ 705].w_fixed =    0.0232210790; grid_points[ 705].w_total = 0.0232210790; grid_points[ 705].i_atom = 0;
    grid_points[ 706].x =  -0.7700633400; grid_points[ 706].y =   0.0000000000; grid_points[ 706].z =   0.0000000000; grid_points[ 706].w_fixed =    0.0819768394; grid_points[ 706].w_total = 0.0819768384; grid_points[ 706].i_atom = 0;
    grid_points[ 707].x =   0.0000000000; grid_points[ 707].y =   0.7700633400; grid_points[ 707].z =   0.0000000000; grid_points[ 707].w_fixed =    0.0819768394; grid_points[ 707].w_total = 0.0819768384; grid_points[ 707].i_atom = 0;
    grid_points[ 708].x =   0.0000000000; grid_points[ 708].y =   0.5445170097; grid_points[ 708].z =   0.5445170097; grid_points[ 708].w_fixed =    0.0655814715; grid_points[ 708].w_total = 0.0655814714; grid_points[ 708].i_atom = 0;
    grid_points[ 709].x =  -0.5445170097; grid_points[ 709].y =   0.5445170097; grid_points[ 709].z =   0.0000000000; grid_points[ 709].w_fixed =    0.0655814715; grid_points[ 709].w_total = 0.0655814713; grid_points[ 709].i_atom = 0;
    grid_points[ 710].x =  -0.4445962766; grid_points[ 710].y =   0.4445962766; grid_points[ 710].z =   0.4445962766; grid_points[ 710].w_fixed =    0.0553343666; grid_points[ 710].w_total = 0.0553343665; grid_points[ 710].i_atom = 0;
    grid_points[ 711].x =  -1.0324394919; grid_points[ 711].y =   0.0000000000; grid_points[ 711].z =   0.0000000000; grid_points[ 711].w_fixed =    0.1890555184; grid_points[ 711].w_total = 0.1890554699; grid_points[ 711].i_atom = 0;
    grid_points[ 712].x =   0.0000000000; grid_points[ 712].y =   1.0324394919; grid_points[ 712].z =   0.0000000000; grid_points[ 712].w_fixed =    0.1890555184; grid_points[ 712].w_total = 0.1890554699; grid_points[ 712].i_atom = 0;
    grid_points[ 713].x =   0.0000000000; grid_points[ 713].y =   0.7300449659; grid_points[ 713].z =   0.7300449659; grid_points[ 713].w_fixed =    0.1512444148; grid_points[ 713].w_total = 0.1512444119; grid_points[ 713].i_atom = 0;
    grid_points[ 714].x =  -0.7300449659; grid_points[ 714].y =   0.7300449659; grid_points[ 714].z =   0.0000000000; grid_points[ 714].w_fixed =    0.1512444148; grid_points[ 714].w_total = 0.1512444075; grid_points[ 714].i_atom = 0;
    grid_points[ 715].x =  -0.5960792185; grid_points[ 715].y =   0.5960792185; grid_points[ 715].z =   0.5960792185; grid_points[ 715].w_fixed =    0.1276124749; grid_points[ 715].w_total = 0.1276124734; grid_points[ 715].i_atom = 0;
    grid_points[ 716].x =   0.0000000000; grid_points[ 716].y =   1.3690014933; grid_points[ 716].z =   0.0000000000; grid_points[ 716].w_fixed =    0.4264811437; grid_points[ 716].w_total = 0.4264788573; grid_points[ 716].i_atom = 0;
    grid_points[ 717].x =   0.0000000000; grid_points[ 717].y =   0.9680302394; grid_points[ 717].z =   0.9680302394; grid_points[ 717].w_fixed =    0.3411849149; grid_points[ 717].w_total = 0.3411847747; grid_points[ 717].i_atom = 0;
    grid_points[ 718].x =  -0.9680302394; grid_points[ 718].y =   0.9680302394; grid_points[ 718].z =   0.0000000000; grid_points[ 718].w_fixed =    0.3411849149; grid_points[ 718].w_total = 0.3411846500; grid_points[ 718].i_atom = 0;
    grid_points[ 719].x =  -0.7903933807; grid_points[ 719].y =   0.7903933807; grid_points[ 719].z =   0.7903933807; grid_points[ 719].w_fixed =    0.2878747720; grid_points[ 719].w_total = 0.2878747232; grid_points[ 719].i_atom = 0;
    grid_points[ 720].x =   0.0000000000; grid_points[ 720].y =   1.8012554250; grid_points[ 720].z =   0.0000000000; grid_points[ 720].w_fixed =    0.2531886212; grid_points[ 720].w_total = 0.2531617203; grid_points[ 720].i_atom = 0;
    grid_points[ 721].x =   0.0000000000; grid_points[ 721].y =   1.2736799257; grid_points[ 721].z =   1.2736799257; grid_points[ 721].w_fixed =    0.4501131044; grid_points[ 721].w_total = 0.4501069125; grid_points[ 721].i_atom = 0;
    grid_points[ 722].x =  -1.0399553045; grid_points[ 722].y =   1.0399553045; grid_points[ 722].z =   1.0399553045; grid_points[ 722].w_fixed =    0.4205799265; grid_points[ 722].w_total = 0.4205787515; grid_points[ 722].i_atom = 0;
    grid_points[ 723].x =  -0.5430989451; grid_points[ 723].y =   0.5430989451; grid_points[ 723].z =   1.6292968354; grid_points[ 723].w_fixed =    0.4022281471; grid_points[ 723].w_total = 0.4021170368; grid_points[ 723].i_atom = 0;
    grid_points[ 724].x =  -0.5430989451; grid_points[ 724].y =   1.6292968354; grid_points[ 724].z =   0.5430989451; grid_points[ 724].w_fixed =    0.4022281471; grid_points[ 724].w_total = 0.4022123636; grid_points[ 724].i_atom = 0;
    grid_points[ 725].x =   0.0000000000; grid_points[ 725].y =   2.3583189787; grid_points[ 725].z =   0.0000000000; grid_points[ 725].w_fixed =    0.5606220456; grid_points[ 725].w_total = 0.5595902681; grid_points[ 725].i_atom = 0;
    grid_points[ 726].x =   0.0000000000; grid_points[ 726].y =   1.6675833421; grid_points[ 726].z =   1.6675833421; grid_points[ 726].w_fixed =    0.9966614143; grid_points[ 726].w_total = 0.9962426478; grid_points[ 726].i_atom = 0;
    grid_points[ 727].x =  -0.7110599262; grid_points[ 727].y =   0.7110599262; grid_points[ 727].z =   2.1331797787; grid_points[ 727].w_fixed =    0.8906323100; grid_points[ 727].w_total = 0.8798148112; grid_points[ 727].i_atom = 0;
    grid_points[ 728].x =  -0.7110599262; grid_points[ 728].y =   2.1331797787; grid_points[ 728].z =   0.7110599262; grid_points[ 728].w_fixed =    0.8906323100; grid_points[ 728].w_total = 0.8900464448; grid_points[ 728].i_atom = 0;
    grid_points[ 729].x =  -0.6652031738; grid_points[ 729].y =   2.1266717541; grid_points[ 729].z =   2.1266717541; grid_points[ 729].w_fixed =    0.9737495508; grid_points[ 729].w_total = 0.9690650946; grid_points[ 729].i_atom = 0;
    grid_points[ 730].x =  -1.2188238289; grid_points[ 730].y =   1.2188238289; grid_points[ 730].z =   2.5528215194; grid_points[ 730].w_fixed =    0.9397324606; grid_points[ 730].w_total = 0.8796756695; grid_points[ 730].i_atom = 0;
    grid_points[ 731].x =  -1.2188238289; grid_points[ 731].y =   2.5528215194; grid_points[ 731].z =   1.2188238289; grid_points[ 731].w_fixed =    0.9397324606; grid_points[ 731].w_total = 0.9361991824; grid_points[ 731].i_atom = 0;
    grid_points[ 732].x =   0.0000000000; grid_points[ 732].y =   1.4734978083; grid_points[ 732].z =   2.7049519350; grid_points[ 732].w_fixed =    0.9494794437; grid_points[ 732].w_total = 0.8169727814; grid_points[ 732].i_atom = 0;
    grid_points[ 733].x =   0.0000000000; grid_points[ 733].y =   2.7049519350; grid_points[ 733].z =   1.4734978083; grid_points[ 733].w_fixed =    0.9494794437; grid_points[ 733].w_total = 0.9435891428; grid_points[ 733].i_atom = 0;
    grid_points[ 734].x =  -0.5702030578; grid_points[ 734].y =   0.5702030578; grid_points[ 734].z =   2.9728265519; grid_points[ 734].w_fixed =    0.8042164697; grid_points[ 734].w_total = 0.4464214906; grid_points[ 734].i_atom = 0;
    grid_points[ 735].x =   0.0000000000; grid_points[ 735].y =   3.0802533600; grid_points[ 735].z =   0.0000000000; grid_points[ 735].w_fixed =    0.3749216611; grid_points[ 735].w_total = 0.3657764259; grid_points[ 735].i_atom = 0;
    grid_points[ 736].x =  -0.5702030578; grid_points[ 736].y =   2.9728265519; grid_points[ 736].z =   0.5702030578; grid_points[ 736].w_fixed =    0.8042164697; grid_points[ 736].w_total = 0.7903823761; grid_points[ 736].i_atom = 0;
    grid_points[ 737].x =  -0.8688367985; grid_points[ 737].y =   2.7776937196; grid_points[ 737].z =   2.7776937196; grid_points[ 737].w_fixed =    2.1793824633; grid_points[ 737].w_total = 2.1192111689; grid_points[ 737].i_atom = 0;
    grid_points[ 738].x =   0.0000000000; grid_points[ 738].y =  12.3210134399; grid_points[ 738].z =   0.0000000000; grid_points[ 738].w_fixed =   26.9943595963; grid_points[ 738].w_total = 0.0000000000; grid_points[ 738].i_atom = 0;
    grid_points[ 739].x =   0.0000000000; grid_points[ 739].y =  11.8583704325; grid_points[ 739].z =  11.8583704325; grid_points[ 739].w_fixed =  424.7625826047; grid_points[ 739].w_total = 0.0000000000; grid_points[ 739].i_atom = 0;
    grid_points[ 740].x =   6.3448636630; grid_points[ 740].y =  -6.3448636630; grid_points[ 740].z =   1.9846144277; grid_points[ 740].w_fixed =   27.8387029598; grid_points[ 740].w_total = 0.0012187121; grid_points[ 740].i_atom = 0;
    grid_points[ 741].x =   6.3448636630; grid_points[ 741].y =   1.9846144277; grid_points[ 741].z =  -6.3448636630; grid_points[ 741].w_fixed =   27.8387029598; grid_points[ 741].w_total = 0.0002691714; grid_points[ 741].i_atom = 0;
    grid_points[ 742].x =   6.6888597417; grid_points[ 742].y =   1.2829568800; grid_points[ 742].z =   1.2829568800; grid_points[ 742].w_fixed =    9.5422169009; grid_points[ 742].w_total = 0.0783957863; grid_points[ 742].i_atom = 0;
    grid_points[ 743].x =   7.1135404262; grid_points[ 743].y =  -7.1135404262; grid_points[ 743].z =  -7.1135404262; grid_points[ 743].w_fixed =   69.0587753800; grid_points[ 743].w_total = 0.0000000000; grid_points[ 743].i_atom = 0;
    grid_points[ 744].x =   8.5066870163; grid_points[ 744].y =  -8.5066870163; grid_points[ 744].z =   2.6608126953; grid_points[ 744].w_fixed =   70.1099676542; grid_points[ 744].w_total = 0.0000000000; grid_points[ 744].i_atom = 0;
    grid_points[ 745].x =   7.0235306998; grid_points[ 745].y =  -7.0235306998; grid_points[ 745].z =  21.0705920995; grid_points[ 745].w_fixed = 1092.3998745831; grid_points[ 745].w_total = 0.0000000000; grid_points[ 745].i_atom = 0;
    grid_points[ 746].x =   8.5066870163; grid_points[ 746].y =   2.6608126953; grid_points[ 746].z =  -8.5066870163; grid_points[ 746].w_fixed =   70.1099676542; grid_points[ 746].w_total = 0.0000000000; grid_points[ 746].i_atom = 0;
    grid_points[ 747].x =   6.9305700600; grid_points[ 747].y =   0.0000000000; grid_points[ 747].z =   0.0000000000; grid_points[ 747].w_fixed =    4.4485333807; grid_points[ 747].w_total = 0.0169320713; grid_points[ 747].i_atom = 0;
    grid_points[ 748].x =   9.1898468013; grid_points[ 748].y =   0.0000000000; grid_points[ 748].z =   0.0000000000; grid_points[ 748].w_fixed =   10.7187035387; grid_points[ 748].w_total = 0.0000000000; grid_points[ 748].i_atom = 0;
    grid_points[ 749].x =   8.8693420267; grid_points[ 749].y =   1.7011843293; grid_points[ 749].z =   1.7011843293; grid_points[ 749].w_fixed =   22.9918908793; grid_points[ 749].w_total = 0.0000018741; grid_points[ 749].i_atom = 0;
    grid_points[ 750].x =   7.0235306998; grid_points[ 750].y =  21.0705920995; grid_points[ 750].z =  -7.0235306998; grid_points[ 750].w_fixed = 1092.3998745831; grid_points[ 750].w_total = 0.0000000000; grid_points[ 750].i_atom = 0;
    grid_points[ 751].x =  16.4716395421; grid_points[ 751].y = -16.4716395421; grid_points[ 751].z =   0.0000000000; grid_points[ 751].w_fixed = 1222.4492551854; grid_points[ 751].w_total = 0.0000000000; grid_points[ 751].i_atom = 0;
    grid_points[ 752].x =  16.4716395421; grid_points[ 752].y =   0.0000000000; grid_points[ 752].z = -16.4716395421; grid_points[ 752].w_fixed = 1222.4492551854; grid_points[ 752].w_total = 0.0000000000; grid_points[ 752].i_atom = 0;
    grid_points[ 753].x =  16.7702682933; grid_points[ 753].y =   0.0000000000; grid_points[ 753].z =   0.0000000000; grid_points[ 753].w_fixed =  238.9289527151; grid_points[ 753].w_total = 0.0000000000; grid_points[ 753].i_atom = 0;
    grid_points[ 754].x = -10.8198077400; grid_points[ 754].y =   0.0000000000; grid_points[ 754].z =  -5.8939912332; grid_points[ 754].w_fixed =   68.3625199466; grid_points[ 754].w_total = 0.0000000000; grid_points[ 754].i_atom = 0;
    grid_points[ 755].x =  -2.7423536150; grid_points[ 755].y =   2.7423536150; grid_points[ 755].z =  -5.7438484187; grid_points[ 755].w_fixed =   11.1501458946; grid_points[ 755].w_total = 0.2040185179; grid_points[ 755].i_atom = 0;
    grid_points[ 756].x =  -1.5919331643; grid_points[ 756].y =   1.5919331643; grid_points[ 756].z =  -3.3342974947; grid_points[ 756].w_fixed =    2.1032476402; grid_points[ 756].w_total = 1.9796616798; grid_points[ 756].i_atom = 0;
    grid_points[ 757].x =  -1.9245685660; grid_points[ 757].y =   0.0000000000; grid_points[ 757].z =  -3.5329984457; grid_points[ 757].w_fixed =    2.1250626994; grid_points[ 757].w_total = 1.7554724438; grid_points[ 757].i_atom = 0;
    grid_points[ 758].x =  -2.0842608672; grid_points[ 758].y =   2.0842608672; grid_points[ 758].z =  -4.3654758527; grid_points[ 758].w_fixed =    4.7844044707; grid_points[ 758].w_total = 2.7579173487; grid_points[ 758].i_atom = 0;
    grid_points[ 759].x =  -2.5197684414; grid_points[ 759].y =   0.0000000000; grid_points[ 759].z =  -4.6256278652; grid_points[ 759].w_fixed =    4.8340287112; grid_points[ 759].w_total = 0.9672553950; grid_points[ 759].i_atom = 0;
    grid_points[ 760].x =  -1.5919331643; grid_points[ 760].y =   1.5919331643; grid_points[ 760].z =   3.3342974947; grid_points[ 760].w_fixed =    2.1032476402; grid_points[ 760].w_total = 1.2794689390; grid_points[ 760].i_atom = 0;
    grid_points[ 761].x =  -2.0842608672; grid_points[ 761].y =   2.0842608672; grid_points[ 761].z =   4.3654758527; grid_points[ 761].w_fixed =    4.7844044707; grid_points[ 761].w_total = 1.0736736243; grid_points[ 761].i_atom = 0;
    grid_points[ 762].x =  -2.7423536150; grid_points[ 762].y =   2.7423536150; grid_points[ 762].z =   5.7438484187; grid_points[ 762].w_fixed =   11.1501458946; grid_points[ 762].w_total = 0.0554063878; grid_points[ 762].i_atom = 0;
    grid_points[ 763].x =  -0.7447550142; grid_points[ 763].y =   0.7447550142; grid_points[ 763].z =  -3.8828754963; grid_points[ 763].w_fixed =    1.7999446257; grid_points[ 763].w_total = 0.6776187779; grid_points[ 763].i_atom = 0;
    grid_points[ 764].x =  -0.9750809686; grid_points[ 764].y =   0.9750809686; grid_points[ 764].z =  -5.0837093106; grid_points[ 764].w_fixed =    4.0944598961; grid_points[ 764].w_total = 0.0118697638; grid_points[ 764].i_atom = 0;
    grid_points[ 765].x =  -0.7447550142; grid_points[ 765].y =   0.7447550142; grid_points[ 765].z =   3.8828754963; grid_points[ 765].w_fixed =    1.7999446257; grid_points[ 765].w_total = 0.0335908066; grid_points[ 765].i_atom = 0;
    grid_points[ 766].x =   0.0000000000; grid_points[ 766].y =   1.9245685660; grid_points[ 766].z =   3.5329984457; grid_points[ 766].w_fixed =    2.1250626994; grid_points[ 766].w_total = 0.7185910872; grid_points[ 766].i_atom = 0;
    grid_points[ 767].x =  -0.9750809686; grid_points[ 767].y =   0.9750809686; grid_points[ 767].z =   5.0837093106; grid_points[ 767].w_fixed =    4.0944598961; grid_points[ 767].w_total = 0.0001490129; grid_points[ 767].i_atom = 0;
    grid_points[ 768].x =   0.0000000000; grid_points[ 768].y =   2.5197684414; grid_points[ 768].z =   4.6256278652; grid_points[ 768].w_fixed =    4.8340287112; grid_points[ 768].w_total = 0.2039773991; grid_points[ 768].i_atom = 0;
    grid_points[ 769].x =   0.0000000000; grid_points[ 769].y =  10.8198077400; grid_points[ 769].z =  -5.8939912332; grid_points[ 769].w_fixed =   68.3625199466; grid_points[ 769].w_total = 0.0000000000; grid_points[ 769].i_atom = 0;
    grid_points[ 770].x =   0.0000000000; grid_points[ 770].y =  10.8198077400; grid_points[ 770].z =   5.8939912332; grid_points[ 770].w_fixed =   68.3625199466; grid_points[ 770].w_total = 0.0000000000; grid_points[ 770].i_atom = 0;
    grid_points[ 771].x =   8.0701458556; grid_points[ 771].y =   0.0000000000; grid_points[ 771].z =  -4.3961380893; grid_points[ 771].w_fixed =   27.1448404564; grid_points[ 771].w_total = 0.0000094033; grid_points[ 771].i_atom = 0;
    grid_points[ 772].x = -30.0795229432; grid_points[ 772].y =  10.0265076477; grid_points[ 772].z = -10.0265076477; grid_points[ 772].w_fixed = 3474.1698452873; grid_points[ 772].w_total = 0.0000000000; grid_points[ 772].i_atom = 0;
    grid_points[ 773].x = -21.0705920995; grid_points[ 773].y =  -7.0235306998; grid_points[ 773].z =   7.0235306998; grid_points[ 773].w_fixed = 1092.3998745831; grid_points[ 773].w_total = 0.0000000000; grid_points[ 773].i_atom = 0;
    grid_points[ 774].x = -19.1993004361; grid_points[ 774].y =  19.1993004361; grid_points[ 774].z = -19.1993004361; grid_points[ 774].w_fixed = 3632.6798826243; grid_points[ 774].w_total = 0.0000000000; grid_points[ 774].i_atom = 0;
    grid_points[ 775].x = -11.8583704325; grid_points[ 775].y =   0.0000000000; grid_points[ 775].z = -11.8583704325; grid_points[ 775].w_fixed =  424.7625826047; grid_points[ 775].w_total = 0.0000000000; grid_points[ 775].i_atom = 0;
    grid_points[ 776].x = -11.8913062075; grid_points[ 776].y =   2.2808122310; grid_points[ 776].z =  -2.2808122310; grid_points[ 776].w_fixed =   57.9035858166; grid_points[ 776].w_total = 0.0000000000; grid_points[ 776].i_atom = 0;
    grid_points[ 777].x = -10.0265076477; grid_points[ 777].y =  10.0265076477; grid_points[ 777].z = -30.0795229432; grid_points[ 777].w_fixed = 3474.1698452873; grid_points[ 777].w_total = 0.0000000000; grid_points[ 777].i_atom = 0;
    grid_points[ 778].x =  -9.6823189135; grid_points[ 778].y =   9.6823189135; grid_points[ 778].z =  -9.6823189135; grid_points[ 778].w_fixed =  396.8927232114; grid_points[ 778].w_total = 0.0000000000; grid_points[ 778].i_atom = 0;
    grid_points[ 779].x = -10.0265076477; grid_points[ 779].y =  30.0795229432; grid_points[ 779].z = -10.0265076477; grid_points[ 779].w_fixed = 3474.1698452873; grid_points[ 779].w_total = 0.0000000000; grid_points[ 779].i_atom = 0;
    grid_points[ 780].x =  -1.2829568800; grid_points[ 780].y =  -6.6888597417; grid_points[ 780].z =  -1.2829568800; grid_points[ 780].w_fixed =    9.5422169009; grid_points[ 780].w_total = 0.0865328484; grid_points[ 780].i_atom = 0;
    grid_points[ 781].x =  -1.7011843293; grid_points[ 781].y =  -8.8693420267; grid_points[ 781].z =  -1.7011843293; grid_points[ 781].w_fixed =   22.9918908793; grid_points[ 781].w_total = 0.0000042871; grid_points[ 781].i_atom = 0;
    grid_points[ 782].x =  -1.9846144277; grid_points[ 782].y =  -6.3448636630; grid_points[ 782].z =   6.3448636630; grid_points[ 782].w_fixed =   27.8387029598; grid_points[ 782].w_total = 0.0001758536; grid_points[ 782].i_atom = 0;
    grid_points[ 783].x =  -2.6608126953; grid_points[ 783].y =  -8.5066870163; grid_points[ 783].z =   8.5066870163; grid_points[ 783].w_fixed =   70.1099676542; grid_points[ 783].w_total = 0.0000000000; grid_points[ 783].i_atom = 0;
    grid_points[ 784].x =  -2.2808122310; grid_points[ 784].y =   2.2808122310; grid_points[ 784].z = -11.8913062075; grid_points[ 784].w_fixed =   57.9035858166; grid_points[ 784].w_total = 0.0000000000; grid_points[ 784].i_atom = 0;
    grid_points[ 785].x =  -1.2736799257; grid_points[ 785].y =   0.0000000000; grid_points[ 785].z =  -1.2736799257; grid_points[ 785].w_fixed =    0.4501131044; grid_points[ 785].w_total = 0.4501105971; grid_points[ 785].i_atom = 0;
    grid_points[ 786].x =  -1.6675833421; grid_points[ 786].y =   0.0000000000; grid_points[ 786].z =  -1.6675833421; grid_points[ 786].w_fixed =    0.9966614143; grid_points[ 786].w_total = 0.9965883126; grid_points[ 786].i_atom = 0;
    grid_points[ 787].x =  -1.3615760972; grid_points[ 787].y =   1.3615760972; grid_points[ 787].z =  -1.3615760972; grid_points[ 787].w_fixed =    0.9312676753; grid_points[ 787].w_total = 0.9312379628; grid_points[ 787].i_atom = 0;
    grid_points[ 788].x =  -1.7783851066; grid_points[ 788].y =   1.7783851066; grid_points[ 788].z =  -1.7783851066; grid_points[ 788].w_fixed =    0.9591496581; grid_points[ 788].w_total = 0.9588589425; grid_points[ 788].i_atom = 0;
    grid_points[ 789].x =  -2.1266717541; grid_points[ 789].y =   0.6652031738; grid_points[ 789].z =  -2.1266717541; grid_points[ 789].w_fixed =    0.9737495508; grid_points[ 789].w_total = 0.9730669536; grid_points[ 789].i_atom = 0;
    grid_points[ 790].x =  -1.4734978083; grid_points[ 790].y =   0.0000000000; grid_points[ 790].z =  -2.7049519350; grid_points[ 790].w_fixed =    0.9494794437; grid_points[ 790].w_total = 0.9420711494; grid_points[ 790].i_atom = 0;
    grid_points[ 791].x =  -2.7049519350; grid_points[ 791].y =   0.0000000000; grid_points[ 791].z =  -1.4734978083; grid_points[ 791].w_fixed =    0.9494794437; grid_points[ 791].w_total = 0.9427515891; grid_points[ 791].i_atom = 0;
    grid_points[ 792].x =  -2.3227887106; grid_points[ 792].y =   2.3227887106; grid_points[ 792].z =  -2.3227887106; grid_points[ 792].w_fixed =    2.1467059397; grid_points[ 792].w_total = 2.1417068903; grid_points[ 792].i_atom = 0;
    grid_points[ 793].x =  -2.7776937196; grid_points[ 793].y =   0.8688367985; grid_points[ 793].z =  -2.7776937196; grid_points[ 793].w_fixed =    2.1793824633; grid_points[ 793].w_total = 2.1646215031; grid_points[ 793].i_atom = 0;
    grid_points[ 794].x =  -1.6292968354; grid_points[ 794].y =   0.5430989451; grid_points[ 794].z =  -0.5430989451; grid_points[ 794].w_fixed =    0.4022281471; grid_points[ 794].w_total = 0.4022099166; grid_points[ 794].i_atom = 0;
    grid_points[ 795].x =  -2.1331797787; grid_points[ 795].y =   0.7110599262; grid_points[ 795].z =  -0.7110599262; grid_points[ 795].w_fixed =    0.8906323100; grid_points[ 795].w_total = 0.8899740741; grid_points[ 795].i_atom = 0;
    grid_points[ 796].x =  -2.9728265519; grid_points[ 796].y =   0.5702030578; grid_points[ 796].z =  -0.5702030578; grid_points[ 796].w_fixed =    0.8042164697; grid_points[ 796].w_total = 0.7899201355; grid_points[ 796].i_atom = 0;
    grid_points[ 797].x =  -2.1266717541; grid_points[ 797].y =   2.1266717541; grid_points[ 797].z =  -0.6652031738; grid_points[ 797].w_fixed =    0.9737495508; grid_points[ 797].w_total = 0.9721447960; grid_points[ 797].i_atom = 0;
    grid_points[ 798].x =  -2.5528215194; grid_points[ 798].y =   1.2188238289; grid_points[ 798].z =  -1.2188238289; grid_points[ 798].w_fixed =    0.9397324606; grid_points[ 798].w_total = 0.9357421258; grid_points[ 798].i_atom = 0;
    grid_points[ 799].x =  -1.2829568800; grid_points[ 799].y =   1.2829568800; grid_points[ 799].z =   6.6888597417; grid_points[ 799].w_fixed =    9.5422169009; grid_points[ 799].w_total = 0.0000003249; grid_points[ 799].i_atom = 0;
    grid_points[ 800].x =  -1.7011843293; grid_points[ 800].y =   1.7011843293; grid_points[ 800].z =   8.8693420267; grid_points[ 800].w_fixed =   22.9918908793; grid_points[ 800].w_total = 0.0000000000; grid_points[ 800].i_atom = 0;
    grid_points[ 801].x =  -3.0411437621; grid_points[ 801].y =   3.0411437621; grid_points[ 801].z =  -3.0411437621; grid_points[ 801].w_fixed =    4.8832621034; grid_points[ 801].w_total = 4.8017896257; grid_points[ 801].i_atom = 0;
    grid_points[ 802].x =  -2.7776937196; grid_points[ 802].y =   2.7776937196; grid_points[ 802].z =  -0.8688367985; grid_points[ 802].w_fixed =    2.1793824633; grid_points[ 802].w_total = 2.1466862309; grid_points[ 802].i_atom = 0;
    grid_points[ 803].x =  -2.2808122310; grid_points[ 803].y =  11.8913062075; grid_points[ 803].z =  -2.2808122310; grid_points[ 803].w_fixed =   57.9035858166; grid_points[ 803].w_total = 0.0000000000; grid_points[ 803].i_atom = 0;
    grid_points[ 804].x =  -0.5430989451; grid_points[ 804].y =   0.5430989451; grid_points[ 804].z =  -1.6292968354; grid_points[ 804].w_fixed =    0.4022281471; grid_points[ 804].w_total = 0.4022275009; grid_points[ 804].i_atom = 0;
    grid_points[ 805].x =  -0.7110599262; grid_points[ 805].y =   0.7110599262; grid_points[ 805].z =  -2.1331797787; grid_points[ 805].w_fixed =    0.8906323100; grid_points[ 805].w_total = 0.8904412859; grid_points[ 805].i_atom = 0;
    grid_points[ 806].x =  -0.5702030578; grid_points[ 806].y =   0.5702030578; grid_points[ 806].z =  -2.9728265519; grid_points[ 806].w_fixed =    0.8042164697; grid_points[ 806].w_total = 0.7676826133; grid_points[ 806].i_atom = 0;
    grid_points[ 807].x =  -0.6652031738; grid_points[ 807].y =   2.1266717541; grid_points[ 807].z =  -2.1266717541; grid_points[ 807].w_fixed =    0.9737495508; grid_points[ 807].w_total = 0.9730669536; grid_points[ 807].i_atom = 0;
    grid_points[ 808].x =  -1.2188238289; grid_points[ 808].y =   1.2188238289; grid_points[ 808].z =  -2.5528215194; grid_points[ 808].w_fixed =    0.9397324606; grid_points[ 808].w_total = 0.9373660293; grid_points[ 808].i_atom = 0;
    grid_points[ 809].x =   0.0000000000; grid_points[ 809].y =   2.7049519350; grid_points[ 809].z =  -1.4734978083; grid_points[ 809].w_fixed =    0.9494794437; grid_points[ 809].w_total = 0.9427515891; grid_points[ 809].i_atom = 0;
    grid_points[ 810].x =  -0.4000533132; grid_points[ 810].y =   0.0000000000; grid_points[ 810].z =  -0.4000533132; grid_points[ 810].w_fixed =    0.0275212788; grid_points[ 810].w_total = 0.0275212788; grid_points[ 810].i_atom = 0;
    grid_points[ 811].x =  -0.3266421624; grid_points[ 811].y =   0.3266421624; grid_points[ 811].z =  -0.3266421624; grid_points[ 811].w_fixed =    0.0232210790; grid_points[ 811].w_total = 0.0232210790; grid_points[ 811].i_atom = 0;
    grid_points[ 812].x =  -0.5445170097; grid_points[ 812].y =   0.0000000000; grid_points[ 812].z =  -0.5445170097; grid_points[ 812].w_fixed =    0.0655814715; grid_points[ 812].w_total = 0.0655814714; grid_points[ 812].i_atom = 0;
    grid_points[ 813].x =  -0.4445962766; grid_points[ 813].y =   0.4445962766; grid_points[ 813].z =  -0.4445962766; grid_points[ 813].w_fixed =    0.0553343666; grid_points[ 813].w_total = 0.0553343665; grid_points[ 813].i_atom = 0;
    grid_points[ 814].x =  -0.7300449659; grid_points[ 814].y =   0.0000000000; grid_points[ 814].z =  -0.7300449659; grid_points[ 814].w_fixed =    0.1512444148; grid_points[ 814].w_total = 0.1512444110; grid_points[ 814].i_atom = 0;
    grid_points[ 815].x =  -0.5960792185; grid_points[ 815].y =   0.5960792185; grid_points[ 815].z =  -0.5960792185; grid_points[ 815].w_fixed =    0.1276124749; grid_points[ 815].w_total = 0.1276124729; grid_points[ 815].i_atom = 0;
    grid_points[ 816].x =  -0.9680302394; grid_points[ 816].y =   0.0000000000; grid_points[ 816].z =  -0.9680302394; grid_points[ 816].w_fixed =    0.3411849149; grid_points[ 816].w_total = 0.3411847831; grid_points[ 816].i_atom = 0;
    grid_points[ 817].x =  -0.7903933807; grid_points[ 817].y =   0.7903933807; grid_points[ 817].z =  -0.7903933807; grid_points[ 817].w_fixed =    0.2878747720; grid_points[ 817].w_total = 0.2878747082; grid_points[ 817].i_atom = 0;
    grid_points[ 818].x =  -1.0399553045; grid_points[ 818].y =   1.0399553045; grid_points[ 818].z =  -1.0399553045; grid_points[ 818].w_fixed =    0.4205799265; grid_points[ 818].w_total = 0.4205787392; grid_points[ 818].i_atom = 0;
    grid_points[ 819].x =  -0.5430989451; grid_points[ 819].y =   1.6292968354; grid_points[ 819].z =  -0.5430989451; grid_points[ 819].w_fixed =    0.4022281471; grid_points[ 819].w_total = 0.4022099166; grid_points[ 819].i_atom = 0;
    grid_points[ 820].x =  -0.7110599262; grid_points[ 820].y =   2.1331797787; grid_points[ 820].z =  -0.7110599262; grid_points[ 820].w_fixed =    0.8906323100; grid_points[ 820].w_total = 0.8899740741; grid_points[ 820].i_atom = 0;
    grid_points[ 821].x =  -1.2188238289; grid_points[ 821].y =   2.5528215194; grid_points[ 821].z =  -1.2188238289; grid_points[ 821].w_fixed =    0.9397324606; grid_points[ 821].w_total = 0.9357421258; grid_points[ 821].i_atom = 0;
    grid_points[ 822].x =  -0.8688367985; grid_points[ 822].y =   2.7776937196; grid_points[ 822].z =  -2.7776937196; grid_points[ 822].w_fixed =    2.1793824633; grid_points[ 822].w_total = 2.1646215031; grid_points[ 822].i_atom = 0;
    grid_points[ 823].x =  -0.5702030578; grid_points[ 823].y =   2.9728265519; grid_points[ 823].z =  -0.5702030578; grid_points[ 823].w_fixed =    0.8042164697; grid_points[ 823].w_total = 0.7899201355; grid_points[ 823].i_atom = 0;
    grid_points[ 824].x =   6.3448636630; grid_points[ 824].y =  -6.3448636630; grid_points[ 824].z =  -1.9846144277; grid_points[ 824].w_fixed =   27.8387029598; grid_points[ 824].w_total = 0.0013808480; grid_points[ 824].i_atom = 0;
    grid_points[ 825].x =   6.6888597417; grid_points[ 825].y =   1.2829568800; grid_points[ 825].z =  -1.2829568800; grid_points[ 825].w_fixed =    9.5422169009; grid_points[ 825].w_total = 0.0865328484; grid_points[ 825].i_atom = 0;
    grid_points[ 826].x =   6.3448636630; grid_points[ 826].y =   1.9846144277; grid_points[ 826].z =   6.3448636630; grid_points[ 826].w_fixed =   27.8387029598; grid_points[ 826].w_total = 0.0001758536; grid_points[ 826].i_atom = 0;
    grid_points[ 827].x =   7.0235306998; grid_points[ 827].y =  -7.0235306998; grid_points[ 827].z = -21.0705920995; grid_points[ 827].w_fixed = 1092.3998745831; grid_points[ 827].w_total = 0.0000000000; grid_points[ 827].i_atom = 0;
    grid_points[ 828].x =   8.5066870163; grid_points[ 828].y =  -8.5066870163; grid_points[ 828].z =  -2.6608126953; grid_points[ 828].w_fixed =   70.1099676542; grid_points[ 828].w_total = 0.0000000000; grid_points[ 828].i_atom = 0;
    grid_points[ 829].x =   7.1135404262; grid_points[ 829].y =  -7.1135404262; grid_points[ 829].z =   7.1135404262; grid_points[ 829].w_fixed =   69.0587753800; grid_points[ 829].w_total = 0.0000000000; grid_points[ 829].i_atom = 0;
    grid_points[ 830].x =   8.8693420267; grid_points[ 830].y =   1.7011843293; grid_points[ 830].z =  -1.7011843293; grid_points[ 830].w_fixed =   22.9918908793; grid_points[ 830].w_total = 0.0000042871; grid_points[ 830].i_atom = 0;
    grid_points[ 831].x =   8.5066870163; grid_points[ 831].y =   2.6608126953; grid_points[ 831].z =   8.5066870163; grid_points[ 831].w_fixed =   70.1099676542; grid_points[ 831].w_total = 0.0000000000; grid_points[ 831].i_atom = 0;
    grid_points[ 832].x =   7.0235306998; grid_points[ 832].y =  21.0705920995; grid_points[ 832].z =   7.0235306998; grid_points[ 832].w_fixed = 1092.3998745831; grid_points[ 832].w_total = 0.0000000000; grid_points[ 832].i_atom = 0;
    grid_points[ 833].x = -10.8198077400; grid_points[ 833].y =  -5.8939912332; grid_points[ 833].z =   0.0000000000; grid_points[ 833].w_fixed =   68.3625199466; grid_points[ 833].w_total = 0.0000000000; grid_points[ 833].i_atom = 0;
    grid_points[ 834].x = -10.8198077400; grid_points[ 834].y =   5.8939912332; grid_points[ 834].z =   0.0000000000; grid_points[ 834].w_fixed =   68.3625199466; grid_points[ 834].w_total = 0.0000000000; grid_points[ 834].i_atom = 0;
    grid_points[ 835].x =  -2.7423536150; grid_points[ 835].y =  -5.7438484187; grid_points[ 835].z =   2.7423536150; grid_points[ 835].w_fixed =   11.1501458946; grid_points[ 835].w_total = 0.7753792790; grid_points[ 835].i_atom = 0;
    grid_points[ 836].x =  -1.5919331643; grid_points[ 836].y =  -3.3342974947; grid_points[ 836].z =   1.5919331643; grid_points[ 836].w_fixed =    2.1032476402; grid_points[ 836].w_total = 2.0194520110; grid_points[ 836].i_atom = 0;
    grid_points[ 837].x =  -1.9245685660; grid_points[ 837].y =  -3.5329984457; grid_points[ 837].z =   0.0000000000; grid_points[ 837].w_fixed =    2.1250626994; grid_points[ 837].w_total = 1.9718171965; grid_points[ 837].i_atom = 0;
    grid_points[ 838].x =  -2.0842608672; grid_points[ 838].y =  -4.3654758527; grid_points[ 838].z =   2.0842608672; grid_points[ 838].w_fixed =    4.7844044707; grid_points[ 838].w_total = 3.3058550217; grid_points[ 838].i_atom = 0;
    grid_points[ 839].x =  -2.5197684414; grid_points[ 839].y =  -4.6256278652; grid_points[ 839].z =   0.0000000000; grid_points[ 839].w_fixed =    4.8340287112; grid_points[ 839].w_total = 2.6697241981; grid_points[ 839].i_atom = 0;
    grid_points[ 840].x =  -1.5919331643; grid_points[ 840].y =   3.3342974947; grid_points[ 840].z =   1.5919331643; grid_points[ 840].w_fixed =    2.1032476402; grid_points[ 840].w_total = 2.0194520110; grid_points[ 840].i_atom = 0;
    grid_points[ 841].x =  -1.9245685660; grid_points[ 841].y =   3.5329984457; grid_points[ 841].z =   0.0000000000; grid_points[ 841].w_fixed =    2.1250626994; grid_points[ 841].w_total = 1.9718171965; grid_points[ 841].i_atom = 0;
    grid_points[ 842].x =  -2.0842608672; grid_points[ 842].y =   4.3654758527; grid_points[ 842].z =   2.0842608672; grid_points[ 842].w_fixed =    4.7844044707; grid_points[ 842].w_total = 3.3058550217; grid_points[ 842].i_atom = 0;
    grid_points[ 843].x =  -2.5197684414; grid_points[ 843].y =   4.6256278652; grid_points[ 843].z =   0.0000000000; grid_points[ 843].w_fixed =    4.8340287112; grid_points[ 843].w_total = 2.6697241981; grid_points[ 843].i_atom = 0;
    grid_points[ 844].x =  -2.7423536150; grid_points[ 844].y =   5.7438484187; grid_points[ 844].z =   2.7423536150; grid_points[ 844].w_fixed =   11.1501458946; grid_points[ 844].w_total = 0.7753792790; grid_points[ 844].i_atom = 0;
    grid_points[ 845].x =   0.0000000000; grid_points[ 845].y =  -5.8939912332; grid_points[ 845].z =  10.8198077400; grid_points[ 845].w_fixed =   68.3625199466; grid_points[ 845].w_total = 0.0000000000; grid_points[ 845].i_atom = 0;
    grid_points[ 846].x =  -0.7447550142; grid_points[ 846].y =  -3.8828754963; grid_points[ 846].z =   0.7447550142; grid_points[ 846].w_fixed =    1.7999446257; grid_points[ 846].w_total = 1.5204246946; grid_points[ 846].i_atom = 0;
    grid_points[ 847].x =  -0.9750809686; grid_points[ 847].y =  -5.0837093106; grid_points[ 847].z =   0.9750809686; grid_points[ 847].w_fixed =    4.0944598961; grid_points[ 847].w_total = 1.3604510062; grid_points[ 847].i_atom = 0;
    grid_points[ 848].x =   0.0000000000; grid_points[ 848].y =   4.0231880620; grid_points[ 848].z =   0.0000000000; grid_points[ 848].w_fixed =    0.8391251042; grid_points[ 848].w_total = 0.6689355801; grid_points[ 848].i_atom = 0;
    grid_points[ 849].x =  -0.7447550142; grid_points[ 849].y =   3.8828754963; grid_points[ 849].z =   0.7447550142; grid_points[ 849].w_fixed =    1.7999446257; grid_points[ 849].w_total = 1.5204246946; grid_points[ 849].i_atom = 0;
    grid_points[ 850].x =   0.0000000000; grid_points[ 850].y =   3.5329984457; grid_points[ 850].z =   1.9245685660; grid_points[ 850].w_fixed =    2.1250626994; grid_points[ 850].w_total = 1.9844543174; grid_points[ 850].i_atom = 0;
    grid_points[ 851].x =   0.0000000000; grid_points[ 851].y =   5.2674155091; grid_points[ 851].z =   0.0000000000; grid_points[ 851].w_fixed =    1.9088165479; grid_points[ 851].w_total = 0.4997305992; grid_points[ 851].i_atom = 0;
    grid_points[ 852].x =  -0.9750809686; grid_points[ 852].y =   5.0837093106; grid_points[ 852].z =   0.9750809686; grid_points[ 852].w_fixed =    4.0944598961; grid_points[ 852].w_total = 1.3604510062; grid_points[ 852].i_atom = 0;
    grid_points[ 853].x =   0.0000000000; grid_points[ 853].y =   4.6256278652; grid_points[ 853].z =   2.5197684414; grid_points[ 853].w_fixed =    4.8340287112; grid_points[ 853].w_total = 2.6802181871; grid_points[ 853].i_atom = 0;
    grid_points[ 854].x =   0.0000000000; grid_points[ 854].y =   5.8939912332; grid_points[ 854].z =  10.8198077400; grid_points[ 854].w_fixed =   68.3625199466; grid_points[ 854].w_total = 0.0000000000; grid_points[ 854].i_atom = 0;
    grid_points[ 855].x =   0.0000000000; grid_points[ 855].y =  23.2944160349; grid_points[ 855].z =   0.0000000000; grid_points[ 855].w_fixed =  687.6277060418; grid_points[ 855].w_total = 0.0000000000; grid_points[ 855].i_atom = 0;
    grid_points[ 856].x =   0.0000000000; grid_points[ 856].y =  33.2541638252; grid_points[ 856].z =   0.0000000000; grid_points[ 856].w_fixed = 2186.8690181113; grid_points[ 856].w_total = 0.0000000000; grid_points[ 856].i_atom = 0;
    grid_points[ 857].x =   8.0701458556; grid_points[ 857].y =  -4.3961380893; grid_points[ 857].z =   0.0000000000; grid_points[ 857].w_fixed =   27.1448404564; grid_points[ 857].w_total = 0.0003003544; grid_points[ 857].i_atom = 0;
    grid_points[ 858].x =   8.0701458556; grid_points[ 858].y =   4.3961380893; grid_points[ 858].z =   0.0000000000; grid_points[ 858].w_fixed =   27.1448404564; grid_points[ 858].w_total = 0.0003003544; grid_points[ 858].i_atom = 0;
    grid_points[ 859].x = -10.2112860776; grid_points[ 859].y =  -4.8752953156; grid_points[ 859].z =  -4.8752953156; grid_points[ 859].w_fixed =   67.6607371617; grid_points[ 859].w_total = 0.0000000000; grid_points[ 859].i_atom = 0;
    grid_points[ 860].x = -10.2112860776; grid_points[ 860].y =  -4.8752953156; grid_points[ 860].z =   4.8752953156; grid_points[ 860].w_fixed =   67.6607371617; grid_points[ 860].w_total = 0.0000000000; grid_points[ 860].i_atom = 0;
    grid_points[ 861].x = -10.2112860776; grid_points[ 861].y =   4.8752953156; grid_points[ 861].z =  -4.8752953156; grid_points[ 861].w_fixed =   67.6607371617; grid_points[ 861].w_total = 0.0000000000; grid_points[ 861].i_atom = 0;
    grid_points[ 862].x = -10.2112860776; grid_points[ 862].y =   4.8752953156; grid_points[ 862].z =   4.8752953156; grid_points[ 862].w_fixed =   67.6607371617; grid_points[ 862].w_total = 0.0000000000; grid_points[ 862].i_atom = 0;
    grid_points[ 863].x =  -1.4967071411; grid_points[ 863].y =  -4.7850114467; grid_points[ 863].z =  -4.7850114467; grid_points[ 863].w_fixed =   11.5537666422; grid_points[ 863].w_total = 2.9255841962; grid_points[ 863].i_atom = 0;
    grid_points[ 864].x =  -1.4967071411; grid_points[ 864].y =  -4.7850114467; grid_points[ 864].z =   4.7850114467; grid_points[ 864].w_fixed =   11.5537666422; grid_points[ 864].w_total = 1.8853671084; grid_points[ 864].i_atom = 0;
    grid_points[ 865].x =  -1.4967071411; grid_points[ 865].y =   4.7850114467; grid_points[ 865].z =  -4.7850114467; grid_points[ 865].w_fixed =   11.5537666422; grid_points[ 865].w_total = 2.9255841962; grid_points[ 865].i_atom = 0;
    grid_points[ 866].x =  -1.4967071411; grid_points[ 866].y =   4.7850114467; grid_points[ 866].z =   4.7850114467; grid_points[ 866].w_fixed =   11.5537666422; grid_points[ 866].w_total = 1.8853671084; grid_points[ 866].i_atom = 0;
    grid_points[ 867].x =  -1.1375367884; grid_points[ 867].y =  -3.6367345381; grid_points[ 867].z =  -3.6367345381; grid_points[ 867].w_fixed =    4.9575936763; grid_points[ 867].w_total = 4.6559094346; grid_points[ 867].i_atom = 0;
    grid_points[ 868].x =  -1.1375367884; grid_points[ 868].y =  -3.6367345381; grid_points[ 868].z =   3.6367345381; grid_points[ 868].w_fixed =    4.9575936763; grid_points[ 868].w_total = 4.5221022069; grid_points[ 868].i_atom = 0;
    grid_points[ 869].x =   0.0000000000; grid_points[ 869].y =  -3.3153700687; grid_points[ 869].z =   6.0861418537; grid_points[ 869].w_fixed =   11.2657961338; grid_points[ 869].w_total = 0.0045977423; grid_points[ 869].i_atom = 0;
    grid_points[ 870].x =  -1.1375367884; grid_points[ 870].y =   3.6367345381; grid_points[ 870].z =  -3.6367345381; grid_points[ 870].w_fixed =    4.9575936763; grid_points[ 870].w_total = 4.6559094346; grid_points[ 870].i_atom = 0;
    grid_points[ 871].x =   0.0000000000; grid_points[ 871].y =   6.0861418537; grid_points[ 871].z =  -3.3153700687; grid_points[ 871].w_fixed =   11.2657961338; grid_points[ 871].w_total = 0.3549590929; grid_points[ 871].i_atom = 0;
    grid_points[ 872].x =  -1.1375367884; grid_points[ 872].y =   3.6367345381; grid_points[ 872].z =   3.6367345381; grid_points[ 872].w_fixed =    4.9575936763; grid_points[ 872].w_total = 4.5221022069; grid_points[ 872].i_atom = 0;
    grid_points[ 873].x =   0.0000000000; grid_points[ 873].y =   3.3153700687; grid_points[ 873].z =   6.0861418537; grid_points[ 873].w_fixed =   11.2657961338; grid_points[ 873].w_total = 0.0045977423; grid_points[ 873].i_atom = 0;
    grid_points[ 874].x =   0.0000000000; grid_points[ 874].y =   6.0861418537; grid_points[ 874].z =   3.3153700687; grid_points[ 874].w_fixed =   11.2657961338; grid_points[ 874].w_total = 0.3471224190; grid_points[ 874].i_atom = 0;
    grid_points[ 875].x =   0.0000000000; grid_points[ 875].y =  23.5142447435; grid_points[ 875].z =  23.5142447435; grid_points[ 875].w_fixed = 3887.7671433090; grid_points[ 875].w_total = 0.0000000000; grid_points[ 875].i_atom = 0;
    grid_points[ 876].x =   7.6162691612; grid_points[ 876].y =  -3.6363256383; grid_points[ 876].z =  -3.6363256383; grid_points[ 876].w_fixed =   26.8661821836; grid_points[ 876].w_total = 0.0010773582; grid_points[ 876].i_atom = 0;
    grid_points[ 877].x =   7.6162691612; grid_points[ 877].y =  -3.6363256383; grid_points[ 877].z =   3.6363256383; grid_points[ 877].w_fixed =   26.8661821836; grid_points[ 877].w_total = 0.0009330210; grid_points[ 877].i_atom = 0;
    grid_points[ 878].x =   7.6162691612; grid_points[ 878].y =   3.6363256383; grid_points[ 878].z =  -3.6363256383; grid_points[ 878].w_fixed =   26.8661821836; grid_points[ 878].w_total = 0.0010773582; grid_points[ 878].i_atom = 0;
    grid_points[ 879].x =   7.6162691612; grid_points[ 879].y =   3.6363256383; grid_points[ 879].z =   3.6363256383; grid_points[ 879].w_fixed =   26.8661821836; grid_points[ 879].w_total = 0.0009330210; grid_points[ 879].i_atom = 0;
    grid_points[ 880].x =  -2.7423536150; grid_points[ 880].y =  -5.7438484187; grid_points[ 880].z =  -2.7423536150; grid_points[ 880].w_fixed =   11.1501458946; grid_points[ 880].w_total = 0.7975816561; grid_points[ 880].i_atom = 0;
    grid_points[ 881].x =  -1.5919331643; grid_points[ 881].y =  -3.3342974947; grid_points[ 881].z =  -1.5919331643; grid_points[ 881].w_fixed =    2.1032476402; grid_points[ 881].w_total = 2.0149126278; grid_points[ 881].i_atom = 0;
    grid_points[ 882].x =  -2.0842608672; grid_points[ 882].y =  -4.3654758527; grid_points[ 882].z =  -2.0842608672; grid_points[ 882].w_fixed =    4.7844044707; grid_points[ 882].w_total = 3.2966811514; grid_points[ 882].i_atom = 0;
    grid_points[ 883].x =  -1.5919331643; grid_points[ 883].y =   3.3342974947; grid_points[ 883].z =  -1.5919331643; grid_points[ 883].w_fixed =    2.1032476402; grid_points[ 883].w_total = 2.0149126278; grid_points[ 883].i_atom = 0;
    grid_points[ 884].x =  -2.0842608672; grid_points[ 884].y =   4.3654758527; grid_points[ 884].z =  -2.0842608672; grid_points[ 884].w_fixed =    4.7844044707; grid_points[ 884].w_total = 3.2966811514; grid_points[ 884].i_atom = 0;
    grid_points[ 885].x =  -2.7423536150; grid_points[ 885].y =   5.7438484187; grid_points[ 885].z =  -2.7423536150; grid_points[ 885].w_fixed =   11.1501458946; grid_points[ 885].w_total = 0.7975816561; grid_points[ 885].i_atom = 0;
    grid_points[ 886].x =  -0.7447550142; grid_points[ 886].y =  -3.8828754963; grid_points[ 886].z =  -0.7447550142; grid_points[ 886].w_fixed =    1.7999446257; grid_points[ 886].w_total = 1.5182141465; grid_points[ 886].i_atom = 0;
    grid_points[ 887].x =  -0.9750809686; grid_points[ 887].y =  -5.0837093106; grid_points[ 887].z =  -0.9750809686; grid_points[ 887].w_fixed =    4.0944598961; grid_points[ 887].w_total = 1.3639957134; grid_points[ 887].i_atom = 0;
    grid_points[ 888].x =   0.0000000000; grid_points[ 888].y =  -4.3961380893; grid_points[ 888].z =   8.0701458556; grid_points[ 888].w_fixed =   27.1448404564; grid_points[ 888].w_total = 0.0000008316; grid_points[ 888].i_atom = 0;
    grid_points[ 889].x =   0.0000000000; grid_points[ 889].y =   3.5329984457; grid_points[ 889].z =  -1.9245685660; grid_points[ 889].w_fixed =    2.1250626994; grid_points[ 889].w_total = 1.9769819275; grid_points[ 889].i_atom = 0;
    grid_points[ 890].x =   0.0000000000; grid_points[ 890].y =   4.6256278652; grid_points[ 890].z =  -2.5197684414; grid_points[ 890].w_fixed =    4.8340287112; grid_points[ 890].w_total = 2.6774240028; grid_points[ 890].i_atom = 0;
    grid_points[ 891].x =  -0.7447550142; grid_points[ 891].y =   3.8828754963; grid_points[ 891].z =  -0.7447550142; grid_points[ 891].w_fixed =    1.7999446257; grid_points[ 891].w_total = 1.5182141465; grid_points[ 891].i_atom = 0;
    grid_points[ 892].x =  -0.9750809686; grid_points[ 892].y =   5.0837093106; grid_points[ 892].z =  -0.9750809686; grid_points[ 892].w_fixed =    4.0944598961; grid_points[ 892].w_total = 1.3639957134; grid_points[ 892].i_atom = 0;
    grid_points[ 893].x =   0.0000000000; grid_points[ 893].y =   4.3961380893; grid_points[ 893].z =   8.0701458556; grid_points[ 893].w_fixed =   27.1448404564; grid_points[ 893].w_total = 0.0000008316; grid_points[ 893].i_atom = 0;
    grid_points[ 894].x = -30.0795229432; grid_points[ 894].y = -10.0265076477; grid_points[ 894].z =  10.0265076477; grid_points[ 894].w_fixed = 3474.1698452873; grid_points[ 894].w_total = 0.0000000000; grid_points[ 894].i_atom = 0;
    grid_points[ 895].x = -19.1993004361; grid_points[ 895].y = -19.1993004361; grid_points[ 895].z =  19.1993004361; grid_points[ 895].w_fixed = 3632.6798826243; grid_points[ 895].w_total = 0.0000000000; grid_points[ 895].i_atom = 0;
    grid_points[ 896].x = -21.0705920995; grid_points[ 896].y =   7.0235306998; grid_points[ 896].z =  -7.0235306998; grid_points[ 896].w_fixed = 1092.3998745831; grid_points[ 896].w_total = 0.0000000000; grid_points[ 896].i_atom = 0;
    grid_points[ 897].x = -10.0265076477; grid_points[ 897].y = -30.0795229432; grid_points[ 897].z =  10.0265076477; grid_points[ 897].w_fixed = 3474.1698452873; grid_points[ 897].w_total = 0.0000000000; grid_points[ 897].i_atom = 0;
    grid_points[ 898].x = -11.8583704325; grid_points[ 898].y = -11.8583704325; grid_points[ 898].z =   0.0000000000; grid_points[ 898].w_fixed =  424.7625826047; grid_points[ 898].w_total = 0.0000000000; grid_points[ 898].i_atom = 0;
    grid_points[ 899].x =  -9.6823189135; grid_points[ 899].y =  -9.6823189135; grid_points[ 899].z =   9.6823189135; grid_points[ 899].w_fixed =  396.8927232114; grid_points[ 899].w_total = 0.0000000000; grid_points[ 899].i_atom = 0;
    grid_points[ 900].x = -10.0265076477; grid_points[ 900].y = -10.0265076477; grid_points[ 900].z =  30.0795229432; grid_points[ 900].w_fixed = 3474.1698452873; grid_points[ 900].w_total = 0.0000000000; grid_points[ 900].i_atom = 0;
    grid_points[ 901].x = -11.8913062075; grid_points[ 901].y =  -2.2808122310; grid_points[ 901].z =   2.2808122310; grid_points[ 901].w_fixed =   57.9035858166; grid_points[ 901].w_total = 0.0000000000; grid_points[ 901].i_atom = 0;
    grid_points[ 902].x = -11.8583704325; grid_points[ 902].y =   0.0000000000; grid_points[ 902].z =  11.8583704325; grid_points[ 902].w_fixed =  424.7625826047; grid_points[ 902].w_total = 0.0000000000; grid_points[ 902].i_atom = 0;
    grid_points[ 903].x =  -2.2808122310; grid_points[ 903].y = -11.8913062075; grid_points[ 903].z =   2.2808122310; grid_points[ 903].w_fixed =   57.9035858166; grid_points[ 903].w_total = 0.0000000000; grid_points[ 903].i_atom = 0;
    grid_points[ 904].x =  -1.2829568800; grid_points[ 904].y =  -1.2829568800; grid_points[ 904].z =  -6.6888597417; grid_points[ 904].w_fixed =    9.5422169009; grid_points[ 904].w_total = 0.0000055082; grid_points[ 904].i_atom = 0;
    grid_points[ 905].x =  -1.7011843293; grid_points[ 905].y =  -1.7011843293; grid_points[ 905].z =  -8.8693420267; grid_points[ 905].w_fixed =   22.9918908793; grid_points[ 905].w_total = 0.0000000000; grid_points[ 905].i_atom = 0;
    grid_points[ 906].x =  -1.2736799257; grid_points[ 906].y =  -1.2736799257; grid_points[ 906].z =   0.0000000000; grid_points[ 906].w_fixed =    0.4501131044; grid_points[ 906].w_total = 0.4501076870; grid_points[ 906].i_atom = 0;
    grid_points[ 907].x =  -1.6675833421; grid_points[ 907].y =  -1.6675833421; grid_points[ 907].z =   0.0000000000; grid_points[ 907].w_fixed =    0.9966614143; grid_points[ 907].w_total = 0.9964922905; grid_points[ 907].i_atom = 0;
    grid_points[ 908].x =  -1.3615760972; grid_points[ 908].y =  -1.3615760972; grid_points[ 908].z =   1.3615760972; grid_points[ 908].w_fixed =    0.9312676753; grid_points[ 908].w_total = 0.9312229918; grid_points[ 908].i_atom = 0;
    grid_points[ 909].x =  -1.7783851066; grid_points[ 909].y =  -1.7783851066; grid_points[ 909].z =   1.7783851066; grid_points[ 909].w_fixed =    0.9591496581; grid_points[ 909].w_total = 0.9586300412; grid_points[ 909].i_atom = 0;
    grid_points[ 910].x =  -2.1266717541; grid_points[ 910].y =  -2.1266717541; grid_points[ 910].z =   0.6652031738; grid_points[ 910].w_fixed =    0.9737495508; grid_points[ 910].w_total = 0.9722784395; grid_points[ 910].i_atom = 0;
    grid_points[ 911].x =  -1.4734978083; grid_points[ 911].y =  -2.7049519350; grid_points[ 911].z =   0.0000000000; grid_points[ 911].w_fixed =    0.9494794437; grid_points[ 911].w_total = 0.9424324035; grid_points[ 911].i_atom = 0;
    grid_points[ 912].x =  -2.7049519350; grid_points[ 912].y =  -1.4734978083; grid_points[ 912].z =   0.0000000000; grid_points[ 912].w_fixed =    0.9494794437; grid_points[ 912].w_total = 0.9424324035; grid_points[ 912].i_atom = 0;
    grid_points[ 913].x =  -2.3227887106; grid_points[ 913].y =  -2.3227887106; grid_points[ 913].z =   2.3227887106; grid_points[ 913].w_fixed =    2.1467059397; grid_points[ 913].w_total = 2.1405170372; grid_points[ 913].i_atom = 0;
    grid_points[ 914].x =  -2.7776937196; grid_points[ 914].y =  -2.7776937196; grid_points[ 914].z =   0.8688367985; grid_points[ 914].w_fixed =    2.1793824633; grid_points[ 914].w_total = 2.1482231938; grid_points[ 914].i_atom = 0;
    grid_points[ 915].x =  -3.0411437621; grid_points[ 915].y =  -3.0411437621; grid_points[ 915].z =   3.0411437621; grid_points[ 915].w_fixed =    4.8832621034; grid_points[ 915].w_total = 4.8025419330; grid_points[ 915].i_atom = 0;
    grid_points[ 916].x =  -2.2808122310; grid_points[ 916].y =  -2.2808122310; grid_points[ 916].z =  11.8913062075; grid_points[ 916].w_fixed =   57.9035858166; grid_points[ 916].w_total = 0.0000000000; grid_points[ 916].i_atom = 0;
    grid_points[ 917].x =  -1.2736799257; grid_points[ 917].y =   0.0000000000; grid_points[ 917].z =   1.2736799257; grid_points[ 917].w_fixed =    0.4501131044; grid_points[ 917].w_total = 0.4501069125; grid_points[ 917].i_atom = 0;
    grid_points[ 918].x =  -1.6292968354; grid_points[ 918].y =  -0.5430989451; grid_points[ 918].z =   0.5430989451; grid_points[ 918].w_fixed =    0.4022281471; grid_points[ 918].w_total = 0.4022123636; grid_points[ 918].i_atom = 0;
    grid_points[ 919].x =  -1.6675833421; grid_points[ 919].y =   0.0000000000; grid_points[ 919].z =   1.6675833421; grid_points[ 919].w_fixed =    0.9966614143; grid_points[ 919].w_total = 0.9962426478; grid_points[ 919].i_atom = 0;
    grid_points[ 920].x =  -2.1331797787; grid_points[ 920].y =  -0.7110599262; grid_points[ 920].z =   0.7110599262; grid_points[ 920].w_fixed =    0.8906323100; grid_points[ 920].w_total = 0.8900464448; grid_points[ 920].i_atom = 0;
    grid_points[ 921].x =  -2.9728265519; grid_points[ 921].y =  -0.5702030578; grid_points[ 921].z =   0.5702030578; grid_points[ 921].w_fixed =    0.8042164697; grid_points[ 921].w_total = 0.7903823761; grid_points[ 921].i_atom = 0;
    grid_points[ 922].x =  -2.1266717541; grid_points[ 922].y =  -0.6652031738; grid_points[ 922].z =   2.1266717541; grid_points[ 922].w_fixed =    0.9737495508; grid_points[ 922].w_total = 0.9690650946; grid_points[ 922].i_atom = 0;
    grid_points[ 923].x =  -2.5528215194; grid_points[ 923].y =  -1.2188238289; grid_points[ 923].z =   1.2188238289; grid_points[ 923].w_fixed =    0.9397324606; grid_points[ 923].w_total = 0.9361991824; grid_points[ 923].i_atom = 0;
    grid_points[ 924].x =  -1.4734978083; grid_points[ 924].y =   0.0000000000; grid_points[ 924].z =   2.7049519350; grid_points[ 924].w_fixed =    0.9494794437; grid_points[ 924].w_total = 0.8169727814; grid_points[ 924].i_atom = 0;
    grid_points[ 925].x =  -2.7049519350; grid_points[ 925].y =   0.0000000000; grid_points[ 925].z =   1.4734978083; grid_points[ 925].w_fixed =    0.9494794437; grid_points[ 925].w_total = 0.9435891428; grid_points[ 925].i_atom = 0;
    grid_points[ 926].x =  -2.7776937196; grid_points[ 926].y =  -0.8688367985; grid_points[ 926].z =   2.7776937196; grid_points[ 926].w_fixed =    2.1793824633; grid_points[ 926].w_total = 2.1192111689; grid_points[ 926].i_atom = 0;
    grid_points[ 927].x =  -1.9846144277; grid_points[ 927].y =   6.3448636630; grid_points[ 927].z =  -6.3448636630; grid_points[ 927].w_fixed =   27.8387029598; grid_points[ 927].w_total = 0.0002691714; grid_points[ 927].i_atom = 0;
    grid_points[ 928].x =  -1.2829568800; grid_points[ 928].y =   6.6888597417; grid_points[ 928].z =   1.2829568800; grid_points[ 928].w_fixed =    9.5422169009; grid_points[ 928].w_total = 0.0783957863; grid_points[ 928].i_atom = 0;
    grid_points[ 929].x =  -2.6608126953; grid_points[ 929].y =   8.5066870163; grid_points[ 929].z =  -8.5066870163; grid_points[ 929].w_fixed =   70.1099676542; grid_points[ 929].w_total = 0.0000000000; grid_points[ 929].i_atom = 0;
    grid_points[ 930].x =  -1.7011843293; grid_points[ 930].y =   8.8693420267; grid_points[ 930].z =   1.7011843293; grid_points[ 930].w_fixed =   22.9918908793; grid_points[ 930].w_total = 0.0000018741; grid_points[ 930].i_atom = 0;
    grid_points[ 931].x =  -0.5430989451; grid_points[ 931].y =  -1.6292968354; grid_points[ 931].z =   0.5430989451; grid_points[ 931].w_fixed =    0.4022281471; grid_points[ 931].w_total = 0.4022123636; grid_points[ 931].i_atom = 0;
    grid_points[ 932].x =  -0.7110599262; grid_points[ 932].y =  -2.1331797787; grid_points[ 932].z =   0.7110599262; grid_points[ 932].w_fixed =    0.8906323100; grid_points[ 932].w_total = 0.8900464448; grid_points[ 932].i_atom = 0;
    grid_points[ 933].x =  -0.5702030578; grid_points[ 933].y =  -2.9728265519; grid_points[ 933].z =   0.5702030578; grid_points[ 933].w_fixed =    0.8042164697; grid_points[ 933].w_total = 0.7903823761; grid_points[ 933].i_atom = 0;
    grid_points[ 934].x =  -0.6652031738; grid_points[ 934].y =  -2.1266717541; grid_points[ 934].z =   2.1266717541; grid_points[ 934].w_fixed =    0.9737495508; grid_points[ 934].w_total = 0.9690650946; grid_points[ 934].i_atom = 0;
    grid_points[ 935].x =  -1.2188238289; grid_points[ 935].y =  -2.5528215194; grid_points[ 935].z =   1.2188238289; grid_points[ 935].w_fixed =    0.9397324606; grid_points[ 935].w_total = 0.9361991824; grid_points[ 935].i_atom = 0;
    grid_points[ 936].x =   0.0000000000; grid_points[ 936].y =  -1.4734978083; grid_points[ 936].z =   2.7049519350; grid_points[ 936].w_fixed =    0.9494794437; grid_points[ 936].w_total = 0.8169727814; grid_points[ 936].i_atom = 0;
    grid_points[ 937].x =  -0.8688367985; grid_points[ 937].y =  -2.7776937196; grid_points[ 937].z =   2.7776937196; grid_points[ 937].w_fixed =    2.1793824633; grid_points[ 937].w_total = 2.1192111689; grid_points[ 937].i_atom = 0;
    grid_points[ 938].x =   0.0000000000; grid_points[ 938].y =   0.0000000000; grid_points[ 938].z =   0.0036626080; grid_points[ 938].w_fixed =    0.0000002129; grid_points[ 938].w_total = 0.0000002129; grid_points[ 938].i_atom = 0;
    grid_points[ 939].x =   0.0000000000; grid_points[ 939].y =   0.0000000000; grid_points[ 939].z =   0.0157155784; grid_points[ 939].w_fixed =    0.0000087099; grid_points[ 939].w_total = 0.0000087099; grid_points[ 939].i_atom = 0;
    grid_points[ 940].x =   0.0000000000; grid_points[ 940].y =   0.0000000000; grid_points[ 940].z =   0.0380278193; grid_points[ 940].w_fixed =    0.0000853157; grid_points[ 940].w_total = 0.0000853157; grid_points[ 940].i_atom = 0;
    grid_points[ 941].x =   0.0000000000; grid_points[ 941].y =   0.0000000000; grid_points[ 941].z =   0.0729054050; grid_points[ 941].w_fixed =    0.0004682262; grid_points[ 941].w_total = 0.0004682262; grid_points[ 941].i_atom = 0;
    grid_points[ 942].x =   0.0000000000; grid_points[ 942].y =   0.0000000000; grid_points[ 942].z =   0.1232101344; grid_points[ 942].w_fixed =    0.0018803520; grid_points[ 942].w_total = 0.0018803520; grid_points[ 942].i_atom = 0;
    grid_points[ 943].x =   0.0000000000; grid_points[ 943].y =   0.0000000000; grid_points[ 943].z =   0.1925158350; grid_points[ 943].w_fixed =    0.0062265394; grid_points[ 943].w_total = 0.0062265394; grid_points[ 943].i_atom = 0;
    grid_points[ 944].x =   0.0000000000; grid_points[ 944].y =   0.0000000000; grid_points[ 944].z =   0.2853164738; grid_points[ 944].w_fixed =    0.0181286111; grid_points[ 944].w_total = 0.0181286111; grid_points[ 944].i_atom = 0;
    grid_points[ 945].x =   0.0000000000; grid_points[ 945].y =   0.0000000000; grid_points[ 945].z =   0.4073062294; grid_points[ 945].w_fixed =    0.0482458934; grid_points[ 945].w_total = 0.0482458934; grid_points[ 945].i_atom = 0;
    grid_points[ 946].x =   0.0000000000; grid_points[ 946].y =   0.0000000000; grid_points[ 946].z =   0.5657608212; grid_points[ 946].w_fixed =    0.0344015986; grid_points[ 946].w_total = 0.0344015986; grid_points[ 946].i_atom = 0;
    grid_points[ 947].x =  -0.4000533132; grid_points[ 947].y =   0.0000000000; grid_points[ 947].z =   0.4000533132; grid_points[ 947].w_fixed =    0.0275212788; grid_points[ 947].w_total = 0.0275212788; grid_points[ 947].i_atom = 0;
    grid_points[ 948].x =  -0.4000533132; grid_points[ 948].y =  -0.4000533132; grid_points[ 948].z =   0.0000000000; grid_points[ 948].w_fixed =    0.0275212788; grid_points[ 948].w_total = 0.0275212788; grid_points[ 948].i_atom = 0;
    grid_points[ 949].x =  -0.3266421624; grid_points[ 949].y =  -0.3266421624; grid_points[ 949].z =   0.3266421624; grid_points[ 949].w_fixed =    0.0232210790; grid_points[ 949].w_total = 0.0232210790; grid_points[ 949].i_atom = 0;
    grid_points[ 950].x =   0.0000000000; grid_points[ 950].y =   0.0000000000; grid_points[ 950].z =   0.7700633400; grid_points[ 950].w_fixed =    0.0819768394; grid_points[ 950].w_total = 0.0819768393; grid_points[ 950].i_atom = 0;
    grid_points[ 951].x =  -0.5445170097; grid_points[ 951].y =   0.0000000000; grid_points[ 951].z =   0.5445170097; grid_points[ 951].w_fixed =    0.0655814715; grid_points[ 951].w_total = 0.0655814714; grid_points[ 951].i_atom = 0;
    grid_points[ 952].x =  -0.5445170097; grid_points[ 952].y =  -0.5445170097; grid_points[ 952].z =   0.0000000000; grid_points[ 952].w_fixed =    0.0655814715; grid_points[ 952].w_total = 0.0655814713; grid_points[ 952].i_atom = 0;
    grid_points[ 953].x =  -0.4445962766; grid_points[ 953].y =  -0.4445962766; grid_points[ 953].z =   0.4445962766; grid_points[ 953].w_fixed =    0.0553343666; grid_points[ 953].w_total = 0.0553343665; grid_points[ 953].i_atom = 0;
    grid_points[ 954].x =   0.0000000000; grid_points[ 954].y =   0.0000000000; grid_points[ 954].z =   1.0324394919; grid_points[ 954].w_fixed =    0.1890555184; grid_points[ 954].w_total = 0.1890554988; grid_points[ 954].i_atom = 0;
    grid_points[ 955].x =  -0.7300449659; grid_points[ 955].y =   0.0000000000; grid_points[ 955].z =   0.7300449659; grid_points[ 955].w_fixed =    0.1512444148; grid_points[ 955].w_total = 0.1512444119; grid_points[ 955].i_atom = 0;
    grid_points[ 956].x =  -0.7300449659; grid_points[ 956].y =  -0.7300449659; grid_points[ 956].z =   0.0000000000; grid_points[ 956].w_fixed =    0.1512444148; grid_points[ 956].w_total = 0.1512444075; grid_points[ 956].i_atom = 0;
    grid_points[ 957].x =  -0.5960792185; grid_points[ 957].y =  -0.5960792185; grid_points[ 957].z =   0.5960792185; grid_points[ 957].w_fixed =    0.1276124749; grid_points[ 957].w_total = 0.1276124734; grid_points[ 957].i_atom = 0;
    grid_points[ 958].x =   0.0000000000; grid_points[ 958].y =   0.0000000000; grid_points[ 958].z =   1.3690014933; grid_points[ 958].w_fixed =    0.4264811437; grid_points[ 958].w_total = 0.4264753435; grid_points[ 958].i_atom = 0;
    grid_points[ 959].x =  -0.9680302394; grid_points[ 959].y =   0.0000000000; grid_points[ 959].z =   0.9680302394; grid_points[ 959].w_fixed =    0.3411849149; grid_points[ 959].w_total = 0.3411847747; grid_points[ 959].i_atom = 0;
    grid_points[ 960].x =  -0.9680302394; grid_points[ 960].y =  -0.9680302394; grid_points[ 960].z =   0.0000000000; grid_points[ 960].w_fixed =    0.3411849149; grid_points[ 960].w_total = 0.3411846500; grid_points[ 960].i_atom = 0;
    grid_points[ 961].x =  -0.7903933807; grid_points[ 961].y =  -0.7903933807; grid_points[ 961].z =   0.7903933807; grid_points[ 961].w_fixed =    0.2878747720; grid_points[ 961].w_total = 0.2878747232; grid_points[ 961].i_atom = 0;
    grid_points[ 962].x =   0.0000000000; grid_points[ 962].y =   0.0000000000; grid_points[ 962].z =   1.8012554250; grid_points[ 962].w_fixed =    0.2531886212; grid_points[ 962].w_total = 0.2528919529; grid_points[ 962].i_atom = 0;
    grid_points[ 963].x =  -1.0399553045; grid_points[ 963].y =  -1.0399553045; grid_points[ 963].z =   1.0399553045; grid_points[ 963].w_fixed =    0.4205799265; grid_points[ 963].w_total = 0.4205787515; grid_points[ 963].i_atom = 0;
    grid_points[ 964].x =  -0.5430989451; grid_points[ 964].y =  -0.5430989451; grid_points[ 964].z =   1.6292968354; grid_points[ 964].w_fixed =    0.4022281471; grid_points[ 964].w_total = 0.4021170368; grid_points[ 964].i_atom = 0;
    grid_points[ 965].x =   0.0000000000; grid_points[ 965].y =   0.0000000000; grid_points[ 965].z =   2.3583189787; grid_points[ 965].w_fixed =    0.5606220456; grid_points[ 965].w_total = 0.5312379547; grid_points[ 965].i_atom = 0;
    grid_points[ 966].x =  -0.7110599262; grid_points[ 966].y =  -0.7110599262; grid_points[ 966].z =   2.1331797787; grid_points[ 966].w_fixed =    0.8906323100; grid_points[ 966].w_total = 0.8798148112; grid_points[ 966].i_atom = 0;
    grid_points[ 967].x =  -1.2188238289; grid_points[ 967].y =  -1.2188238289; grid_points[ 967].z =   2.5528215194; grid_points[ 967].w_fixed =    0.9397324606; grid_points[ 967].w_total = 0.8796756695; grid_points[ 967].i_atom = 0;
    grid_points[ 968].x =   0.0000000000; grid_points[ 968].y =   0.0000000000; grid_points[ 968].z =   3.0802533600; grid_points[ 968].w_fixed =    0.3749216611; grid_points[ 968].w_total = 0.1322139938; grid_points[ 968].i_atom = 0;
    grid_points[ 969].x =  -0.5702030578; grid_points[ 969].y =  -0.5702030578; grid_points[ 969].z =   2.9728265519; grid_points[ 969].w_fixed =    0.8042164697; grid_points[ 969].w_total = 0.4464214906; grid_points[ 969].i_atom = 0;
    grid_points[ 970].x =   0.0000000000; grid_points[ 970].y =   6.9305700600; grid_points[ 970].z =   0.0000000000; grid_points[ 970].w_fixed =    4.4485333807; grid_points[ 970].w_total = 0.0169320713; grid_points[ 970].i_atom = 0;
    grid_points[ 971].x =   0.0000000000; grid_points[ 971].y =   9.1898468013; grid_points[ 971].z =   0.0000000000; grid_points[ 971].w_fixed =   10.7187035387; grid_points[ 971].w_total = 0.0000000000; grid_points[ 971].i_atom = 0;
    grid_points[ 972].x =   0.0000000000; grid_points[ 972].y =  16.7702682933; grid_points[ 972].z =   0.0000000000; grid_points[ 972].w_fixed =  238.9289527151; grid_points[ 972].w_total = 0.0000000000; grid_points[ 972].i_atom = 0;
    grid_points[ 973].x =   6.3448636630; grid_points[ 973].y =  -1.9846144277; grid_points[ 973].z =  -6.3448636630; grid_points[ 973].w_fixed =   27.8387029598; grid_points[ 973].w_total = 0.0002691714; grid_points[ 973].i_atom = 0;
    grid_points[ 974].x =   6.6888597417; grid_points[ 974].y =  -1.2829568800; grid_points[ 974].z =   1.2829568800; grid_points[ 974].w_fixed =    9.5422169009; grid_points[ 974].w_total = 0.0783957863; grid_points[ 974].i_atom = 0;
    grid_points[ 975].x =   6.3448636630; grid_points[ 975].y =   6.3448636630; grid_points[ 975].z =   1.9846144277; grid_points[ 975].w_fixed =   27.8387029598; grid_points[ 975].w_total = 0.0012187121; grid_points[ 975].i_atom = 0;
    grid_points[ 976].x =   7.0235306998; grid_points[ 976].y = -21.0705920995; grid_points[ 976].z =  -7.0235306998; grid_points[ 976].w_fixed = 1092.3998745831; grid_points[ 976].w_total = 0.0000000000; grid_points[ 976].i_atom = 0;
    grid_points[ 977].x =   8.5066870163; grid_points[ 977].y =  -2.6608126953; grid_points[ 977].z =  -8.5066870163; grid_points[ 977].w_fixed =   70.1099676542; grid_points[ 977].w_total = 0.0000000000; grid_points[ 977].i_atom = 0;
    grid_points[ 978].x =   8.8693420267; grid_points[ 978].y =  -1.7011843293; grid_points[ 978].z =   1.7011843293; grid_points[ 978].w_fixed =   22.9918908793; grid_points[ 978].w_total = 0.0000018741; grid_points[ 978].i_atom = 0;
    grid_points[ 979].x =   7.1135404262; grid_points[ 979].y =   7.1135404262; grid_points[ 979].z =  -7.1135404262; grid_points[ 979].w_fixed =   69.0587753800; grid_points[ 979].w_total = 0.0000000000; grid_points[ 979].i_atom = 0;
    grid_points[ 980].x =   8.5066870163; grid_points[ 980].y =   8.5066870163; grid_points[ 980].z =   2.6608126953; grid_points[ 980].w_fixed =   70.1099676542; grid_points[ 980].w_total = 0.0000000000; grid_points[ 980].i_atom = 0;
    grid_points[ 981].x =   7.0235306998; grid_points[ 981].y =   7.0235306998; grid_points[ 981].z =  21.0705920995; grid_points[ 981].w_fixed = 1092.3998745831; grid_points[ 981].w_total = 0.0000000000; grid_points[ 981].i_atom = 0;
    grid_points[ 982].x =  16.4716395421; grid_points[ 982].y =  16.4716395421; grid_points[ 982].z =   0.0000000000; grid_points[ 982].w_fixed = 1222.4492551854; grid_points[ 982].w_total = 0.0000000000; grid_points[ 982].i_atom = 0;
    grid_points[ 983].x = -10.8198077400; grid_points[ 983].y =   0.0000000000; grid_points[ 983].z =   5.8939912332; grid_points[ 983].w_fixed =   68.3625199466; grid_points[ 983].w_total = 0.0000000000; grid_points[ 983].i_atom = 0;
    grid_points[ 984].x =  -2.7423536150; grid_points[ 984].y =  -2.7423536150; grid_points[ 984].z =  -5.7438484187; grid_points[ 984].w_fixed =   11.1501458946; grid_points[ 984].w_total = 0.2040185179; grid_points[ 984].i_atom = 0;
    grid_points[ 985].x =  -1.5919331643; grid_points[ 985].y =  -1.5919331643; grid_points[ 985].z =  -3.3342974947; grid_points[ 985].w_fixed =    2.1032476402; grid_points[ 985].w_total = 1.9796616798; grid_points[ 985].i_atom = 0;
    grid_points[ 986].x =  -2.0842608672; grid_points[ 986].y =  -2.0842608672; grid_points[ 986].z =  -4.3654758527; grid_points[ 986].w_fixed =    4.7844044707; grid_points[ 986].w_total = 2.7579173487; grid_points[ 986].i_atom = 0;
    grid_points[ 987].x =  -1.5919331643; grid_points[ 987].y =  -1.5919331643; grid_points[ 987].z =   3.3342974947; grid_points[ 987].w_fixed =    2.1032476402; grid_points[ 987].w_total = 1.2794689390; grid_points[ 987].i_atom = 0;
    grid_points[ 988].x =  -2.0842608672; grid_points[ 988].y =  -2.0842608672; grid_points[ 988].z =   4.3654758527; grid_points[ 988].w_fixed =    4.7844044707; grid_points[ 988].w_total = 1.0736736243; grid_points[ 988].i_atom = 0;
    grid_points[ 989].x =  -2.7423536150; grid_points[ 989].y =  -2.7423536150; grid_points[ 989].z =   5.7438484187; grid_points[ 989].w_fixed =   11.1501458946; grid_points[ 989].w_total = 0.0554063878; grid_points[ 989].i_atom = 0;
    grid_points[ 990].x =  -1.9245685660; grid_points[ 990].y =   0.0000000000; grid_points[ 990].z =   3.5329984457; grid_points[ 990].w_fixed =    2.1250626994; grid_points[ 990].w_total = 0.7185910872; grid_points[ 990].i_atom = 0;
    grid_points[ 991].x =  -2.5197684414; grid_points[ 991].y =   0.0000000000; grid_points[ 991].z =   4.6256278652; grid_points[ 991].w_fixed =    4.8340287112; grid_points[ 991].w_total = 0.2039773991; grid_points[ 991].i_atom = 0;
    grid_points[ 992].x =   0.0000000000; grid_points[ 992].y =  -1.9245685660; grid_points[ 992].z =   3.5329984457; grid_points[ 992].w_fixed =    2.1250626994; grid_points[ 992].w_total = 0.7185910872; grid_points[ 992].i_atom = 0;
    grid_points[ 993].x =   0.0000000000; grid_points[ 993].y =  -2.5197684414; grid_points[ 993].z =   4.6256278652; grid_points[ 993].w_fixed =    4.8340287112; grid_points[ 993].w_total = 0.2039773991; grid_points[ 993].i_atom = 0;
    grid_points[ 994].x =  -0.7447550142; grid_points[ 994].y =  -0.7447550142; grid_points[ 994].z =  -3.8828754963; grid_points[ 994].w_fixed =    1.7999446257; grid_points[ 994].w_total = 0.6776187779; grid_points[ 994].i_atom = 0;
    grid_points[ 995].x =  -0.9750809686; grid_points[ 995].y =  -0.9750809686; grid_points[ 995].z =  -5.0837093106; grid_points[ 995].w_fixed =    4.0944598961; grid_points[ 995].w_total = 0.0118697638; grid_points[ 995].i_atom = 0;
    grid_points[ 996].x =   0.0000000000; grid_points[ 996].y =   0.0000000000; grid_points[ 996].z =   4.0231880620; grid_points[ 996].w_fixed =    0.8391251042; grid_points[ 996].w_total = 0.0000017193; grid_points[ 996].i_atom = 0;
    grid_points[ 997].x =  -0.7447550142; grid_points[ 997].y =  -0.7447550142; grid_points[ 997].z =   3.8828754963; grid_points[ 997].w_fixed =    1.7999446257; grid_points[ 997].w_total = 0.0335908066; grid_points[ 997].i_atom = 0;
    grid_points[ 998].x =  -0.9750809686; grid_points[ 998].y =  -0.9750809686; grid_points[ 998].z =   5.0837093106; grid_points[ 998].w_fixed =    4.0944598961; grid_points[ 998].w_total = 0.0001490129; grid_points[ 998].i_atom = 0;
    grid_points[ 999].x =   0.0000000000; grid_points[ 999].y =   8.0701458556; grid_points[ 999].z =  -4.3961380893; grid_points[ 999].w_fixed =   27.1448404564; grid_points[ 999].w_total = 0.0000094033; grid_points[ 999].i_atom = 0;
    grid_points[1000].x =   0.0000000000; grid_points[1000].y =   8.0701458556; grid_points[1000].z =   4.3961380893; grid_points[1000].w_fixed =   27.1448404564; grid_points[1000].w_total = 0.0000023325; grid_points[1000].i_atom = 0;
    grid_points[1001].x =   8.0701458556; grid_points[1001].y =   0.0000000000; grid_points[1001].z =   4.3961380893; grid_points[1001].w_fixed =   27.1448404564; grid_points[1001].w_total = 0.0000023325; grid_points[1001].i_atom = 0;
    grid_points[1002].x = -30.0795229432; grid_points[1002].y = -10.0265076477; grid_points[1002].z = -10.0265076477; grid_points[1002].w_fixed = 3474.1698452873; grid_points[1002].w_total = 0.0000000000; grid_points[1002].i_atom = 0;
    grid_points[1003].x = -19.1993004361; grid_points[1003].y = -19.1993004361; grid_points[1003].z = -19.1993004361; grid_points[1003].w_fixed = 3632.6798826243; grid_points[1003].w_total = 0.0000000000; grid_points[1003].i_atom = 0;
    grid_points[1004].x = -21.0705920995; grid_points[1004].y =   7.0235306998; grid_points[1004].z =   7.0235306998; grid_points[1004].w_fixed = 1092.3998745831; grid_points[1004].w_total = 0.0000000000; grid_points[1004].i_atom = 0;
    grid_points[1005].x = -10.0265076477; grid_points[1005].y = -30.0795229432; grid_points[1005].z = -10.0265076477; grid_points[1005].w_fixed = 3474.1698452873; grid_points[1005].w_total = 0.0000000000; grid_points[1005].i_atom = 0;
    grid_points[1006].x = -10.0265076477; grid_points[1006].y = -10.0265076477; grid_points[1006].z = -30.0795229432; grid_points[1006].w_fixed = 3474.1698452873; grid_points[1006].w_total = 0.0000000000; grid_points[1006].i_atom = 0;
    grid_points[1007].x =  -9.6823189135; grid_points[1007].y =  -9.6823189135; grid_points[1007].z =  -9.6823189135; grid_points[1007].w_fixed =  396.8927232114; grid_points[1007].w_total = 0.0000000000; grid_points[1007].i_atom = 0;
    grid_points[1008].x = -11.8913062075; grid_points[1008].y =  -2.2808122310; grid_points[1008].z =  -2.2808122310; grid_points[1008].w_fixed =   57.9035858166; grid_points[1008].w_total = 0.0000000000; grid_points[1008].i_atom = 0;
    grid_points[1009].x =  -2.2808122310; grid_points[1009].y = -11.8913062075; grid_points[1009].z =  -2.2808122310; grid_points[1009].w_fixed =   57.9035858166; grid_points[1009].w_total = 0.0000000000; grid_points[1009].i_atom = 0;
    grid_points[1010].x =  -2.2808122310; grid_points[1010].y =  -2.2808122310; grid_points[1010].z = -11.8913062075; grid_points[1010].w_fixed =   57.9035858166; grid_points[1010].w_total = 0.0000000000; grid_points[1010].i_atom = 0;
    grid_points[1011].x =  -1.3615760972; grid_points[1011].y =  -1.3615760972; grid_points[1011].z =  -1.3615760972; grid_points[1011].w_fixed =    0.9312676753; grid_points[1011].w_total = 0.9312379628; grid_points[1011].i_atom = 0;
    grid_points[1012].x =  -1.7783851066; grid_points[1012].y =  -1.7783851066; grid_points[1012].z =  -1.7783851066; grid_points[1012].w_fixed =    0.9591496581; grid_points[1012].w_total = 0.9588589425; grid_points[1012].i_atom = 0;
    grid_points[1013].x =  -2.3227887106; grid_points[1013].y =  -2.3227887106; grid_points[1013].z =  -2.3227887106; grid_points[1013].w_fixed =    2.1467059397; grid_points[1013].w_total = 2.1417068903; grid_points[1013].i_atom = 0;
    grid_points[1014].x =  -3.0411437621; grid_points[1014].y =  -3.0411437621; grid_points[1014].z =  -3.0411437621; grid_points[1014].w_fixed =    4.8832621034; grid_points[1014].w_total = 4.8017896257; grid_points[1014].i_atom = 0;
    grid_points[1015].x =  -2.1266717541; grid_points[1015].y =  -2.1266717541; grid_points[1015].z =  -0.6652031738; grid_points[1015].w_fixed =    0.9737495508; grid_points[1015].w_total = 0.9721447960; grid_points[1015].i_atom = 0;
    grid_points[1016].x =  -2.7776937196; grid_points[1016].y =  -2.7776937196; grid_points[1016].z =  -0.8688367985; grid_points[1016].w_fixed =    2.1793824633; grid_points[1016].w_total = 2.1466862309; grid_points[1016].i_atom = 0;
    grid_points[1017].x =  -1.2829568800; grid_points[1017].y =  -1.2829568800; grid_points[1017].z =   6.6888597417; grid_points[1017].w_fixed =    9.5422169009; grid_points[1017].w_total = 0.0000003249; grid_points[1017].i_atom = 0;
    grid_points[1018].x =  -1.7011843293; grid_points[1018].y =  -1.7011843293; grid_points[1018].z =   8.8693420267; grid_points[1018].w_fixed =   22.9918908793; grid_points[1018].w_total = 0.0000000000; grid_points[1018].i_atom = 0;
    grid_points[1019].x =  -2.1266717541; grid_points[1019].y =  -0.6652031738; grid_points[1019].z =  -2.1266717541; grid_points[1019].w_fixed =    0.9737495508; grid_points[1019].w_total = 0.9730669536; grid_points[1019].i_atom = 0;
    grid_points[1020].x =  -2.7776937196; grid_points[1020].y =  -0.8688367985; grid_points[1020].z =  -2.7776937196; grid_points[1020].w_fixed =    2.1793824633; grid_points[1020].w_total = 2.1646215031; grid_points[1020].i_atom = 0;
    grid_points[1021].x =  -1.6292968354; grid_points[1021].y =  -0.5430989451; grid_points[1021].z =  -0.5430989451; grid_points[1021].w_fixed =    0.4022281471; grid_points[1021].w_total = 0.4022099166; grid_points[1021].i_atom = 0;
    grid_points[1022].x =  -2.1331797787; grid_points[1022].y =  -0.7110599262; grid_points[1022].z =  -0.7110599262; grid_points[1022].w_fixed =    0.8906323100; grid_points[1022].w_total = 0.8899740741; grid_points[1022].i_atom = 0;
    grid_points[1023].x =  -2.9728265519; grid_points[1023].y =  -0.5702030578; grid_points[1023].z =  -0.5702030578; grid_points[1023].w_fixed =    0.8042164697; grid_points[1023].w_total = 0.7899201355; grid_points[1023].i_atom = 0;
    grid_points[1024].x =  -2.5528215194; grid_points[1024].y =  -1.2188238289; grid_points[1024].z =  -1.2188238289; grid_points[1024].w_fixed =    0.9397324606; grid_points[1024].w_total = 0.9357421258; grid_points[1024].i_atom = 0;
    grid_points[1025].x =  -1.2829568800; grid_points[1025].y =   6.6888597417; grid_points[1025].z =  -1.2829568800; grid_points[1025].w_fixed =    9.5422169009; grid_points[1025].w_total = 0.0865328484; grid_points[1025].i_atom = 0;
    grid_points[1026].x =  -1.9846144277; grid_points[1026].y =   6.3448636630; grid_points[1026].z =   6.3448636630; grid_points[1026].w_fixed =   27.8387029598; grid_points[1026].w_total = 0.0001758536; grid_points[1026].i_atom = 0;
    grid_points[1027].x =  -1.7011843293; grid_points[1027].y =   8.8693420267; grid_points[1027].z =  -1.7011843293; grid_points[1027].w_fixed =   22.9918908793; grid_points[1027].w_total = 0.0000042871; grid_points[1027].i_atom = 0;
    grid_points[1028].x =  -2.6608126953; grid_points[1028].y =   8.5066870163; grid_points[1028].z =   8.5066870163; grid_points[1028].w_fixed =   70.1099676542; grid_points[1028].w_total = 0.0000000000; grid_points[1028].i_atom = 0;
    grid_points[1029].x =  -0.6652031738; grid_points[1029].y =  -2.1266717541; grid_points[1029].z =  -2.1266717541; grid_points[1029].w_fixed =    0.9737495508; grid_points[1029].w_total = 0.9730669536; grid_points[1029].i_atom = 0;
    grid_points[1030].x =  -0.8688367985; grid_points[1030].y =  -2.7776937196; grid_points[1030].z =  -2.7776937196; grid_points[1030].w_fixed =    2.1793824633; grid_points[1030].w_total = 2.1646215031; grid_points[1030].i_atom = 0;
    grid_points[1031].x =  -0.5430989451; grid_points[1031].y =  -1.6292968354; grid_points[1031].z =  -0.5430989451; grid_points[1031].w_fixed =    0.4022281471; grid_points[1031].w_total = 0.4022099166; grid_points[1031].i_atom = 0;
    grid_points[1032].x =  -0.7110599262; grid_points[1032].y =  -2.1331797787; grid_points[1032].z =  -0.7110599262; grid_points[1032].w_fixed =    0.8906323100; grid_points[1032].w_total = 0.8899740741; grid_points[1032].i_atom = 0;
    grid_points[1033].x =  -0.5702030578; grid_points[1033].y =  -2.9728265519; grid_points[1033].z =  -0.5702030578; grid_points[1033].w_fixed =    0.8042164697; grid_points[1033].w_total = 0.7899201355; grid_points[1033].i_atom = 0;
    grid_points[1034].x =  -1.2188238289; grid_points[1034].y =  -2.5528215194; grid_points[1034].z =  -1.2188238289; grid_points[1034].w_fixed =    0.9397324606; grid_points[1034].w_total = 0.9357421258; grid_points[1034].i_atom = 0;
    grid_points[1035].x =  -0.5430989451; grid_points[1035].y =  -0.5430989451; grid_points[1035].z =  -1.6292968354; grid_points[1035].w_fixed =    0.4022281471; grid_points[1035].w_total = 0.4022275009; grid_points[1035].i_atom = 0;
    grid_points[1036].x =  -0.7110599262; grid_points[1036].y =  -0.7110599262; grid_points[1036].z =  -2.1331797787; grid_points[1036].w_fixed =    0.8906323100; grid_points[1036].w_total = 0.8904412859; grid_points[1036].i_atom = 0;
    grid_points[1037].x =  -0.5702030578; grid_points[1037].y =  -0.5702030578; grid_points[1037].z =  -2.9728265519; grid_points[1037].w_fixed =    0.8042164697; grid_points[1037].w_total = 0.7676826133; grid_points[1037].i_atom = 0;
    grid_points[1038].x =  -1.2188238289; grid_points[1038].y =  -1.2188238289; grid_points[1038].z =  -2.5528215194; grid_points[1038].w_fixed =    0.9397324606; grid_points[1038].w_total = 0.9373660293; grid_points[1038].i_atom = 0;
    grid_points[1039].x =  -0.3266421624; grid_points[1039].y =  -0.3266421624; grid_points[1039].z =  -0.3266421624; grid_points[1039].w_fixed =    0.0232210790; grid_points[1039].w_total = 0.0232210790; grid_points[1039].i_atom = 0;
    grid_points[1040].x =  -0.4445962766; grid_points[1040].y =  -0.4445962766; grid_points[1040].z =  -0.4445962766; grid_points[1040].w_fixed =    0.0553343666; grid_points[1040].w_total = 0.0553343665; grid_points[1040].i_atom = 0;
    grid_points[1041].x =  -0.5960792185; grid_points[1041].y =  -0.5960792185; grid_points[1041].z =  -0.5960792185; grid_points[1041].w_fixed =    0.1276124749; grid_points[1041].w_total = 0.1276124729; grid_points[1041].i_atom = 0;
    grid_points[1042].x =  -0.7903933807; grid_points[1042].y =  -0.7903933807; grid_points[1042].z =  -0.7903933807; grid_points[1042].w_fixed =    0.2878747720; grid_points[1042].w_total = 0.2878747082; grid_points[1042].i_atom = 0;
    grid_points[1043].x =  -1.0399553045; grid_points[1043].y =  -1.0399553045; grid_points[1043].z =  -1.0399553045; grid_points[1043].w_fixed =    0.4205799265; grid_points[1043].w_total = 0.4205787392; grid_points[1043].i_atom = 0;
    grid_points[1044].x =   0.0000000000; grid_points[1044].y =  16.4716395421; grid_points[1044].z =  16.4716395421; grid_points[1044].w_fixed = 1222.4492551854; grid_points[1044].w_total = 0.0000000000; grid_points[1044].i_atom = 0;
    grid_points[1045].x =   6.6888597417; grid_points[1045].y =  -1.2829568800; grid_points[1045].z =  -1.2829568800; grid_points[1045].w_fixed =    9.5422169009; grid_points[1045].w_total = 0.0865328484; grid_points[1045].i_atom = 0;
    grid_points[1046].x =   6.3448636630; grid_points[1046].y =  -1.9846144277; grid_points[1046].z =   6.3448636630; grid_points[1046].w_fixed =   27.8387029598; grid_points[1046].w_total = 0.0001758536; grid_points[1046].i_atom = 0;
    grid_points[1047].x =   6.3448636630; grid_points[1047].y =   6.3448636630; grid_points[1047].z =  -1.9846144277; grid_points[1047].w_fixed =   27.8387029598; grid_points[1047].w_total = 0.0013808480; grid_points[1047].i_atom = 0;
    grid_points[1048].x =   7.0235306998; grid_points[1048].y = -21.0705920995; grid_points[1048].z =   7.0235306998; grid_points[1048].w_fixed = 1092.3998745831; grid_points[1048].w_total = 0.0000000000; grid_points[1048].i_atom = 0;
    grid_points[1049].x =   8.8693420267; grid_points[1049].y =  -1.7011843293; grid_points[1049].z =  -1.7011843293; grid_points[1049].w_fixed =   22.9918908793; grid_points[1049].w_total = 0.0000042871; grid_points[1049].i_atom = 0;
    grid_points[1050].x =   8.5066870163; grid_points[1050].y =  -2.6608126953; grid_points[1050].z =   8.5066870163; grid_points[1050].w_fixed =   70.1099676542; grid_points[1050].w_total = 0.0000000000; grid_points[1050].i_atom = 0;
    grid_points[1051].x =   7.0235306998; grid_points[1051].y =   7.0235306998; grid_points[1051].z = -21.0705920995; grid_points[1051].w_fixed = 1092.3998745831; grid_points[1051].w_total = 0.0000000000; grid_points[1051].i_atom = 0;
    grid_points[1052].x =   8.5066870163; grid_points[1052].y =   8.5066870163; grid_points[1052].z =  -2.6608126953; grid_points[1052].w_fixed =   70.1099676542; grid_points[1052].w_total = 0.0000000000; grid_points[1052].i_atom = 0;
    grid_points[1053].x =   7.1135404262; grid_points[1053].y =   7.1135404262; grid_points[1053].z =   7.1135404262; grid_points[1053].w_fixed =   69.0587753800; grid_points[1053].w_total = 0.0000000000; grid_points[1053].i_atom = 0;
    grid_points[1054].x =  16.4716395421; grid_points[1054].y =   0.0000000000; grid_points[1054].z =  16.4716395421; grid_points[1054].w_fixed = 1222.4492551854; grid_points[1054].w_total = 0.0000000000; grid_points[1054].i_atom = 0;
    grid_points[1055].x =  -7.5660149734; grid_points[1055].y =   2.5220049911; grid_points[1055].z =   1.6731867016; grid_points[1055].w_fixed =   55.2891360098; grid_points[1055].w_total = 0.0000000000; grid_points[1055].i_atom = 1;
    grid_points[1056].x =   2.5220049911; grid_points[1056].y =  -7.5660149734; grid_points[1056].z =   1.6731867016; grid_points[1056].w_fixed =   55.2891360098; grid_points[1056].w_total = 0.0000000000; grid_points[1056].i_atom = 1;
    grid_points[1057].x =   0.3227069453; grid_points[1057].y =   0.3227069453; grid_points[1057].z =   2.5127177700; grid_points[1057].w_fixed =    0.1518581277; grid_points[1057].w_total = 0.0160452730; grid_points[1057].i_atom = 1;
    grid_points[1058].x =   0.4279052608; grid_points[1058].y =   0.4279052608; grid_points[1058].z =   1.9642529008; grid_points[1058].w_fixed =    0.3659008737; grid_points[1058].w_total = 0.0015586987; grid_points[1058].i_atom = 1;
    grid_points[1059].x =   1.5959473017; grid_points[1059].y =   0.4991974941; grid_points[1059].z =   2.5992443910; grid_points[1059].w_fixed =    0.4430347113; grid_points[1059].w_total = 0.0354740456; grid_points[1059].i_atom = 1;
    grid_points[1060].x =   0.4991974941; grid_points[1060].y =   1.5959473017; grid_points[1060].z =   2.5992443910; grid_points[1060].w_fixed =    0.4430347113; grid_points[1060].w_total = 0.0354740456; grid_points[1060].i_atom = 1;
    grid_points[1061].x =   0.9146585961; grid_points[1061].y =   0.9146585961; grid_points[1061].z =   2.2794430080; grid_points[1061].w_fixed =    0.4275576805; grid_points[1061].w_total = 0.0104656786; grid_points[1061].i_atom = 1;
    grid_points[1062].x =   1.7892954446; grid_points[1062].y =   1.7892954446; grid_points[1062].z =   2.4058962481; grid_points[1062].w_fixed =    1.0990251469; grid_points[1062].w_total = 0.0123810957; grid_points[1062].i_atom = 1;
    grid_points[1063].x =   0.5737012360; grid_points[1063].y =   0.5737012360; grid_points[1063].z =   1.2041269412; grid_points[1063].w_fixed =    0.9214976165; grid_points[1063].w_total = 0.0000026829; grid_points[1063].i_atom = 1;
    grid_points[1064].x =   2.1397188201; grid_points[1064].y =   0.6692841749; grid_points[1064].z =   2.0554728727; grid_points[1064].w_fixed =    1.1157541830; grid_points[1064].w_total = 0.0031940004; grid_points[1064].i_atom = 1;
    grid_points[1065].x =   0.6692841749; grid_points[1065].y =   2.1397188201; grid_points[1065].z =   2.0554728727; grid_points[1065].w_fixed =    1.1157541830; grid_points[1065].w_total = 0.0031940004; grid_points[1065].i_atom = 1;
    grid_points[1066].x =   1.2263012757; grid_points[1066].y =   1.2263012757; grid_points[1066].z =   1.6267086916; grid_points[1066].w_fixed =    1.0767762850; grid_points[1066].w_total = 0.0002888020; grid_points[1066].i_atom = 1;
    grid_points[1067].x =   2.4354299108; grid_points[1067].y =   2.4354299108; grid_points[1067].z =   1.7597617820; grid_points[1067].w_fixed =    6.3162875542; grid_points[1067].w_total = 0.0001328651; grid_points[1067].i_atom = 1;
    grid_points[1068].x =   1.2718617903; grid_points[1068].y =   1.2718617903; grid_points[1068].z =   0.3796063218; grid_points[1068].w_fixed =    6.0406797362; grid_points[1068].w_total = 0.0000000021; grid_points[1068].i_atom = 1;
    grid_points[1069].x =   0.2452657651; grid_points[1069].y =   0.2452657651; grid_points[1069].z =   2.9164672649; grid_points[1069].w_fixed =    0.0651606456; grid_points[1069].w_total = 0.0288023961; grid_points[1069].i_atom = 1;
    grid_points[1070].x =   0.5242619359; grid_points[1070].y =   0.5242619359; grid_points[1070].z =   3.0971272145; grid_points[1070].w_fixed =    0.0761406613; grid_points[1070].w_total = 0.0448544141; grid_points[1070].i_atom = 1;
    grid_points[1071].x =   1.2035918363; grid_points[1071].y =   0.3764723484; grid_points[1071].z =   2.9915998565; grid_points[1071].w_fixed =    0.1838706234; grid_points[1071].w_total = 0.0694127716; grid_points[1071].i_atom = 1;
    grid_points[1072].x =   0.3764723484; grid_points[1072].y =   1.2035918363; grid_points[1072].z =   2.9915998565; grid_points[1072].w_fixed =    0.1838706234; grid_points[1072].w_total = 0.0694127716; grid_points[1072].i_atom = 1;
    grid_points[1073].x =   0.6897944676; grid_points[1073].y =   0.6897944676; grid_points[1073].z =   2.7504200046; grid_points[1073].w_fixed =    0.1774472638; grid_points[1073].w_total = 0.0390890156; grid_points[1073].i_atom = 1;
    grid_points[1074].x =   1.3345778006; grid_points[1074].y =   1.3345778006; grid_points[1074].z =   2.8606138921; grid_points[1074].w_fixed =    0.4363920800; grid_points[1074].w_total = 0.0730945598; grid_points[1074].i_atom = 1;
    grid_points[1075].x =   1.2263012757; grid_points[1075].y =   2.5684830011; grid_points[1075].z =   2.9688904170; grid_points[1075].w_fixed =    1.0767762850; grid_points[1075].w_total = 0.0592143299; grid_points[1075].i_atom = 1;
    grid_points[1076].x =   2.5684830011; grid_points[1076].y =   1.2263012757; grid_points[1076].z =   2.9688904170; grid_points[1076].w_fixed =    1.0767762850; grid_points[1076].w_total = 0.0592143299; grid_points[1076].i_atom = 1;
    grid_points[1077].x =   1.7666549613; grid_points[1077].y =   1.7666549613; grid_points[1077].z =   9.4951565767; grid_points[1077].w_fixed =   17.3848280115; grid_points[1077].w_total = 0.0000000000; grid_points[1077].i_atom = 1;
    grid_points[1078].x =   2.5220049911; grid_points[1078].y =   2.5220049911; grid_points[1078].z =  11.7612066662; grid_points[1078].w_fixed =   55.2891360098; grid_points[1078].w_total = 0.0000000000; grid_points[1078].i_atom = 1;
    grid_points[1079].x =   0.0000000000; grid_points[1079].y =  -8.3645442750; grid_points[1079].z =   4.1951916927; grid_points[1079].w_fixed =   34.8025870819; grid_points[1079].w_total = 0.0000000000; grid_points[1079].i_atom = 1;
    grid_points[1080].x =   2.5220049911; grid_points[1080].y =   2.5220049911; grid_points[1080].z =  -3.3708232807; grid_points[1080].w_fixed =   55.2891360098; grid_points[1080].w_total = 0.0000000000; grid_points[1080].i_atom = 1;
    grid_points[1081].x =   0.1006269070; grid_points[1081].y =   0.1006269070; grid_points[1081].z =   4.1951916927; grid_points[1081].w_fixed =    0.0004379831; grid_points[1081].w_total = 0.0004379831; grid_points[1081].i_atom = 1;
    grid_points[1082].x =   0.0821615255; grid_points[1082].y =   0.0821615255; grid_points[1082].z =   4.2773532183; grid_points[1082].w_fixed =    0.0003695483; grid_points[1082].w_total = 0.0003695483; grid_points[1082].i_atom = 1;
    grid_points[1083].x =   0.0821615255; grid_points[1083].y =   0.0821615255; grid_points[1083].z =   4.1130301672; grid_points[1083].w_fixed =    0.0003695483; grid_points[1083].w_total = 0.0003695482; grid_points[1083].i_atom = 1;
    grid_points[1084].x =   0.1369644012; grid_points[1084].y =   0.1369644012; grid_points[1084].z =   4.1951916927; grid_points[1084].w_fixed =    0.0010436861; grid_points[1084].w_total = 0.0010436861; grid_points[1084].i_atom = 1;
    grid_points[1085].x =   0.1118309653; grid_points[1085].y =   0.1118309653; grid_points[1085].z =   4.3070226580; grid_points[1085].w_fixed =    0.0008806102; grid_points[1085].w_total = 0.0008806101; grid_points[1085].i_atom = 1;
    grid_points[1086].x =   0.1118309653; grid_points[1086].y =   0.1118309653; grid_points[1086].z =   4.0833607275; grid_points[1086].w_fixed =    0.0008806102; grid_points[1086].w_total = 0.0008806095; grid_points[1086].i_atom = 1;
    grid_points[1087].x =   0.1836309423; grid_points[1087].y =   0.1836309423; grid_points[1087].z =   4.1951916927; grid_points[1087].w_fixed =    0.0024069557; grid_points[1087].w_total = 0.0024069553; grid_points[1087].i_atom = 1;
    grid_points[1088].x =   0.1499340366; grid_points[1088].y =   0.1499340366; grid_points[1088].z =   4.3451257293; grid_points[1088].w_fixed =    0.0020308689; grid_points[1088].w_total = 0.0020308678; grid_points[1088].i_atom = 1;
    grid_points[1089].x =   0.1499340366; grid_points[1089].y =   0.1499340366; grid_points[1089].z =   4.0452576562; grid_points[1089].w_fixed =    0.0020308689; grid_points[1089].w_total = 0.0020308551; grid_points[1089].i_atom = 1;
    grid_points[1090].x =   0.2434922688; grid_points[1090].y =   0.2434922688; grid_points[1090].z =   4.1951916927; grid_points[1090].w_fixed =    0.0054297343; grid_points[1090].w_total = 0.0054297264; grid_points[1090].i_atom = 1;
    grid_points[1091].x =   0.1988106050; grid_points[1091].y =   0.1988106050; grid_points[1091].z =   4.3940022977; grid_points[1091].w_fixed =    0.0045813383; grid_points[1091].w_total = 0.0045813179; grid_points[1091].i_atom = 1;
    grid_points[1092].x =   0.1988106050; grid_points[1092].y =   0.1988106050; grid_points[1092].z =   3.9963810878; grid_points[1092].w_fixed =    0.0045813383; grid_points[1092].w_total = 0.0045810919; grid_points[1092].i_atom = 1;
    grid_points[1093].x =   0.3203734782; grid_points[1093].y =   0.3203734782; grid_points[1093].z =   4.1951916927; grid_points[1093].w_fixed =    0.0071632550; grid_points[1093].w_total = 0.0071631809; grid_points[1093].i_atom = 1;
    grid_points[1094].x =   0.2615838496; grid_points[1094].y =   0.2615838496; grid_points[1094].z =   4.4567755423; grid_points[1094].w_fixed =    0.0066932539; grid_points[1094].w_total = 0.0066930222; grid_points[1094].i_atom = 1;
    grid_points[1095].x =   0.2615838496; grid_points[1095].y =   0.2615838496; grid_points[1095].z =   3.9336078431; grid_points[1095].w_fixed =    0.0066932539; grid_points[1095].w_total = 0.0066907381; grid_points[1095].i_atom = 1;
    grid_points[1096].x =   0.1366077101; grid_points[1096].y =   0.1366077101; grid_points[1096].z =   4.6050148231; grid_points[1096].w_fixed =    0.0064011973; grid_points[1096].w_total = 0.0064000836; grid_points[1096].i_atom = 1;
    grid_points[1097].x =   0.1366077101; grid_points[1097].y =   0.1366077101; grid_points[1097].z =   3.7853685624; grid_points[1097].w_fixed =    0.0064011973; grid_points[1097].w_total = 0.0063905489; grid_points[1097].i_atom = 1;
    grid_points[1098].x =   0.1366077101; grid_points[1098].y =   0.4098231304; grid_points[1098].z =   4.3317994029; grid_points[1098].w_fixed =    0.0064011973; grid_points[1098].w_total = 0.0064011513; grid_points[1098].i_atom = 1;
    grid_points[1099].x =   0.1366077101; grid_points[1099].y =   0.4098231304; grid_points[1099].z =   4.0585839826; grid_points[1099].w_fixed =    0.0064011973; grid_points[1099].w_total = 0.0064006801; grid_points[1099].i_atom = 1;
    grid_points[1100].x =   0.4098231304; grid_points[1100].y =   0.1366077101; grid_points[1100].z =   4.3317994029; grid_points[1100].w_fixed =    0.0064011973; grid_points[1100].w_total = 0.0064011513; grid_points[1100].i_atom = 1;
    grid_points[1101].x =   0.4098231304; grid_points[1101].y =   0.1366077101; grid_points[1101].z =   4.0585839826; grid_points[1101].w_fixed =    0.0064011973; grid_points[1101].w_total = 0.0064006801; grid_points[1101].i_atom = 1;
    grid_points[1102].x =   0.4194534787; grid_points[1102].y =   0.4194534787; grid_points[1102].z =   4.1951916927; grid_points[1102].w_fixed =    0.0158612132; grid_points[1102].w_total = 0.0158601692; grid_points[1102].i_atom = 1;
    grid_points[1103].x =   0.3424823312; grid_points[1103].y =   0.3424823312; grid_points[1103].z =   4.5376740239; grid_points[1103].w_fixed =    0.0148205147; grid_points[1103].w_total = 0.0148168993; grid_points[1103].i_atom = 1;
    grid_points[1104].x =   0.3424823312; grid_points[1104].y =   0.3424823312; grid_points[1104].z =   3.8527093615; grid_points[1104].w_fixed =    0.0148205147; grid_points[1104].w_total = 0.0147863672; grid_points[1104].i_atom = 1;
    grid_points[1105].x =   0.1788555643; grid_points[1105].y =   0.1788555643; grid_points[1105].z =   4.7317583855; grid_points[1105].w_fixed =    0.0141738295; grid_points[1105].w_total = 0.0141564043; grid_points[1105].i_atom = 1;
    grid_points[1106].x =   0.1788555643; grid_points[1106].y =   0.1788555643; grid_points[1106].z =   3.6586249999; grid_points[1106].w_fixed =    0.0141738295; grid_points[1106].w_total = 0.0140351874; grid_points[1106].i_atom = 1;
    grid_points[1107].x =   0.1788555643; grid_points[1107].y =   0.5365666928; grid_points[1107].z =   4.3740472570; grid_points[1107].w_fixed =    0.0141738295; grid_points[1107].w_total = 0.0141731261; grid_points[1107].i_atom = 1;
    grid_points[1108].x =   0.1788555643; grid_points[1108].y =   0.5365666928; grid_points[1108].z =   4.0163361285; grid_points[1108].w_fixed =    0.0141738295; grid_points[1108].w_total = 0.0141666532; grid_points[1108].i_atom = 1;
    grid_points[1109].x =   0.5365666928; grid_points[1109].y =   0.1788555643; grid_points[1109].z =   4.3740472570; grid_points[1109].w_fixed =    0.0141738295; grid_points[1109].w_total = 0.0141731261; grid_points[1109].i_atom = 1;
    grid_points[1110].x =   0.5365666928; grid_points[1110].y =   0.1788555643; grid_points[1110].z =   4.0163361285; grid_points[1110].w_fixed =    0.0141738295; grid_points[1110].w_total = 0.0141666532; grid_points[1110].i_atom = 1;
    grid_points[1111].x =   0.4473238612; grid_points[1111].y =   0.4473238612; grid_points[1111].z =   4.6425155539; grid_points[1111].w_fixed =    0.0152642382; grid_points[1111].w_total = 0.0152405851; grid_points[1111].i_atom = 1;
    grid_points[1112].x =   0.4473238612; grid_points[1112].y =   0.4473238612; grid_points[1112].z =   3.7478678316; grid_points[1112].w_fixed =    0.0152642382; grid_points[1112].w_total = 0.0150785295; grid_points[1112].i_atom = 1;
    grid_points[1113].x =   0.1434253090; grid_points[1113].y =   0.1434253090; grid_points[1113].z =   4.9429578806; grid_points[1113].w_fixed =    0.0127985780; grid_points[1113].w_total = 0.0126709482; grid_points[1113].i_atom = 1;
    grid_points[1114].x =   0.1434253090; grid_points[1114].y =   0.1434253090; grid_points[1114].z =   3.4474255048; grid_points[1114].w_fixed =    0.0127985780; grid_points[1114].w_total = 0.0120391551; grid_points[1114].i_atom = 1;
    grid_points[1115].x =   0.1434253090; grid_points[1115].y =   0.7477661879; grid_points[1115].z =   4.3386170017; grid_points[1115].w_fixed =    0.0127985780; grid_points[1115].w_total = 0.0127960888; grid_points[1115].i_atom = 1;
    grid_points[1116].x =   0.1434253090; grid_points[1116].y =   0.7477661879; grid_points[1116].z =   4.0517663837; grid_points[1116].w_fixed =    0.0127985780; grid_points[1116].w_total = 0.0127814266; grid_points[1116].i_atom = 1;
    grid_points[1117].x =   0.7477661879; grid_points[1117].y =   0.1434253090; grid_points[1117].z =   4.3386170017; grid_points[1117].w_fixed =    0.0127985780; grid_points[1117].w_total = 0.0127960888; grid_points[1117].i_atom = 1;
    grid_points[1118].x =   0.7477661879; grid_points[1118].y =   0.1434253090; grid_points[1118].z =   4.0517663837; grid_points[1118].w_fixed =    0.0127985780; grid_points[1118].w_total = 0.0127814266; grid_points[1118].i_atom = 1;
    grid_points[1119].x =   0.5349297050; grid_points[1119].y =   0.5349297050; grid_points[1119].z =   4.3625127365; grid_points[1119].w_fixed =    0.0154965859; grid_points[1119].w_total = 0.0154933439; grid_points[1119].i_atom = 1;
    grid_points[1120].x =   0.5349297050; grid_points[1120].y =   0.5349297050; grid_points[1120].z =   4.0278706490; grid_points[1120].w_fixed =    0.0154965859; grid_points[1120].w_total = 0.0154711868; grid_points[1120].i_atom = 1;
    grid_points[1121].x =   0.5349297050; grid_points[1121].y =   0.1673210437; grid_points[1121].z =   4.7301213978; grid_points[1121].w_fixed =    0.0154965859; grid_points[1121].w_total = 0.0154531108; grid_points[1121].i_atom = 1;
    grid_points[1122].x =   0.5349297050; grid_points[1122].y =   0.1673210437; grid_points[1122].z =   3.6602619877; grid_points[1122].w_fixed =    0.0154965859; grid_points[1122].w_total = 0.0151815178; grid_points[1122].i_atom = 1;
    grid_points[1123].x =   0.1673210437; grid_points[1123].y =   0.5349297050; grid_points[1123].z =   4.7301213978; grid_points[1123].w_fixed =    0.0154965859; grid_points[1123].w_total = 0.0154531108; grid_points[1123].i_atom = 1;
    grid_points[1124].x =   0.1673210437; grid_points[1124].y =   0.5349297050; grid_points[1124].z =   3.6602619877; grid_points[1124].w_fixed =    0.0154965859; grid_points[1124].w_total = 0.0151815178; grid_points[1124].i_atom = 1;
    grid_points[1125].x =   0.3065753189; grid_points[1125].y =   0.3065753189; grid_points[1125].z =   4.8373124430; grid_points[1125].w_fixed =    0.0149552262; grid_points[1125].w_total = 0.0148734574; grid_points[1125].i_atom = 1;
    grid_points[1126].x =   0.3065753189; grid_points[1126].y =   0.3065753189; grid_points[1126].z =   3.5530709425; grid_points[1126].w_fixed =    0.0149552262; grid_points[1126].w_total = 0.0144183652; grid_points[1126].i_atom = 1;
    grid_points[1127].x =   0.3065753189; grid_points[1127].y =   0.6421207503; grid_points[1127].z =   4.5017670117; grid_points[1127].w_fixed =    0.0149552262; grid_points[1127].w_total = 0.0149471131; grid_points[1127].i_atom = 1;
    grid_points[1128].x =   0.3065753189; grid_points[1128].y =   0.6421207503; grid_points[1128].z =   3.8886163738; grid_points[1128].w_fixed =    0.0149552262; grid_points[1128].w_total = 0.0148835083; grid_points[1128].i_atom = 1;
    grid_points[1129].x =   0.6421207503; grid_points[1129].y =   0.3065753189; grid_points[1129].z =   4.5017670117; grid_points[1129].w_fixed =    0.0149552262; grid_points[1129].w_total = 0.0149471131; grid_points[1129].i_atom = 1;
    grid_points[1130].x =   0.6421207503; grid_points[1130].y =   0.3065753189; grid_points[1130].z =   3.8886163738; grid_points[1130].w_fixed =    0.0149552262; grid_points[1130].w_total = 0.0148835083; grid_points[1130].i_atom = 1;
    grid_points[1131].x =   0.3706344180; grid_points[1131].y =   0.6803866830; grid_points[1131].z =   4.1951916927; grid_points[1131].w_fixed =    0.0151103430; grid_points[1131].w_total = 0.0151046974; grid_points[1131].i_atom = 1;
    grid_points[1132].x =   0.6803866830; grid_points[1132].y =   0.3706344180; grid_points[1132].z =   4.1951916927; grid_points[1132].w_fixed =    0.0151103430; grid_points[1132].w_total = 0.0151046974; grid_points[1132].i_atom = 1;
    grid_points[1133].x =   0.5842597370; grid_points[1133].y =   0.5842597370; grid_points[1133].z =   4.7794514298; grid_points[1133].w_fixed =    0.0341634180; grid_points[1133].w_total = 0.0338660424; grid_points[1133].i_atom = 1;
    grid_points[1134].x =   0.5842597370; grid_points[1134].y =   0.5842597370; grid_points[1134].z =   3.6109319557; grid_points[1134].w_fixed =    0.0341634180; grid_points[1134].w_total = 0.0323366953; grid_points[1134].i_atom = 1;
    grid_points[1135].x =   0.1873310158; grid_points[1135].y =   0.1873310158; grid_points[1135].z =   5.1718658973; grid_points[1135].w_fixed =    0.0286449390; grid_points[1135].w_total = 0.0271394907; grid_points[1135].i_atom = 1;
    grid_points[1136].x =   0.1873310158; grid_points[1136].y =   0.1873310158; grid_points[1136].z =   3.2185174881; grid_points[1136].w_fixed =    0.0286449390; grid_points[1136].w_total = 0.0223716243; grid_points[1136].i_atom = 1;
    grid_points[1137].x =   0.1873310158; grid_points[1137].y =   0.9766742046; grid_points[1137].z =   4.3825227086; grid_points[1137].w_fixed =    0.0286449390; grid_points[1137].w_total = 0.0286147757; grid_points[1137].i_atom = 1;
    grid_points[1138].x =   0.1873310158; grid_points[1138].y =   0.9766742046; grid_points[1138].z =   4.0078606769; grid_points[1138].w_fixed =    0.0286449390; grid_points[1138].w_total = 0.0284587127; grid_points[1138].i_atom = 1;
    grid_points[1139].x =   0.9766742046; grid_points[1139].y =   0.1873310158; grid_points[1139].z =   4.3825227086; grid_points[1139].w_fixed =    0.0286449390; grid_points[1139].w_total = 0.0286147757; grid_points[1139].i_atom = 1;
    grid_points[1140].x =   0.9766742046; grid_points[1140].y =   0.1873310158; grid_points[1140].z =   4.0078606769; grid_points[1140].w_fixed =    0.0286449390; grid_points[1140].w_total = 0.0284587127; grid_points[1140].i_atom = 1;
    grid_points[1141].x =   0.6986836963; grid_points[1141].y =   0.6986836963; grid_points[1141].z =   4.4137334641; grid_points[1141].w_fixed =    0.0346834435; grid_points[1141].w_total = 0.0346436110; grid_points[1141].i_atom = 1;
    grid_points[1142].x =   0.6986836963; grid_points[1142].y =   0.6986836963; grid_points[1142].z =   3.9766499214; grid_points[1142].w_fixed =    0.0346834435; grid_points[1142].w_total = 0.0344092150; grid_points[1142].i_atom = 1;
    grid_points[1143].x =   0.6986836963; grid_points[1143].y =   0.2185417714; grid_points[1143].z =   4.8938753891; grid_points[1143].w_fixed =    0.0346834435; grid_points[1143].w_total = 0.0341432423; grid_points[1143].i_atom = 1;
    grid_points[1144].x =   0.6986836963; grid_points[1144].y =   0.2185417714; grid_points[1144].z =   3.4965079964; grid_points[1144].w_fixed =    0.0346834435; grid_points[1144].w_total = 0.0317117114; grid_points[1144].i_atom = 1;
    grid_points[1145].x =   0.2185417714; grid_points[1145].y =   0.6986836963; grid_points[1145].z =   4.8938753891; grid_points[1145].w_fixed =    0.0346834435; grid_points[1145].w_total = 0.0341432423; grid_points[1145].i_atom = 1;
    grid_points[1146].x =   0.2185417714; grid_points[1146].y =   0.6986836963; grid_points[1146].z =   3.4965079964; grid_points[1146].w_fixed =    0.0346834435; grid_points[1146].w_total = 0.0317117114; grid_points[1146].i_atom = 1;
    grid_points[1147].x =   0.4004249063; grid_points[1147].y =   0.4004249063; grid_points[1147].z =   5.0338800196; grid_points[1147].w_fixed =    0.0334718077; grid_points[1147].w_total = 0.0324778009; grid_points[1147].i_atom = 1;
    grid_points[1148].x =   0.4004249063; grid_points[1148].y =   0.4004249063; grid_points[1148].z =   3.3565033658; grid_points[1148].w_fixed =    0.0334718077; grid_points[1148].w_total = 0.0287086336; grid_points[1148].i_atom = 1;
    grid_points[1149].x =   0.4004249063; grid_points[1149].y =   0.8386883269; grid_points[1149].z =   4.5956165991; grid_points[1149].w_fixed =    0.0334718077; grid_points[1149].w_total = 0.0333691081; grid_points[1149].i_atom = 1;
    grid_points[1150].x =   0.4004249063; grid_points[1150].y =   0.8386883269; grid_points[1150].z =   3.7947667864; grid_points[1150].w_fixed =    0.0334718077; grid_points[1150].w_total = 0.0327270168; grid_points[1150].i_atom = 1;
    grid_points[1151].x =   0.8386883269; grid_points[1151].y =   0.4004249063; grid_points[1151].z =   4.5956165991; grid_points[1151].w_fixed =    0.0334718077; grid_points[1151].w_total = 0.0333691081; grid_points[1151].i_atom = 1;
    grid_points[1152].x =   0.8386883269; grid_points[1152].y =   0.4004249063; grid_points[1152].z =   3.7947667864; grid_points[1152].w_fixed =    0.0334718077; grid_points[1152].w_total = 0.0327270168; grid_points[1152].i_atom = 1;
    grid_points[1153].x =   0.4840939338; grid_points[1153].y =   0.8886683207; grid_points[1153].z =   4.1951916927; grid_points[1153].w_fixed =    0.0338189800; grid_points[1153].w_total = 0.0337557378; grid_points[1153].i_atom = 1;
    grid_points[1154].x =   0.8886683207; grid_points[1154].y =   0.4840939338; grid_points[1154].z =   4.1951916927; grid_points[1154].w_fixed =    0.0338189800; grid_points[1154].w_total = 0.0337557378; grid_points[1154].i_atom = 1;
    grid_points[1155].x =   0.7649502715; grid_points[1155].y =   0.7649502715; grid_points[1155].z =   4.9601419642; grid_points[1155].w_fixed =    0.0777139157; grid_points[1155].w_total = 0.0744471256; grid_points[1155].i_atom = 1;
    grid_points[1156].x =   0.7649502715; grid_points[1156].y =   0.7649502715; grid_points[1156].z =   3.4302414213; grid_points[1156].w_fixed =    0.0777139157; grid_points[1156].w_total = 0.0632755179; grid_points[1156].i_atom = 1;
    grid_points[1157].x =   0.2452657651; grid_points[1157].y =   0.2452657651; grid_points[1157].z =   5.4739161206; grid_points[1157].w_fixed =    0.0651606456; grid_points[1157].w_total = 0.0509705226; grid_points[1157].i_atom = 1;
    grid_points[1158].x =   0.2452657651; grid_points[1158].y =   1.2787244278; grid_points[1158].z =   4.4404574578; grid_points[1158].w_fixed =    0.0651606456; grid_points[1158].w_total = 0.0648319178; grid_points[1158].i_atom = 1;
    grid_points[1159].x =   0.2452657651; grid_points[1159].y =   1.2787244278; grid_points[1159].z =   3.9499259276; grid_points[1159].w_fixed =    0.0651606456; grid_points[1159].w_total = 0.0634126219; grid_points[1159].i_atom = 1;
    grid_points[1160].x =   1.2787244278; grid_points[1160].y =   0.2452657651; grid_points[1160].z =   4.4404574578; grid_points[1160].w_fixed =    0.0651606456; grid_points[1160].w_total = 0.0648319178; grid_points[1160].i_atom = 1;
    grid_points[1161].x =   1.2787244278; grid_points[1161].y =   0.2452657651; grid_points[1161].z =   3.9499259276; grid_points[1161].w_fixed =    0.0651606456; grid_points[1161].w_total = 0.0634126219; grid_points[1161].i_atom = 1;
    grid_points[1162].x =   0.9147614482; grid_points[1162].y =   0.9147614482; grid_points[1162].z =   4.4813205782; grid_points[1162].w_fixed =    0.0788968540; grid_points[1162].w_total = 0.0784571533; grid_points[1162].i_atom = 1;
    grid_points[1163].x =   0.9147614482; grid_points[1163].y =   0.9147614482; grid_points[1163].z =   3.9090628073; grid_points[1163].w_fixed =    0.0788968540; grid_points[1163].w_total = 0.0763502659; grid_points[1163].i_atom = 1;
    grid_points[1164].x =   0.9147614482; grid_points[1164].y =   0.2861288854; grid_points[1164].z =   5.1099531410; grid_points[1164].w_fixed =    0.0788968540; grid_points[1164].w_total = 0.0731312984; grid_points[1164].i_atom = 1;
    grid_points[1165].x =   0.9147614482; grid_points[1165].y =   0.2861288854; grid_points[1165].z =   3.2804302445; grid_points[1165].w_fixed =    0.0788968540; grid_points[1165].w_total = 0.0571063680; grid_points[1165].i_atom = 1;
    grid_points[1166].x =   0.2861288854; grid_points[1166].y =   0.9147614482; grid_points[1166].z =   5.1099531410; grid_points[1166].w_fixed =    0.0788968540; grid_points[1166].w_total = 0.0731312984; grid_points[1166].i_atom = 1;
    grid_points[1167].x =   0.2861288854; grid_points[1167].y =   0.9147614482; grid_points[1167].z =   3.2804302445; grid_points[1167].w_fixed =    0.0788968540; grid_points[1167].w_total = 0.0571063680; grid_points[1167].i_atom = 1;
    grid_points[1168].x =   0.5242619359; grid_points[1168].y =   0.5242619359; grid_points[1168].z =   5.2932561710; grid_points[1168].w_fixed =    0.0761406613; grid_points[1168].w_total = 0.0660727372; grid_points[1168].i_atom = 1;
    grid_points[1169].x =   0.5242619359; grid_points[1169].y =   1.0980644783; grid_points[1169].z =   4.7194536287; grid_points[1169].w_fixed =    0.0761406613; grid_points[1169].w_total = 0.0749832017; grid_points[1169].i_atom = 1;
    grid_points[1170].x =   0.5242619359; grid_points[1170].y =   1.0980644783; grid_points[1170].z =   3.6709297568; grid_points[1170].w_fixed =    0.0761406613; grid_points[1170].w_total = 0.0696636336; grid_points[1170].i_atom = 1;
    grid_points[1171].x =   1.0980644783; grid_points[1171].y =   0.5242619359; grid_points[1171].z =   4.7194536287; grid_points[1171].w_fixed =    0.0761406613; grid_points[1171].w_total = 0.0749832017; grid_points[1171].i_atom = 1;
    grid_points[1172].x =   1.0980644783; grid_points[1172].y =   0.5242619359; grid_points[1172].z =   3.6709297568; grid_points[1172].w_fixed =    0.0761406613; grid_points[1172].w_total = 0.0696636336; grid_points[1172].i_atom = 1;
    grid_points[1173].x =   0.6338067859; grid_points[1173].y =   1.1635014876; grid_points[1173].z =   4.1951916927; grid_points[1173].w_fixed =    0.0769303986; grid_points[1173].w_total = 0.0763061913; grid_points[1173].i_atom = 1;
    grid_points[1174].x =   1.1635014876; grid_points[1174].y =   0.6338067859; grid_points[1174].z =   4.1951916927; grid_points[1174].w_fixed =    0.0769303986; grid_points[1174].w_total = 0.0763061913; grid_points[1174].i_atom = 1;
    grid_points[1175].x =   1.0064786876; grid_points[1175].y =   1.0064786876; grid_points[1175].z =   5.2016703804; grid_points[1175].w_fixed =    0.1811137632; grid_points[1175].w_total = 0.1504983152; grid_points[1175].i_atom = 1;
    grid_points[1176].x =   1.0064786876; grid_points[1176].y =   1.0064786876; grid_points[1176].z =   3.1887130051; grid_points[1176].w_fixed =    0.1811137632; grid_points[1176].w_total = 0.0948100401; grid_points[1176].i_atom = 1;
    grid_points[1177].x =   0.3227069453; grid_points[1177].y =   0.3227069453; grid_points[1177].z =   5.8776656155; grid_points[1177].w_fixed =    0.1518581277; grid_points[1177].w_total = 0.0603864433; grid_points[1177].i_atom = 1;
    grid_points[1178].x =   0.3227069453; grid_points[1178].y =   1.6824739228; grid_points[1178].z =   4.5178986380; grid_points[1178].w_fixed =    0.1518581277; grid_points[1178].w_total = 0.1485902764; grid_points[1178].i_atom = 1;
    grid_points[1179].x =   0.3227069453; grid_points[1179].y =   1.6824739228; grid_points[1179].z =   3.8724847475; grid_points[1179].w_fixed =    0.1518581277; grid_points[1179].w_total = 0.1376185851; grid_points[1179].i_atom = 1;
    grid_points[1180].x =   1.6824739228; grid_points[1180].y =   0.3227069453; grid_points[1180].z =   4.5178986380; grid_points[1180].w_fixed =    0.1518581277; grid_points[1180].w_total = 0.1485902764; grid_points[1180].i_atom = 1;
    grid_points[1181].x =   1.6824739228; grid_points[1181].y =   0.3227069453; grid_points[1181].z =   3.8724847475; grid_points[1181].w_fixed =    0.1518581277; grid_points[1181].w_total = 0.1376185851; grid_points[1181].i_atom = 1;
    grid_points[1182].x =   1.2035918363; grid_points[1182].y =   1.2035918363; grid_points[1182].z =   4.5716640411; grid_points[1182].w_fixed =    0.1838706234; grid_points[1182].w_total = 0.1794730010; grid_points[1182].i_atom = 1;
    grid_points[1183].x =   1.2035918363; grid_points[1183].y =   1.2035918363; grid_points[1183].z =   3.8187193444; grid_points[1183].w_fixed =    0.1838706234; grid_points[1183].w_total = 0.1635857757; grid_points[1183].i_atom = 1;
    grid_points[1184].x =   1.2035918363; grid_points[1184].y =   0.3764723484; grid_points[1184].z =   5.3987835290; grid_points[1184].w_fixed =    0.1838706234; grid_points[1184].w_total = 0.1337193515; grid_points[1184].i_atom = 1;
    grid_points[1185].x =   0.3764723484; grid_points[1185].y =   1.2035918363; grid_points[1185].z =   5.3987835290; grid_points[1185].w_fixed =    0.1838706234; grid_points[1185].w_total = 0.1337193515; grid_points[1185].i_atom = 1;
    grid_points[1186].x =   0.6897944676; grid_points[1186].y =   0.6897944676; grid_points[1186].z =   5.6399633809; grid_points[1186].w_fixed =    0.1774472638; grid_points[1186].w_total = 0.1007712847; grid_points[1186].i_atom = 1;
    grid_points[1187].x =   0.6897944676; grid_points[1187].y =   1.4447716881; grid_points[1187].z =   4.8849861603; grid_points[1187].w_fixed =    0.1774472638; grid_points[1187].w_total = 0.1658131153; grid_points[1187].i_atom = 1;
    grid_points[1188].x =   0.6897944676; grid_points[1188].y =   1.4447716881; grid_points[1188].z =   3.5053972252; grid_points[1188].w_fixed =    0.1774472638; grid_points[1188].w_total = 0.1314016342; grid_points[1188].i_atom = 1;
    grid_points[1189].x =   1.4447716881; grid_points[1189].y =   0.6897944676; grid_points[1189].z =   4.8849861603; grid_points[1189].w_fixed =    0.1774472638; grid_points[1189].w_total = 0.1658131153; grid_points[1189].i_atom = 1;
    grid_points[1190].x =   1.4447716881; grid_points[1190].y =   0.6897944676; grid_points[1190].z =   3.5053972252; grid_points[1190].w_fixed =    0.1774472638; grid_points[1190].w_total = 0.1314016342; grid_points[1190].i_atom = 1;
    grid_points[1191].x =   0.8339274406; grid_points[1191].y =   1.5308700368; grid_points[1191].z =   4.1951916927; grid_points[1191].w_fixed =    0.1792877615; grid_points[1191].w_total = 0.1738017119; grid_points[1191].i_atom = 1;
    grid_points[1192].x =   1.5308700368; grid_points[1192].y =   0.8339274406; grid_points[1192].z =   4.1951916927; grid_points[1192].w_fixed =    0.1792877615; grid_points[1192].w_total = 0.1738017119; grid_points[1192].i_atom = 1;
    grid_points[1193].x =   1.3345778006; grid_points[1193].y =   1.3345778006; grid_points[1193].z =   5.5297694934; grid_points[1193].w_fixed =    0.4363920800; grid_points[1193].w_total = 0.2072585535; grid_points[1193].i_atom = 1;
    grid_points[1194].x =   0.4279052608; grid_points[1194].y =   2.2309387920; grid_points[1194].z =   4.6230969535; grid_points[1194].w_fixed =    0.3659008737; grid_points[1194].w_total = 0.3324753846; grid_points[1194].i_atom = 1;
    grid_points[1195].x =   0.4279052608; grid_points[1195].y =   2.2309387920; grid_points[1195].z =   3.7672864320; grid_points[1195].w_fixed =    0.3659008737; grid_points[1195].w_total = 0.2578929044; grid_points[1195].i_atom = 1;
    grid_points[1196].x =   2.2309387920; grid_points[1196].y =   0.4279052608; grid_points[1196].z =   4.6230969535; grid_points[1196].w_fixed =    0.3659008737; grid_points[1196].w_total = 0.3324753846; grid_points[1196].i_atom = 1;
    grid_points[1197].x =   2.2309387920; grid_points[1197].y =   0.4279052608; grid_points[1197].z =   3.7672864320; grid_points[1197].w_fixed =    0.3659008737; grid_points[1197].w_total = 0.2578929044; grid_points[1197].i_atom = 1;
    grid_points[1198].x =   1.5959473017; grid_points[1198].y =   1.5959473017; grid_points[1198].z =   4.6943891868; grid_points[1198].w_fixed =    0.4430347113; grid_points[1198].w_total = 0.3987026542; grid_points[1198].i_atom = 1;
    grid_points[1199].x =   1.5959473017; grid_points[1199].y =   1.5959473017; grid_points[1199].z =   3.6959941987; grid_points[1199].w_fixed =    0.4430347113; grid_points[1199].w_total = 0.2967672396; grid_points[1199].i_atom = 1;
    grid_points[1200].x =   1.5959473017; grid_points[1200].y =   0.4991974941; grid_points[1200].z =   5.7911389945; grid_points[1200].w_fixed =    0.4430347113; grid_points[1200].w_total = 0.1354214943; grid_points[1200].i_atom = 1;
    grid_points[1201].x =   0.4991974941; grid_points[1201].y =   1.5959473017; grid_points[1201].z =   5.7911389945; grid_points[1201].w_fixed =    0.4430347113; grid_points[1201].w_total = 0.1354214943; grid_points[1201].i_atom = 1;
    grid_points[1202].x =   0.9146585961; grid_points[1202].y =   0.9146585961; grid_points[1202].z =   6.1109403775; grid_points[1202].w_fixed =    0.4275576805; grid_points[1202].w_total = 0.0649370027; grid_points[1202].i_atom = 1;
    grid_points[1203].x =   0.9146585961; grid_points[1203].y =   1.9157486847; grid_points[1203].z =   5.1098502889; grid_points[1203].w_fixed =    0.4275576805; grid_points[1203].w_total = 0.3174678581; grid_points[1203].i_atom = 1;
    grid_points[1204].x =   0.9146585961; grid_points[1204].y =   1.9157486847; grid_points[1204].z =   3.2805330966; grid_points[1204].w_fixed =    0.4275576805; grid_points[1204].w_total = 0.1674035783; grid_points[1204].i_atom = 1;
    grid_points[1205].x =   1.9157486847; grid_points[1205].y =   0.9146585961; grid_points[1205].z =   5.1098502889; grid_points[1205].w_fixed =    0.4275576805; grid_points[1205].w_total = 0.3174678581; grid_points[1205].i_atom = 1;
    grid_points[1206].x =   1.9157486847; grid_points[1206].y =   0.9146585961; grid_points[1206].z =   3.2805330966; grid_points[1206].w_fixed =    0.4275576805; grid_points[1206].w_total = 0.1674035783; grid_points[1206].i_atom = 1;
    grid_points[1207].x =   1.1057770654; grid_points[1207].y =   2.0299139882; grid_points[1207].z =   4.1951916927; grid_points[1207].w_fixed =    0.4319923443; grid_points[1207].w_total = 0.3838311380; grid_points[1207].i_atom = 1;
    grid_points[1208].x =   2.0299139882; grid_points[1208].y =   1.1057770654; grid_points[1208].z =   4.1951916927; grid_points[1208].w_fixed =    0.4319923443; grid_points[1208].w_total = 0.3838311380; grid_points[1208].i_atom = 1;
    grid_points[1209].x =   1.7892954446; grid_points[1209].y =   1.7892954446; grid_points[1209].z =   5.9844871374; grid_points[1209].w_fixed =    1.0990251469; grid_points[1209].w_total = 0.0671861197; grid_points[1209].i_atom = 1;
    grid_points[1210].x =   2.1397188201; grid_points[1210].y =   2.1397188201; grid_points[1210].z =   4.8644758676; grid_points[1210].w_fixed =    1.1157541830; grid_points[1210].w_total = 0.5622597816; grid_points[1210].i_atom = 1;
    grid_points[1211].x =   2.1397188201; grid_points[1211].y =   2.1397188201; grid_points[1211].z =   3.5259075179; grid_points[1211].w_fixed =    1.1157541830; grid_points[1211].w_total = 0.2172616206; grid_points[1211].i_atom = 1;
    grid_points[1212].x =   1.2263012757; grid_points[1212].y =   2.5684830011; grid_points[1212].z =   5.4214929684; grid_points[1212].w_fixed =    1.0767762850; grid_points[1212].w_total = 0.2278386018; grid_points[1212].i_atom = 1;
    grid_points[1213].x =   2.5684830011; grid_points[1213].y =   1.2263012757; grid_points[1213].z =   5.4214929684; grid_points[1213].w_fixed =    1.0767762850; grid_points[1213].w_total = 0.2278386018; grid_points[1213].i_atom = 1;
    grid_points[1214].x =   1.4825376722; grid_points[1214].y =   2.7215467321; grid_points[1214].z =   4.1951916927; grid_points[1214].w_fixed =    1.0879446980; grid_points[1214].w_total = 0.5523500480; grid_points[1214].i_atom = 1;
    grid_points[1215].x =   2.7215467321; grid_points[1215].y =   1.4825376722; grid_points[1215].z =   4.1951916927; grid_points[1215].w_fixed =    1.0879446980; grid_points[1215].w_total = 0.5523500480; grid_points[1215].i_atom = 1;
    grid_points[1216].x =   0.5737012360; grid_points[1216].y =   2.9910647516; grid_points[1216].z =   4.7688929288; grid_points[1216].w_fixed =    0.9214976165; grid_points[1216].w_total = 0.4792995334; grid_points[1216].i_atom = 1;
    grid_points[1217].x =   0.5737012360; grid_points[1217].y =   2.9910647516; grid_points[1217].z =   3.6214904567; grid_points[1217].w_fixed =    0.9214976165; grid_points[1217].w_total = 0.1977834187; grid_points[1217].i_atom = 1;
    grid_points[1218].x =   2.9910647516; grid_points[1218].y =   0.5737012360; grid_points[1218].z =   4.7688929288; grid_points[1218].w_fixed =    0.9214976165; grid_points[1218].w_total = 0.4792995334; grid_points[1218].i_atom = 1;
    grid_points[1219].x =   2.9910647516; grid_points[1219].y =   0.5737012360; grid_points[1219].z =   3.6214904567; grid_points[1219].w_fixed =    0.9214976165; grid_points[1219].w_total = 0.1977834187; grid_points[1219].i_atom = 1;
    grid_points[1220].x =   2.9827802928; grid_points[1220].y =   2.9827802928; grid_points[1220].z =   4.1951916927; grid_points[1220].w_fixed =    6.7598180997; grid_points[1220].w_total = 0.0643528246; grid_points[1220].i_atom = 1;
    grid_points[1221].x =  -7.5660149734; grid_points[1221].y =   2.5220049911; grid_points[1221].z =   6.7171966839; grid_points[1221].w_fixed =   55.2891360098; grid_points[1221].w_total = 0.0000000000; grid_points[1221].i_atom = 1;
    grid_points[1222].x =   2.5220049911; grid_points[1222].y =  -7.5660149734; grid_points[1222].z =   6.7171966839; grid_points[1222].w_fixed =   55.2891360098; grid_points[1222].w_total = 0.0000000000; grid_points[1222].i_atom = 1;
    grid_points[1223].x =   1.7666549613; grid_points[1223].y =   1.7666549613; grid_points[1223].z =  -1.1047731912; grid_points[1223].w_fixed =   17.3848280115; grid_points[1223].w_total = 0.0000000000; grid_points[1223].i_atom = 1;
    grid_points[1224].x =   0.4279052608; grid_points[1224].y =   0.4279052608; grid_points[1224].z =   6.4261304847; grid_points[1224].w_fixed =    0.3659008737; grid_points[1224].w_total = 0.0221696713; grid_points[1224].i_atom = 1;
    grid_points[1225].x =   2.1397188201; grid_points[1225].y =   0.6692841749; grid_points[1225].z =   6.3349105128; grid_points[1225].w_fixed =    1.1157541830; grid_points[1225].w_total = 0.0256935379; grid_points[1225].i_atom = 1;
    grid_points[1226].x =   0.6692841749; grid_points[1226].y =   2.1397188201; grid_points[1226].z =   6.3349105128; grid_points[1226].w_fixed =    1.1157541830; grid_points[1226].w_total = 0.0256935379; grid_points[1226].i_atom = 1;
    grid_points[1227].x =   2.4354299108; grid_points[1227].y =   2.4354299108; grid_points[1227].z =   6.6306216035; grid_points[1227].w_fixed =    6.3162875542; grid_points[1227].w_total = 0.0011258350; grid_points[1227].i_atom = 1;
    grid_points[1228].x =   0.5737012360; grid_points[1228].y =   0.5737012360; grid_points[1228].z =   7.1862564443; grid_points[1228].w_fixed =    0.9214976165; grid_points[1228].w_total = 0.0005834114; grid_points[1228].i_atom = 1;
    grid_points[1229].x =   1.2263012757; grid_points[1229].y =   1.2263012757; grid_points[1229].z =   6.7636746939; grid_points[1229].w_fixed =    1.0767762850; grid_points[1229].w_total = 0.0058282274; grid_points[1229].i_atom = 1;
    grid_points[1230].x =   1.2718617903; grid_points[1230].y =   1.2718617903; grid_points[1230].z =   8.0107770637; grid_points[1230].w_fixed =    6.0406797362; grid_points[1230].w_total = 0.0000010065; grid_points[1230].i_atom = 1;
    grid_points[1231].x =   1.7666549613; grid_points[1231].y =  -5.2999648839; grid_points[1231].z =   2.4285367314; grid_points[1231].w_fixed =   17.3848280115; grid_points[1231].w_total = 0.0000000010; grid_points[1231].i_atom = 1;
    grid_points[1232].x =   0.0000000000; grid_points[1232].y =  -4.1431731363; grid_points[1232].z =   0.0520185564; grid_points[1232].w_fixed =   19.4544786654; grid_points[1232].w_total = 0.0000000000; grid_points[1232].i_atom = 1;
    grid_points[1233].x =   1.2718617903; grid_points[1233].y =  -3.8155853710; grid_points[1233].z =   2.9233299024; grid_points[1233].w_fixed =    6.0406797362; grid_points[1233].w_total = 0.0032694313; grid_points[1233].i_atom = 1;
    grid_points[1234].x =   1.7666549613; grid_points[1234].y =   5.2999648839; grid_points[1234].z =   2.4285367314; grid_points[1234].w_fixed =   17.3848280115; grid_points[1234].w_total = 0.0000000010; grid_points[1234].i_atom = 1;
    grid_points[1235].x =   1.2718617903; grid_points[1235].y =   3.8155853710; grid_points[1235].z =   2.9233299024; grid_points[1235].w_fixed =    6.0406797362; grid_points[1235].w_total = 0.0032694313; grid_points[1235].i_atom = 1;
    grid_points[1236].x =   0.0000000000; grid_points[1236].y =  -5.8593316407; grid_points[1236].z =   4.1951916927; grid_points[1236].w_fixed =   10.9431442493; grid_points[1236].w_total = 0.0000000054; grid_points[1236].i_atom = 1;
    grid_points[1237].x =   1.7666549613; grid_points[1237].y =  -5.2999648839; grid_points[1237].z =   5.9618466541; grid_points[1237].w_fixed =   17.3848280115; grid_points[1237].w_total = 0.0000000058; grid_points[1237].i_atom = 1;
    grid_points[1238].x =   0.0000000000; grid_points[1238].y =  -4.2182883437; grid_points[1238].z =   4.1951916927; grid_points[1238].w_fixed =    3.8023976811; grid_points[1238].w_total = 0.0216466137; grid_points[1238].i_atom = 1;
    grid_points[1239].x =   1.2718617903; grid_points[1239].y =  -3.8155853710; grid_points[1239].z =   5.4670534831; grid_points[1239].w_fixed =    6.0406797362; grid_points[1239].w_total = 0.0167747560; grid_points[1239].i_atom = 1;
    grid_points[1240].x =   1.2718617903; grid_points[1240].y =   3.8155853710; grid_points[1240].z =   5.4670534831; grid_points[1240].w_fixed =    6.0406797362; grid_points[1240].w_total = 0.0167747560; grid_points[1240].i_atom = 1;
    grid_points[1241].x =   1.7666549613; grid_points[1241].y =   5.2999648839; grid_points[1241].z =   5.9618466541; grid_points[1241].w_fixed =   17.3848280115; grid_points[1241].w_total = 0.0000000058; grid_points[1241].i_atom = 1;
    grid_points[1242].x =   0.0000000000; grid_points[1242].y =  -5.9146259784; grid_points[1242].z =  -1.7194342857; grid_points[1242].w_fixed =   61.8712659234; grid_points[1242].w_total = 0.0000000000; grid_points[1242].i_atom = 1;
    grid_points[1243].x =  -7.5660149734; grid_points[1243].y =  -2.5220049911; grid_points[1243].z =   1.6731867016; grid_points[1243].w_fixed =   55.2891360098; grid_points[1243].w_total = 0.0000000000; grid_points[1243].i_atom = 1;
    grid_points[1244].x =   0.4991974941; grid_points[1244].y =  -1.5959473017; grid_points[1244].z =   2.5992443910; grid_points[1244].w_fixed =    0.4430347113; grid_points[1244].w_total = 0.0354740456; grid_points[1244].i_atom = 1;
    grid_points[1245].x =   1.7892954446; grid_points[1245].y =  -1.7892954446; grid_points[1245].z =   2.4058962481; grid_points[1245].w_fixed =    1.0990251469; grid_points[1245].w_total = 0.0123810957; grid_points[1245].i_atom = 1;
    grid_points[1246].x =   0.6692841749; grid_points[1246].y =  -2.1397188201; grid_points[1246].z =   2.0554728727; grid_points[1246].w_fixed =    1.1157541830; grid_points[1246].w_total = 0.0031940004; grid_points[1246].i_atom = 1;
    grid_points[1247].x =   0.0000000000; grid_points[1247].y =  -1.4825376722; grid_points[1247].z =   1.4736449606; grid_points[1247].w_fixed =    1.0879446980; grid_points[1247].w_total = 0.0000849175; grid_points[1247].i_atom = 1;
    grid_points[1248].x =   0.0000000000; grid_points[1248].y =  -2.7215467321; grid_points[1248].z =   2.7126540206; grid_points[1248].w_fixed =    1.0879446980; grid_points[1248].w_total = 0.0291805697; grid_points[1248].i_atom = 1;
    grid_points[1249].x =   0.0000000000; grid_points[1249].y =  -2.9827802928; grid_points[1249].z =   1.2124113999; grid_points[1249].w_fixed =    6.7598180997; grid_points[1249].w_total = 0.0000096089; grid_points[1249].i_atom = 1;
    grid_points[1250].x =   2.4354299108; grid_points[1250].y =  -2.4354299108; grid_points[1250].z =   1.7597617820; grid_points[1250].w_fixed =    6.3162875542; grid_points[1250].w_total = 0.0001328651; grid_points[1250].i_atom = 1;
    grid_points[1251].x =   1.2718617903; grid_points[1251].y =  -1.2718617903; grid_points[1251].z =   0.3796063218; grid_points[1251].w_fixed =    6.0406797362; grid_points[1251].w_total = 0.0000000021; grid_points[1251].i_atom = 1;
    grid_points[1252].x =   1.3345778006; grid_points[1252].y =  -1.3345778006; grid_points[1252].z =   2.8606138921; grid_points[1252].w_fixed =    0.4363920800; grid_points[1252].w_total = 0.0730945598; grid_points[1252].i_atom = 1;
    grid_points[1253].x =   1.2263012757; grid_points[1253].y =  -2.5684830011; grid_points[1253].z =   2.9688904170; grid_points[1253].w_fixed =    1.0767762850; grid_points[1253].w_total = 0.0592143299; grid_points[1253].i_atom = 1;
    grid_points[1254].x =   1.7666549613; grid_points[1254].y =  -1.7666549613; grid_points[1254].z =   9.4951565767; grid_points[1254].w_fixed =   17.3848280115; grid_points[1254].w_total = 0.0000000000; grid_points[1254].i_atom = 1;
    grid_points[1255].x =   2.5220049911; grid_points[1255].y =  -2.5220049911; grid_points[1255].z =  11.7612066662; grid_points[1255].w_fixed =   55.2891360098; grid_points[1255].w_total = 0.0000000000; grid_points[1255].i_atom = 1;
    grid_points[1256].x =   0.3227069453; grid_points[1256].y =  -0.3227069453; grid_points[1256].z =   2.5127177700; grid_points[1256].w_fixed =    0.1518581277; grid_points[1256].w_total = 0.0160452730; grid_points[1256].i_atom = 1;
    grid_points[1257].x =   0.8339274406; grid_points[1257].y =   0.0000000000; grid_points[1257].z =   2.6643216559; grid_points[1257].w_fixed =    0.1792877615; grid_points[1257].w_total = 0.0310706777; grid_points[1257].i_atom = 1;
    grid_points[1258].x =   0.4279052608; grid_points[1258].y =  -0.4279052608; grid_points[1258].z =   1.9642529008; grid_points[1258].w_fixed =    0.3659008737; grid_points[1258].w_total = 0.0015586987; grid_points[1258].i_atom = 1;
    grid_points[1259].x =   1.5959473017; grid_points[1259].y =  -0.4991974941; grid_points[1259].z =   2.5992443910; grid_points[1259].w_fixed =    0.4430347113; grid_points[1259].w_total = 0.0354740456; grid_points[1259].i_atom = 1;
    grid_points[1260].x =   0.9146585961; grid_points[1260].y =  -0.9146585961; grid_points[1260].z =   2.2794430080; grid_points[1260].w_fixed =    0.4275576805; grid_points[1260].w_total = 0.0104656786; grid_points[1260].i_atom = 1;
    grid_points[1261].x =   1.1057770654; grid_points[1261].y =   0.0000000000; grid_points[1261].z =   2.1652777045; grid_points[1261].w_fixed =    0.4319923443; grid_points[1261].w_total = 0.0061142765; grid_points[1261].i_atom = 1;
    grid_points[1262].x =   0.5737012360; grid_points[1262].y =  -0.5737012360; grid_points[1262].z =   1.2041269412; grid_points[1262].w_fixed =    0.9214976165; grid_points[1262].w_total = 0.0000026829; grid_points[1262].i_atom = 1;
    grid_points[1263].x =   2.1397188201; grid_points[1263].y =  -0.6692841749; grid_points[1263].z =   2.0554728727; grid_points[1263].w_fixed =    1.1157541830; grid_points[1263].w_total = 0.0031940004; grid_points[1263].i_atom = 1;
    grid_points[1264].x =   1.2263012757; grid_points[1264].y =  -1.2263012757; grid_points[1264].z =   1.6267086916; grid_points[1264].w_fixed =    1.0767762850; grid_points[1264].w_total = 0.0002888020; grid_points[1264].i_atom = 1;
    grid_points[1265].x =   1.4825376722; grid_points[1265].y =   0.0000000000; grid_points[1265].z =   1.4736449606; grid_points[1265].w_fixed =    1.0879446980; grid_points[1265].w_total = 0.0000849175; grid_points[1265].i_atom = 1;
    grid_points[1266].x =   2.7215467321; grid_points[1266].y =   0.0000000000; grid_points[1266].z =   2.7126540206; grid_points[1266].w_fixed =    1.0879446980; grid_points[1266].w_total = 0.0291805697; grid_points[1266].i_atom = 1;
    grid_points[1267].x =   0.2452657651; grid_points[1267].y =  -0.2452657651; grid_points[1267].z =   2.9164672649; grid_points[1267].w_fixed =    0.0651606456; grid_points[1267].w_total = 0.0288023961; grid_points[1267].i_atom = 1;
    grid_points[1268].x =   0.5242619359; grid_points[1268].y =  -0.5242619359; grid_points[1268].z =   3.0971272145; grid_points[1268].w_fixed =    0.0761406613; grid_points[1268].w_total = 0.0448544141; grid_points[1268].i_atom = 1;
    grid_points[1269].x =   0.6338067859; grid_points[1269].y =   0.0000000000; grid_points[1269].z =   3.0316902052; grid_points[1269].w_fixed =    0.0769303986; grid_points[1269].w_total = 0.0412765861; grid_points[1269].i_atom = 1;
    grid_points[1270].x =   1.2035918363; grid_points[1270].y =  -0.3764723484; grid_points[1270].z =   2.9915998565; grid_points[1270].w_fixed =    0.1838706234; grid_points[1270].w_total = 0.0694127716; grid_points[1270].i_atom = 1;
    grid_points[1271].x =   0.3764723484; grid_points[1271].y =  -1.2035918363; grid_points[1271].z =   2.9915998565; grid_points[1271].w_fixed =    0.1838706234; grid_points[1271].w_total = 0.0694127716; grid_points[1271].i_atom = 1;
    grid_points[1272].x =   0.6897944676; grid_points[1272].y =  -0.6897944676; grid_points[1272].z =   2.7504200046; grid_points[1272].w_fixed =    0.1774472638; grid_points[1272].w_total = 0.0390890156; grid_points[1272].i_atom = 1;
    grid_points[1273].x =   2.0299139882; grid_points[1273].y =   0.0000000000; grid_points[1273].z =   3.0894146273; grid_points[1273].w_fixed =    0.4319923443; grid_points[1273].w_total = 0.1180180438; grid_points[1273].i_atom = 1;
    grid_points[1274].x =   2.5684830011; grid_points[1274].y =  -1.2263012757; grid_points[1274].z =   2.9688904170; grid_points[1274].w_fixed =    1.0767762850; grid_points[1274].w_total = 0.0592143299; grid_points[1274].i_atom = 1;
    grid_points[1275].x =   2.5220049911; grid_points[1275].y =   7.5660149734; grid_points[1275].z =   1.6731867016; grid_points[1275].w_fixed =   55.2891360098; grid_points[1275].w_total = 0.0000000000; grid_points[1275].i_atom = 1;
    grid_points[1276].x =   2.9827802928; grid_points[1276].y =   0.0000000000; grid_points[1276].z =   1.2124113999; grid_points[1276].w_fixed =    6.7598180997; grid_points[1276].w_total = 0.0000096089; grid_points[1276].i_atom = 1;
    grid_points[1277].x =  -8.3645442750; grid_points[1277].y =   0.0000000000; grid_points[1277].z =   4.1951916927; grid_points[1277].w_fixed =   34.8025870819; grid_points[1277].w_total = 0.0000000000; grid_points[1277].i_atom = 1;
    grid_points[1278].x =   2.5220049911; grid_points[1278].y =  -2.5220049911; grid_points[1278].z =  -3.3708232807; grid_points[1278].w_fixed =   55.2891360098; grid_points[1278].w_total = 0.0000000000; grid_points[1278].i_atom = 1;
    grid_points[1279].x =   0.2452657651; grid_points[1279].y =  -1.2787244278; grid_points[1279].z =   4.4404574578; grid_points[1279].w_fixed =    0.0651606456; grid_points[1279].w_total = 0.0648319178; grid_points[1279].i_atom = 1;
    grid_points[1280].x =   0.2452657651; grid_points[1280].y =  -1.2787244278; grid_points[1280].z =   3.9499259276; grid_points[1280].w_fixed =    0.0651606456; grid_points[1280].w_total = 0.0634126219; grid_points[1280].i_atom = 1;
    grid_points[1281].x =   0.3227069453; grid_points[1281].y =  -1.6824739228; grid_points[1281].z =   4.5178986380; grid_points[1281].w_fixed =    0.1518581277; grid_points[1281].w_total = 0.1485902764; grid_points[1281].i_atom = 1;
    grid_points[1282].x =   0.3227069453; grid_points[1282].y =  -1.6824739228; grid_points[1282].z =   3.8724847475; grid_points[1282].w_fixed =    0.1518581277; grid_points[1282].w_total = 0.1376185851; grid_points[1282].i_atom = 1;
    grid_points[1283].x =   0.6897944676; grid_points[1283].y =  -1.4447716881; grid_points[1283].z =   4.8849861603; grid_points[1283].w_fixed =    0.1774472638; grid_points[1283].w_total = 0.1658131153; grid_points[1283].i_atom = 1;
    grid_points[1284].x =   0.6897944676; grid_points[1284].y =  -1.4447716881; grid_points[1284].z =   3.5053972252; grid_points[1284].w_fixed =    0.1774472638; grid_points[1284].w_total = 0.1314016342; grid_points[1284].i_atom = 1;
    grid_points[1285].x =   0.8339274406; grid_points[1285].y =  -1.5308700368; grid_points[1285].z =   4.1951916927; grid_points[1285].w_fixed =    0.1792877615; grid_points[1285].w_total = 0.1738017119; grid_points[1285].i_atom = 1;
    grid_points[1286].x =   1.3345778006; grid_points[1286].y =  -1.3345778006; grid_points[1286].z =   5.5297694934; grid_points[1286].w_fixed =    0.4363920800; grid_points[1286].w_total = 0.2072585535; grid_points[1286].i_atom = 1;
    grid_points[1287].x =   0.4279052608; grid_points[1287].y =  -2.2309387920; grid_points[1287].z =   4.6230969535; grid_points[1287].w_fixed =    0.3659008737; grid_points[1287].w_total = 0.3324753846; grid_points[1287].i_atom = 1;
    grid_points[1288].x =   0.4279052608; grid_points[1288].y =  -2.2309387920; grid_points[1288].z =   3.7672864320; grid_points[1288].w_fixed =    0.3659008737; grid_points[1288].w_total = 0.2578929044; grid_points[1288].i_atom = 1;
    grid_points[1289].x =   1.5959473017; grid_points[1289].y =  -1.5959473017; grid_points[1289].z =   4.6943891868; grid_points[1289].w_fixed =    0.4430347113; grid_points[1289].w_total = 0.3987026542; grid_points[1289].i_atom = 1;
    grid_points[1290].x =   1.5959473017; grid_points[1290].y =  -1.5959473017; grid_points[1290].z =   3.6959941987; grid_points[1290].w_fixed =    0.4430347113; grid_points[1290].w_total = 0.2967672396; grid_points[1290].i_atom = 1;
    grid_points[1291].x =   0.4991974941; grid_points[1291].y =  -1.5959473017; grid_points[1291].z =   5.7911389945; grid_points[1291].w_fixed =    0.4430347113; grid_points[1291].w_total = 0.1354214943; grid_points[1291].i_atom = 1;
    grid_points[1292].x =   0.9146585961; grid_points[1292].y =  -1.9157486847; grid_points[1292].z =   5.1098502889; grid_points[1292].w_fixed =    0.4275576805; grid_points[1292].w_total = 0.3174678581; grid_points[1292].i_atom = 1;
    grid_points[1293].x =   0.9146585961; grid_points[1293].y =  -1.9157486847; grid_points[1293].z =   3.2805330966; grid_points[1293].w_fixed =    0.4275576805; grid_points[1293].w_total = 0.1674035783; grid_points[1293].i_atom = 1;
    grid_points[1294].x =   1.1057770654; grid_points[1294].y =  -2.0299139882; grid_points[1294].z =   4.1951916927; grid_points[1294].w_fixed =    0.4319923443; grid_points[1294].w_total = 0.3838311380; grid_points[1294].i_atom = 1;
    grid_points[1295].x =   1.7892954446; grid_points[1295].y =  -1.7892954446; grid_points[1295].z =   5.9844871374; grid_points[1295].w_fixed =    1.0990251469; grid_points[1295].w_total = 0.0671861197; grid_points[1295].i_atom = 1;
    grid_points[1296].x =   0.5737012360; grid_points[1296].y =  -2.9910647516; grid_points[1296].z =   4.7688929288; grid_points[1296].w_fixed =    0.9214976165; grid_points[1296].w_total = 0.4792995334; grid_points[1296].i_atom = 1;
    grid_points[1297].x =   0.5737012360; grid_points[1297].y =  -2.9910647516; grid_points[1297].z =   3.6214904567; grid_points[1297].w_fixed =    0.9214976165; grid_points[1297].w_total = 0.1977834187; grid_points[1297].i_atom = 1;
    grid_points[1298].x =   2.1397188201; grid_points[1298].y =  -2.1397188201; grid_points[1298].z =   4.8644758676; grid_points[1298].w_fixed =    1.1157541830; grid_points[1298].w_total = 0.5622597816; grid_points[1298].i_atom = 1;
    grid_points[1299].x =   2.1397188201; grid_points[1299].y =  -2.1397188201; grid_points[1299].z =   3.5259075179; grid_points[1299].w_fixed =    1.1157541830; grid_points[1299].w_total = 0.2172616206; grid_points[1299].i_atom = 1;
    grid_points[1300].x =   1.2263012757; grid_points[1300].y =  -2.5684830011; grid_points[1300].z =   5.4214929684; grid_points[1300].w_fixed =    1.0767762850; grid_points[1300].w_total = 0.2278386018; grid_points[1300].i_atom = 1;
    grid_points[1301].x =   1.4825376722; grid_points[1301].y =  -2.7215467321; grid_points[1301].z =   4.1951916927; grid_points[1301].w_fixed =    1.0879446980; grid_points[1301].w_total = 0.5523500480; grid_points[1301].i_atom = 1;
    grid_points[1302].x =   2.7215467321; grid_points[1302].y =  -1.4825376722; grid_points[1302].z =   4.1951916927; grid_points[1302].w_fixed =    1.0879446980; grid_points[1302].w_total = 0.5523500480; grid_points[1302].i_atom = 1;
    grid_points[1303].x =   0.0009212695; grid_points[1303].y =   0.0000000000; grid_points[1303].z =   4.1951916927; grid_points[1303].w_fixed =    0.0000000034; grid_points[1303].w_total = 0.0000000034; grid_points[1303].i_atom = 1;
    grid_points[1304].x =   0.0039529982; grid_points[1304].y =   0.0000000000; grid_points[1304].z =   4.1951916927; grid_points[1304].w_fixed =    0.0000001386; grid_points[1304].w_total = 0.0000001386; grid_points[1304].i_atom = 1;
    grid_points[1305].x =   0.0095652797; grid_points[1305].y =   0.0000000000; grid_points[1305].z =   4.1951916927; grid_points[1305].w_fixed =    0.0000013577; grid_points[1305].w_total = 0.0000013577; grid_points[1305].i_atom = 1;
    grid_points[1306].x =   0.0183381693; grid_points[1306].y =   0.0000000000; grid_points[1306].z =   4.1951916927; grid_points[1306].w_fixed =    0.0000074515; grid_points[1306].w_total = 0.0000074515; grid_points[1306].i_atom = 1;
    grid_points[1307].x =   0.0309915062; grid_points[1307].y =   0.0000000000; grid_points[1307].z =   4.1951916927; grid_points[1307].w_fixed =    0.0000299246; grid_points[1307].w_total = 0.0000299246; grid_points[1307].i_atom = 1;
    grid_points[1308].x =   0.0484242284; grid_points[1308].y =   0.0000000000; grid_points[1308].z =   4.1951916927; grid_points[1308].w_fixed =    0.0000990913; grid_points[1308].w_total = 0.0000990913; grid_points[1308].i_atom = 1;
    grid_points[1309].x =   0.0717667204; grid_points[1309].y =   0.0000000000; grid_points[1309].z =   4.1951916927; grid_points[1309].w_fixed =    0.0002885050; grid_points[1309].w_total = 0.0002885050; grid_points[1309].i_atom = 1;
    grid_points[1310].x =   0.1024512602; grid_points[1310].y =   0.0000000000; grid_points[1310].z =   4.1951916927; grid_points[1310].w_fixed =    0.0007678018; grid_points[1310].w_total = 0.0007678018; grid_points[1310].i_atom = 1;
    grid_points[1311].x =   0.1423079366; grid_points[1311].y =   0.0000000000; grid_points[1311].z =   4.1951916927; grid_points[1311].w_fixed =    0.0005474789; grid_points[1311].w_total = 0.0005474789; grid_points[1311].i_atom = 1;
    grid_points[1312].x =   0.1006269070; grid_points[1312].y =   0.0000000000; grid_points[1312].z =   4.2958185997; grid_points[1312].w_fixed =    0.0004379831; grid_points[1312].w_total = 0.0004379831; grid_points[1312].i_atom = 1;
    grid_points[1313].x =   0.1006269070; grid_points[1313].y =   0.0000000000; grid_points[1313].z =   4.0945647857; grid_points[1313].w_fixed =    0.0004379831; grid_points[1313].w_total = 0.0004379831; grid_points[1313].i_atom = 1;
    grid_points[1314].x =   0.1006269070; grid_points[1314].y =  -0.1006269070; grid_points[1314].z =   4.1951916927; grid_points[1314].w_fixed =    0.0004379831; grid_points[1314].w_total = 0.0004379831; grid_points[1314].i_atom = 1;
    grid_points[1315].x =   0.0821615255; grid_points[1315].y =  -0.0821615255; grid_points[1315].z =   4.2773532183; grid_points[1315].w_fixed =    0.0003695483; grid_points[1315].w_total = 0.0003695483; grid_points[1315].i_atom = 1;
    grid_points[1316].x =   0.0821615255; grid_points[1316].y =  -0.0821615255; grid_points[1316].z =   4.1130301672; grid_points[1316].w_fixed =    0.0003695483; grid_points[1316].w_total = 0.0003695482; grid_points[1316].i_atom = 1;
    grid_points[1317].x =   0.1936969137; grid_points[1317].y =   0.0000000000; grid_points[1317].z =   4.1951916927; grid_points[1317].w_fixed =    0.0013046077; grid_points[1317].w_total = 0.0013046076; grid_points[1317].i_atom = 1;
    grid_points[1318].x =   0.1369644012; grid_points[1318].y =   0.0000000000; grid_points[1318].z =   4.3321560939; grid_points[1318].w_fixed =    0.0010436861; grid_points[1318].w_total = 0.0010436860; grid_points[1318].i_atom = 1;
    grid_points[1319].x =   0.1369644012; grid_points[1319].y =   0.0000000000; grid_points[1319].z =   4.0582272915; grid_points[1319].w_fixed =    0.0010436861; grid_points[1319].w_total = 0.0010436847; grid_points[1319].i_atom = 1;
    grid_points[1320].x =   0.1369644012; grid_points[1320].y =  -0.1369644012; grid_points[1320].z =   4.1951916927; grid_points[1320].w_fixed =    0.0010436861; grid_points[1320].w_total = 0.0010436861; grid_points[1320].i_atom = 1;
    grid_points[1321].x =   0.1118309653; grid_points[1321].y =  -0.1118309653; grid_points[1321].z =   4.3070226580; grid_points[1321].w_fixed =    0.0008806102; grid_points[1321].w_total = 0.0008806101; grid_points[1321].i_atom = 1;
    grid_points[1322].x =   0.1118309653; grid_points[1322].y =  -0.1118309653; grid_points[1322].z =   4.0833607275; grid_points[1322].w_fixed =    0.0008806102; grid_points[1322].w_total = 0.0008806095; grid_points[1322].i_atom = 1;
    grid_points[1323].x =   0.2596933691; grid_points[1323].y =   0.0000000000; grid_points[1323].z =   4.1951916927; grid_points[1323].w_fixed =    0.0030086947; grid_points[1323].w_total = 0.0030086941; grid_points[1323].i_atom = 1;
    grid_points[1324].x =   0.1836309423; grid_points[1324].y =   0.0000000000; grid_points[1324].z =   4.3788226351; grid_points[1324].w_fixed =    0.0024069557; grid_points[1324].w_total = 0.0024069533; grid_points[1324].i_atom = 1;
    grid_points[1325].x =   0.1836309423; grid_points[1325].y =   0.0000000000; grid_points[1325].z =   4.0115607504; grid_points[1325].w_fixed =    0.0024069557; grid_points[1325].w_total = 0.0024069249; grid_points[1325].i_atom = 1;
    grid_points[1326].x =   0.1836309423; grid_points[1326].y =  -0.1836309423; grid_points[1326].z =   4.1951916927; grid_points[1326].w_fixed =    0.0024069557; grid_points[1326].w_total = 0.0024069553; grid_points[1326].i_atom = 1;
    grid_points[1327].x =   0.1499340366; grid_points[1327].y =  -0.1499340366; grid_points[1327].z =   4.3451257293; grid_points[1327].w_fixed =    0.0020308689; grid_points[1327].w_total = 0.0020308678; grid_points[1327].i_atom = 1;
    grid_points[1328].x =   0.1499340366; grid_points[1328].y =  -0.1499340366; grid_points[1328].z =   4.0452576562; grid_points[1328].w_fixed =    0.0020308689; grid_points[1328].w_total = 0.0020308551; grid_points[1328].i_atom = 1;
    grid_points[1329].x =   0.3443500689; grid_points[1329].y =   0.0000000000; grid_points[1329].z =   4.1951916927; grid_points[1329].w_fixed =    0.0067871679; grid_points[1329].w_total = 0.0067871580; grid_points[1329].i_atom = 1;
    grid_points[1330].x =   0.2434922688; grid_points[1330].y =   0.0000000000; grid_points[1330].z =   4.4386839615; grid_points[1330].w_fixed =    0.0054297343; grid_points[1330].w_total = 0.0054296870; grid_points[1330].i_atom = 1;
    grid_points[1331].x =   0.2434922688; grid_points[1331].y =   0.0000000000; grid_points[1331].z =   3.9516994239; grid_points[1331].w_fixed =    0.0054297343; grid_points[1331].w_total = 0.0054291838; grid_points[1331].i_atom = 1;
    grid_points[1332].x =   0.2434922688; grid_points[1332].y =  -0.2434922688; grid_points[1332].z =   4.1951916927; grid_points[1332].w_fixed =    0.0054297343; grid_points[1332].w_total = 0.0054297264; grid_points[1332].i_atom = 1;
    grid_points[1333].x =   0.1988106050; grid_points[1333].y =  -0.1988106050; grid_points[1333].z =   4.3940022977; grid_points[1333].w_fixed =    0.0045813383; grid_points[1333].w_total = 0.0045813179; grid_points[1333].i_atom = 1;
    grid_points[1334].x =   0.1988106050; grid_points[1334].y =  -0.1988106050; grid_points[1334].z =   3.9963810878; grid_points[1334].w_fixed =    0.0045813383; grid_points[1334].w_total = 0.0045810919; grid_points[1334].i_atom = 1;
    grid_points[1335].x =   0.4530765180; grid_points[1335].y =   0.0000000000; grid_points[1335].z =   4.1951916927; grid_points[1335].w_fixed =    0.0040293310; grid_points[1335].w_total = 0.0040292893; grid_points[1335].i_atom = 1;
    grid_points[1336].x =   0.3203734782; grid_points[1336].y =   0.0000000000; grid_points[1336].z =   4.5155651710; grid_points[1336].w_fixed =    0.0071632550; grid_points[1336].w_total = 0.0071627669; grid_points[1336].i_atom = 1;
    grid_points[1337].x =   0.3203734782; grid_points[1337].y =   0.0000000000; grid_points[1337].z =   3.8748182145; grid_points[1337].w_fixed =    0.0071632550; grid_points[1337].w_total = 0.0071582150; grid_points[1337].i_atom = 1;
    grid_points[1338].x =   0.3203734782; grid_points[1338].y =  -0.3203734782; grid_points[1338].z =   4.1951916927; grid_points[1338].w_fixed =    0.0071632550; grid_points[1338].w_total = 0.0071631809; grid_points[1338].i_atom = 1;
    grid_points[1339].x =   0.2615838496; grid_points[1339].y =  -0.2615838496; grid_points[1339].z =   4.4567755423; grid_points[1339].w_fixed =    0.0066932539; grid_points[1339].w_total = 0.0066930222; grid_points[1339].i_atom = 1;
    grid_points[1340].x =   0.2615838496; grid_points[1340].y =  -0.2615838496; grid_points[1340].z =   3.9336078431; grid_points[1340].w_fixed =    0.0066932539; grid_points[1340].w_total = 0.0066907381; grid_points[1340].i_atom = 1;
    grid_points[1341].x =   0.1366077101; grid_points[1341].y =  -0.1366077101; grid_points[1341].z =   4.6050148231; grid_points[1341].w_fixed =    0.0064011973; grid_points[1341].w_total = 0.0064000836; grid_points[1341].i_atom = 1;
    grid_points[1342].x =   0.1366077101; grid_points[1342].y =  -0.1366077101; grid_points[1342].z =   3.7853685624; grid_points[1342].w_fixed =    0.0064011973; grid_points[1342].w_total = 0.0063905489; grid_points[1342].i_atom = 1;
    grid_points[1343].x =   0.1366077101; grid_points[1343].y =  -0.4098231304; grid_points[1343].z =   4.3317994029; grid_points[1343].w_fixed =    0.0064011973; grid_points[1343].w_total = 0.0064011513; grid_points[1343].i_atom = 1;
    grid_points[1344].x =   0.1366077101; grid_points[1344].y =  -0.4098231304; grid_points[1344].z =   4.0585839826; grid_points[1344].w_fixed =    0.0064011973; grid_points[1344].w_total = 0.0064006801; grid_points[1344].i_atom = 1;
    grid_points[1345].x =   0.4098231304; grid_points[1345].y =  -0.1366077101; grid_points[1345].z =   4.3317994029; grid_points[1345].w_fixed =    0.0064011973; grid_points[1345].w_total = 0.0064011513; grid_points[1345].i_atom = 1;
    grid_points[1346].x =   0.4098231304; grid_points[1346].y =  -0.1366077101; grid_points[1346].z =   4.0585839826; grid_points[1346].w_fixed =    0.0064011973; grid_points[1346].w_total = 0.0064006801; grid_points[1346].i_atom = 1;
    grid_points[1347].x =   0.5931967983; grid_points[1347].y =   0.0000000000; grid_points[1347].z =   4.1951916927; grid_points[1347].w_fixed =    0.0089219324; grid_points[1347].w_total = 0.0089213452; grid_points[1347].i_atom = 1;
    grid_points[1348].x =   0.4194534787; grid_points[1348].y =   0.0000000000; grid_points[1348].z =   4.6146451714; grid_points[1348].w_fixed =    0.0158612132; grid_points[1348].w_total = 0.0158535719; grid_points[1348].i_atom = 1;
    grid_points[1349].x =   0.4194534787; grid_points[1349].y =   0.0000000000; grid_points[1349].z =   3.7757382141; grid_points[1349].w_fixed =    0.0158612132; grid_points[1349].w_total = 0.0157937824; grid_points[1349].i_atom = 1;
    grid_points[1350].x =   0.4194534787; grid_points[1350].y =  -0.4194534787; grid_points[1350].z =   4.1951916927; grid_points[1350].w_fixed =    0.0158612132; grid_points[1350].w_total = 0.0158601692; grid_points[1350].i_atom = 1;
    grid_points[1351].x =   0.3424823312; grid_points[1351].y =  -0.3424823312; grid_points[1351].z =   4.5376740239; grid_points[1351].w_fixed =    0.0148205147; grid_points[1351].w_total = 0.0148168993; grid_points[1351].i_atom = 1;
    grid_points[1352].x =   0.3424823312; grid_points[1352].y =  -0.3424823312; grid_points[1352].z =   3.8527093615; grid_points[1352].w_fixed =    0.0148205147; grid_points[1352].w_total = 0.0147863672; grid_points[1352].i_atom = 1;
    grid_points[1353].x =   0.1788555643; grid_points[1353].y =  -0.1788555643; grid_points[1353].z =   4.7317583855; grid_points[1353].w_fixed =    0.0141738295; grid_points[1353].w_total = 0.0141564043; grid_points[1353].i_atom = 1;
    grid_points[1354].x =   0.1788555643; grid_points[1354].y =  -0.1788555643; grid_points[1354].z =   3.6586249999; grid_points[1354].w_fixed =    0.0141738295; grid_points[1354].w_total = 0.0140351874; grid_points[1354].i_atom = 1;
    grid_points[1355].x =   0.1788555643; grid_points[1355].y =  -0.5365666928; grid_points[1355].z =   4.3740472570; grid_points[1355].w_fixed =    0.0141738295; grid_points[1355].w_total = 0.0141731261; grid_points[1355].i_atom = 1;
    grid_points[1356].x =   0.1788555643; grid_points[1356].y =  -0.5365666928; grid_points[1356].z =   4.0163361285; grid_points[1356].w_fixed =    0.0141738295; grid_points[1356].w_total = 0.0141666532; grid_points[1356].i_atom = 1;
    grid_points[1357].x =   0.5365666928; grid_points[1357].y =  -0.1788555643; grid_points[1357].z =   4.3740472570; grid_points[1357].w_fixed =    0.0141738295; grid_points[1357].w_total = 0.0141731261; grid_points[1357].i_atom = 1;
    grid_points[1358].x =   0.5365666928; grid_points[1358].y =  -0.1788555643; grid_points[1358].z =   4.0163361285; grid_points[1358].w_fixed =    0.0141738295; grid_points[1358].w_total = 0.0141666532; grid_points[1358].i_atom = 1;
    grid_points[1359].x =   0.7747876550; grid_points[1359].y =   0.0000000000; grid_points[1359].z =   4.1951916927; grid_points[1359].w_fixed =    0.0059666325; grid_points[1359].w_total = 0.0059644033; grid_points[1359].i_atom = 1;
    grid_points[1360].x =   0.4473238612; grid_points[1360].y =  -0.4473238612; grid_points[1360].z =   4.6425155539; grid_points[1360].w_fixed =    0.0152642382; grid_points[1360].w_total = 0.0152405851; grid_points[1360].i_atom = 1;
    grid_points[1361].x =   0.4473238612; grid_points[1361].y =  -0.4473238612; grid_points[1361].z =   3.7478678316; grid_points[1361].w_fixed =    0.0152642382; grid_points[1361].w_total = 0.0150785295; grid_points[1361].i_atom = 1;
    grid_points[1362].x =   0.1434253090; grid_points[1362].y =  -0.1434253090; grid_points[1362].z =   4.9429578806; grid_points[1362].w_fixed =    0.0127985780; grid_points[1362].w_total = 0.0126709482; grid_points[1362].i_atom = 1;
    grid_points[1363].x =   0.1434253090; grid_points[1363].y =  -0.1434253090; grid_points[1363].z =   3.4474255048; grid_points[1363].w_fixed =    0.0127985780; grid_points[1363].w_total = 0.0120391551; grid_points[1363].i_atom = 1;
    grid_points[1364].x =   0.1434253090; grid_points[1364].y =  -0.7477661879; grid_points[1364].z =   4.3386170017; grid_points[1364].w_fixed =    0.0127985780; grid_points[1364].w_total = 0.0127960888; grid_points[1364].i_atom = 1;
    grid_points[1365].x =   0.1434253090; grid_points[1365].y =  -0.7477661879; grid_points[1365].z =   4.0517663837; grid_points[1365].w_fixed =    0.0127985780; grid_points[1365].w_total = 0.0127814266; grid_points[1365].i_atom = 1;
    grid_points[1366].x =   0.7477661879; grid_points[1366].y =  -0.1434253090; grid_points[1366].z =   4.3386170017; grid_points[1366].w_fixed =    0.0127985780; grid_points[1366].w_total = 0.0127960888; grid_points[1366].i_atom = 1;
    grid_points[1367].x =   0.7477661879; grid_points[1367].y =  -0.1434253090; grid_points[1367].z =   4.0517663837; grid_points[1367].w_fixed =    0.0127985780; grid_points[1367].w_total = 0.0127814266; grid_points[1367].i_atom = 1;
    grid_points[1368].x =   0.5349297050; grid_points[1368].y =  -0.5349297050; grid_points[1368].z =   4.3625127365; grid_points[1368].w_fixed =    0.0154965859; grid_points[1368].w_total = 0.0154933439; grid_points[1368].i_atom = 1;
    grid_points[1369].x =   0.5349297050; grid_points[1369].y =  -0.5349297050; grid_points[1369].z =   4.0278706490; grid_points[1369].w_fixed =    0.0154965859; grid_points[1369].w_total = 0.0154711868; grid_points[1369].i_atom = 1;
    grid_points[1370].x =   0.5349297050; grid_points[1370].y =  -0.1673210437; grid_points[1370].z =   4.7301213978; grid_points[1370].w_fixed =    0.0154965859; grid_points[1370].w_total = 0.0154531108; grid_points[1370].i_atom = 1;
    grid_points[1371].x =   0.5349297050; grid_points[1371].y =  -0.1673210437; grid_points[1371].z =   3.6602619877; grid_points[1371].w_fixed =    0.0154965859; grid_points[1371].w_total = 0.0151815178; grid_points[1371].i_atom = 1;
    grid_points[1372].x =   0.1673210437; grid_points[1372].y =  -0.5349297050; grid_points[1372].z =   4.7301213978; grid_points[1372].w_fixed =    0.0154965859; grid_points[1372].w_total = 0.0154531108; grid_points[1372].i_atom = 1;
    grid_points[1373].x =   0.1673210437; grid_points[1373].y =  -0.5349297050; grid_points[1373].z =   3.6602619877; grid_points[1373].w_fixed =    0.0154965859; grid_points[1373].w_total = 0.0151815178; grid_points[1373].i_atom = 1;
    grid_points[1374].x =   0.3065753189; grid_points[1374].y =  -0.3065753189; grid_points[1374].z =   4.8373124430; grid_points[1374].w_fixed =    0.0149552262; grid_points[1374].w_total = 0.0148734574; grid_points[1374].i_atom = 1;
    grid_points[1375].x =   0.3065753189; grid_points[1375].y =  -0.3065753189; grid_points[1375].z =   3.5530709425; grid_points[1375].w_fixed =    0.0149552262; grid_points[1375].w_total = 0.0144183652; grid_points[1375].i_atom = 1;
    grid_points[1376].x =   0.3065753189; grid_points[1376].y =  -0.6421207503; grid_points[1376].z =   4.5017670117; grid_points[1376].w_fixed =    0.0149552262; grid_points[1376].w_total = 0.0149471131; grid_points[1376].i_atom = 1;
    grid_points[1377].x =   0.3065753189; grid_points[1377].y =  -0.6421207503; grid_points[1377].z =   3.8886163738; grid_points[1377].w_fixed =    0.0149552262; grid_points[1377].w_total = 0.0148835083; grid_points[1377].i_atom = 1;
    grid_points[1378].x =   0.6421207503; grid_points[1378].y =  -0.3065753189; grid_points[1378].z =   4.5017670117; grid_points[1378].w_fixed =    0.0149552262; grid_points[1378].w_total = 0.0149471131; grid_points[1378].i_atom = 1;
    grid_points[1379].x =   0.6421207503; grid_points[1379].y =  -0.3065753189; grid_points[1379].z =   3.8886163738; grid_points[1379].w_fixed =    0.0149552262; grid_points[1379].w_total = 0.0148835083; grid_points[1379].i_atom = 1;
    grid_points[1380].x =   0.3706344180; grid_points[1380].y =  -0.6803866830; grid_points[1380].z =   4.1951916927; grid_points[1380].w_fixed =    0.0151103430; grid_points[1380].w_total = 0.0151046974; grid_points[1380].i_atom = 1;
    grid_points[1381].x =   0.6803866830; grid_points[1381].y =  -0.3706344180; grid_points[1381].z =   4.1951916927; grid_points[1381].w_fixed =    0.0151103430; grid_points[1381].w_total = 0.0151046974; grid_points[1381].i_atom = 1;
    grid_points[1382].x =   0.3706344180; grid_points[1382].y =   0.0000000000; grid_points[1382].z =   4.8755783758; grid_points[1382].w_fixed =    0.0151103430; grid_points[1382].w_total = 0.0150070014; grid_points[1382].i_atom = 1;
    grid_points[1383].x =   0.3706344180; grid_points[1383].y =   0.0000000000; grid_points[1383].z =   3.5148050097; grid_points[1383].w_fixed =    0.0151103430; grid_points[1383].w_total = 0.0144554995; grid_points[1383].i_atom = 1;
    grid_points[1384].x =   0.6803866830; grid_points[1384].y =   0.0000000000; grid_points[1384].z =   4.5658261108; grid_points[1384].w_fixed =    0.0151103430; grid_points[1384].w_total = 0.0150969464; grid_points[1384].i_atom = 1;
    grid_points[1385].x =   0.6803866830; grid_points[1385].y =   0.0000000000; grid_points[1385].z =   3.8245572747; grid_points[1385].w_fixed =    0.0151103430; grid_points[1385].w_total = 0.0149977013; grid_points[1385].i_atom = 1;
    grid_points[1386].x =   1.0119675493; grid_points[1386].y =   0.0000000000; grid_points[1386].z =   4.1951916927; grid_points[1386].w_fixed =    0.0133541260; grid_points[1386].w_total = 0.0133291503; grid_points[1386].i_atom = 1;
    grid_points[1387].x =   0.5842597370; grid_points[1387].y =  -0.5842597370; grid_points[1387].z =   4.7794514298; grid_points[1387].w_fixed =    0.0341634180; grid_points[1387].w_total = 0.0338660424; grid_points[1387].i_atom = 1;
    grid_points[1388].x =   0.5842597370; grid_points[1388].y =  -0.5842597370; grid_points[1388].z =   3.6109319557; grid_points[1388].w_fixed =    0.0341634180; grid_points[1388].w_total = 0.0323366953; grid_points[1388].i_atom = 1;
    grid_points[1389].x =   0.1873310158; grid_points[1389].y =  -0.1873310158; grid_points[1389].z =   5.1718658973; grid_points[1389].w_fixed =    0.0286449390; grid_points[1389].w_total = 0.0271394907; grid_points[1389].i_atom = 1;
    grid_points[1390].x =   0.1873310158; grid_points[1390].y =  -0.1873310158; grid_points[1390].z =   3.2185174881; grid_points[1390].w_fixed =    0.0286449390; grid_points[1390].w_total = 0.0223716243; grid_points[1390].i_atom = 1;
    grid_points[1391].x =   0.1873310158; grid_points[1391].y =  -0.9766742046; grid_points[1391].z =   4.3825227086; grid_points[1391].w_fixed =    0.0286449390; grid_points[1391].w_total = 0.0286147757; grid_points[1391].i_atom = 1;
    grid_points[1392].x =   0.1873310158; grid_points[1392].y =  -0.9766742046; grid_points[1392].z =   4.0078606769; grid_points[1392].w_fixed =    0.0286449390; grid_points[1392].w_total = 0.0284587127; grid_points[1392].i_atom = 1;
    grid_points[1393].x =   0.9766742046; grid_points[1393].y =  -0.1873310158; grid_points[1393].z =   4.3825227086; grid_points[1393].w_fixed =    0.0286449390; grid_points[1393].w_total = 0.0286147757; grid_points[1393].i_atom = 1;
    grid_points[1394].x =   0.9766742046; grid_points[1394].y =  -0.1873310158; grid_points[1394].z =   4.0078606769; grid_points[1394].w_fixed =    0.0286449390; grid_points[1394].w_total = 0.0284587127; grid_points[1394].i_atom = 1;
    grid_points[1395].x =   0.6986836963; grid_points[1395].y =  -0.6986836963; grid_points[1395].z =   4.4137334641; grid_points[1395].w_fixed =    0.0346834435; grid_points[1395].w_total = 0.0346436110; grid_points[1395].i_atom = 1;
    grid_points[1396].x =   0.6986836963; grid_points[1396].y =  -0.6986836963; grid_points[1396].z =   3.9766499214; grid_points[1396].w_fixed =    0.0346834435; grid_points[1396].w_total = 0.0344092150; grid_points[1396].i_atom = 1;
    grid_points[1397].x =   0.6986836963; grid_points[1397].y =  -0.2185417714; grid_points[1397].z =   4.8938753891; grid_points[1397].w_fixed =    0.0346834435; grid_points[1397].w_total = 0.0341432423; grid_points[1397].i_atom = 1;
    grid_points[1398].x =   0.6986836963; grid_points[1398].y =  -0.2185417714; grid_points[1398].z =   3.4965079964; grid_points[1398].w_fixed =    0.0346834435; grid_points[1398].w_total = 0.0317117114; grid_points[1398].i_atom = 1;
    grid_points[1399].x =   0.2185417714; grid_points[1399].y =  -0.6986836963; grid_points[1399].z =   4.8938753891; grid_points[1399].w_fixed =    0.0346834435; grid_points[1399].w_total = 0.0341432423; grid_points[1399].i_atom = 1;
    grid_points[1400].x =   0.2185417714; grid_points[1400].y =  -0.6986836963; grid_points[1400].z =   3.4965079964; grid_points[1400].w_fixed =    0.0346834435; grid_points[1400].w_total = 0.0317117114; grid_points[1400].i_atom = 1;
    grid_points[1401].x =   0.4004249063; grid_points[1401].y =  -0.4004249063; grid_points[1401].z =   5.0338800196; grid_points[1401].w_fixed =    0.0334718077; grid_points[1401].w_total = 0.0324778009; grid_points[1401].i_atom = 1;
    grid_points[1402].x =   0.4004249063; grid_points[1402].y =  -0.4004249063; grid_points[1402].z =   3.3565033658; grid_points[1402].w_fixed =    0.0334718077; grid_points[1402].w_total = 0.0287086336; grid_points[1402].i_atom = 1;
    grid_points[1403].x =   0.4004249063; grid_points[1403].y =  -0.8386883269; grid_points[1403].z =   4.5956165991; grid_points[1403].w_fixed =    0.0334718077; grid_points[1403].w_total = 0.0333691081; grid_points[1403].i_atom = 1;
    grid_points[1404].x =   0.4004249063; grid_points[1404].y =  -0.8386883269; grid_points[1404].z =   3.7947667864; grid_points[1404].w_fixed =    0.0334718077; grid_points[1404].w_total = 0.0327270168; grid_points[1404].i_atom = 1;
    grid_points[1405].x =   0.8386883269; grid_points[1405].y =  -0.4004249063; grid_points[1405].z =   4.5956165991; grid_points[1405].w_fixed =    0.0334718077; grid_points[1405].w_total = 0.0333691081; grid_points[1405].i_atom = 1;
    grid_points[1406].x =   0.8386883269; grid_points[1406].y =  -0.4004249063; grid_points[1406].z =   3.7947667864; grid_points[1406].w_fixed =    0.0334718077; grid_points[1406].w_total = 0.0327270168; grid_points[1406].i_atom = 1;
    grid_points[1407].x =   0.4840939338; grid_points[1407].y =  -0.8886683207; grid_points[1407].z =   4.1951916927; grid_points[1407].w_fixed =    0.0338189800; grid_points[1407].w_total = 0.0337557378; grid_points[1407].i_atom = 1;
    grid_points[1408].x =   0.8886683207; grid_points[1408].y =  -0.4840939338; grid_points[1408].z =   4.1951916927; grid_points[1408].w_fixed =    0.0338189800; grid_points[1408].w_total = 0.0337557378; grid_points[1408].i_atom = 1;
    grid_points[1409].x =   0.4840939338; grid_points[1409].y =   0.0000000000; grid_points[1409].z =   5.0838600134; grid_points[1409].w_fixed =    0.0338189800; grid_points[1409].w_total = 0.0325750639; grid_points[1409].i_atom = 1;
    grid_points[1410].x =   0.4840939338; grid_points[1410].y =   0.0000000000; grid_points[1410].z =   3.3065233720; grid_points[1410].w_fixed =    0.0338189800; grid_points[1410].w_total = 0.0281497735; grid_points[1410].i_atom = 1;
    grid_points[1411].x =   0.8886683207; grid_points[1411].y =   0.0000000000; grid_points[1411].z =   4.6792856265; grid_points[1411].w_fixed =    0.0338189800; grid_points[1411].w_total = 0.0336495489; grid_points[1411].i_atom = 1;
    grid_points[1412].x =   0.8886683207; grid_points[1412].y =   0.0000000000; grid_points[1412].z =   3.7110977590; grid_points[1412].w_fixed =    0.0338189800; grid_points[1412].w_total = 0.0326753125; grid_points[1412].i_atom = 1;
    grid_points[1413].x =   1.3249327354; grid_points[1413].y =   0.0000000000; grid_points[1413].z =   4.1951916927; grid_points[1413].w_fixed =    0.0303775643; grid_points[1413].w_total = 0.0301307226; grid_points[1413].i_atom = 1;
    grid_points[1414].x =   0.7649502715; grid_points[1414].y =  -0.7649502715; grid_points[1414].z =   4.9601419642; grid_points[1414].w_fixed =    0.0777139157; grid_points[1414].w_total = 0.0744471256; grid_points[1414].i_atom = 1;
    grid_points[1415].x =   0.7649502715; grid_points[1415].y =  -0.7649502715; grid_points[1415].z =   3.4302414213; grid_points[1415].w_fixed =    0.0777139157; grid_points[1415].w_total = 0.0632755179; grid_points[1415].i_atom = 1;
    grid_points[1416].x =   0.2452657651; grid_points[1416].y =  -0.2452657651; grid_points[1416].z =   5.4739161206; grid_points[1416].w_fixed =    0.0651606456; grid_points[1416].w_total = 0.0509705226; grid_points[1416].i_atom = 1;
    grid_points[1417].x =   1.2787244278; grid_points[1417].y =  -0.2452657651; grid_points[1417].z =   4.4404574578; grid_points[1417].w_fixed =    0.0651606456; grid_points[1417].w_total = 0.0648319178; grid_points[1417].i_atom = 1;
    grid_points[1418].x =   1.2787244278; grid_points[1418].y =  -0.2452657651; grid_points[1418].z =   3.9499259276; grid_points[1418].w_fixed =    0.0651606456; grid_points[1418].w_total = 0.0634126219; grid_points[1418].i_atom = 1;
    grid_points[1419].x =   0.9147614482; grid_points[1419].y =  -0.9147614482; grid_points[1419].z =   4.4813205782; grid_points[1419].w_fixed =    0.0788968540; grid_points[1419].w_total = 0.0784571533; grid_points[1419].i_atom = 1;
    grid_points[1420].x =   0.9147614482; grid_points[1420].y =  -0.9147614482; grid_points[1420].z =   3.9090628073; grid_points[1420].w_fixed =    0.0788968540; grid_points[1420].w_total = 0.0763502659; grid_points[1420].i_atom = 1;
    grid_points[1421].x =   0.9147614482; grid_points[1421].y =  -0.2861288854; grid_points[1421].z =   5.1099531410; grid_points[1421].w_fixed =    0.0788968540; grid_points[1421].w_total = 0.0731312984; grid_points[1421].i_atom = 1;
    grid_points[1422].x =   0.9147614482; grid_points[1422].y =  -0.2861288854; grid_points[1422].z =   3.2804302445; grid_points[1422].w_fixed =    0.0788968540; grid_points[1422].w_total = 0.0571063680; grid_points[1422].i_atom = 1;
    grid_points[1423].x =   0.2861288854; grid_points[1423].y =  -0.9147614482; grid_points[1423].z =   5.1099531410; grid_points[1423].w_fixed =    0.0788968540; grid_points[1423].w_total = 0.0731312984; grid_points[1423].i_atom = 1;
    grid_points[1424].x =   0.2861288854; grid_points[1424].y =  -0.9147614482; grid_points[1424].z =   3.2804302445; grid_points[1424].w_fixed =    0.0788968540; grid_points[1424].w_total = 0.0571063680; grid_points[1424].i_atom = 1;
    grid_points[1425].x =   0.5242619359; grid_points[1425].y =  -0.5242619359; grid_points[1425].z =   5.2932561710; grid_points[1425].w_fixed =    0.0761406613; grid_points[1425].w_total = 0.0660727372; grid_points[1425].i_atom = 1;
    grid_points[1426].x =   0.5242619359; grid_points[1426].y =  -1.0980644783; grid_points[1426].z =   4.7194536287; grid_points[1426].w_fixed =    0.0761406613; grid_points[1426].w_total = 0.0749832017; grid_points[1426].i_atom = 1;
    grid_points[1427].x =   0.5242619359; grid_points[1427].y =  -1.0980644783; grid_points[1427].z =   3.6709297568; grid_points[1427].w_fixed =    0.0761406613; grid_points[1427].w_total = 0.0696636336; grid_points[1427].i_atom = 1;
    grid_points[1428].x =   1.0980644783; grid_points[1428].y =  -0.5242619359; grid_points[1428].z =   4.7194536287; grid_points[1428].w_fixed =    0.0761406613; grid_points[1428].w_total = 0.0749832017; grid_points[1428].i_atom = 1;
    grid_points[1429].x =   1.0980644783; grid_points[1429].y =  -0.5242619359; grid_points[1429].z =   3.6709297568; grid_points[1429].w_fixed =    0.0761406613; grid_points[1429].w_total = 0.0696636336; grid_points[1429].i_atom = 1;
    grid_points[1430].x =   0.6338067859; grid_points[1430].y =  -1.1635014876; grid_points[1430].z =   4.1951916927; grid_points[1430].w_fixed =    0.0769303986; grid_points[1430].w_total = 0.0763061913; grid_points[1430].i_atom = 1;
    grid_points[1431].x =   1.1635014876; grid_points[1431].y =  -0.6338067859; grid_points[1431].z =   4.1951916927; grid_points[1431].w_fixed =    0.0769303986; grid_points[1431].w_total = 0.0763061913; grid_points[1431].i_atom = 1;
    grid_points[1432].x =   0.6338067859; grid_points[1432].y =   0.0000000000; grid_points[1432].z =   5.3586931803; grid_points[1432].w_fixed =    0.0769303986; grid_points[1432].w_total = 0.0646229177; grid_points[1432].i_atom = 1;
    grid_points[1433].x =   1.1635014876; grid_points[1433].y =   0.0000000000; grid_points[1433].z =   4.8289984786; grid_points[1433].w_fixed =    0.0769303986; grid_points[1433].w_total = 0.0750352153; grid_points[1433].i_atom = 1;
    grid_points[1434].x =   1.1635014876; grid_points[1434].y =   0.0000000000; grid_points[1434].z =   3.5613849069; grid_points[1434].w_fixed =    0.0769303986; grid_points[1434].w_total = 0.0673622975; grid_points[1434].i_atom = 1;
    grid_points[1435].x =   1.7432722237; grid_points[1435].y =   0.0000000000; grid_points[1435].z =   4.1951916927; grid_points[1435].w_fixed =    0.0707954930; grid_points[1435].w_total = 0.0686117539; grid_points[1435].i_atom = 1;
    grid_points[1436].x =   1.0064786876; grid_points[1436].y =  -1.0064786876; grid_points[1436].z =   5.2016703804; grid_points[1436].w_fixed =    0.1811137632; grid_points[1436].w_total = 0.1504983152; grid_points[1436].i_atom = 1;
    grid_points[1437].x =   1.0064786876; grid_points[1437].y =  -1.0064786876; grid_points[1437].z =   3.1887130051; grid_points[1437].w_fixed =    0.1811137632; grid_points[1437].w_total = 0.0948100401; grid_points[1437].i_atom = 1;
    grid_points[1438].x =   0.3227069453; grid_points[1438].y =  -0.3227069453; grid_points[1438].z =   5.8776656155; grid_points[1438].w_fixed =    0.1518581277; grid_points[1438].w_total = 0.0603864433; grid_points[1438].i_atom = 1;
    grid_points[1439].x =   1.6824739228; grid_points[1439].y =  -0.3227069453; grid_points[1439].z =   4.5178986380; grid_points[1439].w_fixed =    0.1518581277; grid_points[1439].w_total = 0.1485902764; grid_points[1439].i_atom = 1;
    grid_points[1440].x =   1.6824739228; grid_points[1440].y =  -0.3227069453; grid_points[1440].z =   3.8724847475; grid_points[1440].w_fixed =    0.1518581277; grid_points[1440].w_total = 0.1376185851; grid_points[1440].i_atom = 1;
    grid_points[1441].x =   1.2035918363; grid_points[1441].y =  -1.2035918363; grid_points[1441].z =   4.5716640411; grid_points[1441].w_fixed =    0.1838706234; grid_points[1441].w_total = 0.1794730010; grid_points[1441].i_atom = 1;
    grid_points[1442].x =   1.2035918363; grid_points[1442].y =  -1.2035918363; grid_points[1442].z =   3.8187193444; grid_points[1442].w_fixed =    0.1838706234; grid_points[1442].w_total = 0.1635857757; grid_points[1442].i_atom = 1;
    grid_points[1443].x =   1.2035918363; grid_points[1443].y =  -0.3764723484; grid_points[1443].z =   5.3987835290; grid_points[1443].w_fixed =    0.1838706234; grid_points[1443].w_total = 0.1337193515; grid_points[1443].i_atom = 1;
    grid_points[1444].x =   0.3764723484; grid_points[1444].y =  -1.2035918363; grid_points[1444].z =   5.3987835290; grid_points[1444].w_fixed =    0.1838706234; grid_points[1444].w_total = 0.1337193515; grid_points[1444].i_atom = 1;
    grid_points[1445].x =   0.6897944676; grid_points[1445].y =  -0.6897944676; grid_points[1445].z =   5.6399633809; grid_points[1445].w_fixed =    0.1774472638; grid_points[1445].w_total = 0.1007712847; grid_points[1445].i_atom = 1;
    grid_points[1446].x =   1.4447716881; grid_points[1446].y =  -0.6897944676; grid_points[1446].z =   4.8849861603; grid_points[1446].w_fixed =    0.1774472638; grid_points[1446].w_total = 0.1658131153; grid_points[1446].i_atom = 1;
    grid_points[1447].x =   1.4447716881; grid_points[1447].y =  -0.6897944676; grid_points[1447].z =   3.5053972252; grid_points[1447].w_fixed =    0.1774472638; grid_points[1447].w_total = 0.1314016342; grid_points[1447].i_atom = 1;
    grid_points[1448].x =   1.5308700368; grid_points[1448].y =  -0.8339274406; grid_points[1448].z =   4.1951916927; grid_points[1448].w_fixed =    0.1792877615; grid_points[1448].w_total = 0.1738017119; grid_points[1448].i_atom = 1;
    grid_points[1449].x =   0.8339274406; grid_points[1449].y =   0.0000000000; grid_points[1449].z =   5.7260617296; grid_points[1449].w_fixed =    0.1792877615; grid_points[1449].w_total = 0.0907093839; grid_points[1449].i_atom = 1;
    grid_points[1450].x =   1.5308700368; grid_points[1450].y =   0.0000000000; grid_points[1450].z =   5.0291191333; grid_points[1450].w_fixed =    0.1792877615; grid_points[1450].w_total = 0.1606309925; grid_points[1450].i_atom = 1;
    grid_points[1451].x =   1.5308700368; grid_points[1451].y =   0.0000000000; grid_points[1451].z =   3.3612642521; grid_points[1451].w_fixed =    0.1792877615; grid_points[1451].w_total = 0.1158377808; grid_points[1451].i_atom = 1;
    grid_points[1452].x =   2.3115565574; grid_points[1452].y =   0.0000000000; grid_points[1452].z =   4.1951916927; grid_points[1452].w_fixed =    0.1705811414; grid_points[1452].w_total = 0.1510163491; grid_points[1452].i_atom = 1;
    grid_points[1453].x =   2.2309387920; grid_points[1453].y =  -0.4279052608; grid_points[1453].z =   4.6230969535; grid_points[1453].w_fixed =    0.3659008737; grid_points[1453].w_total = 0.3324753846; grid_points[1453].i_atom = 1;
    grid_points[1454].x =   2.2309387920; grid_points[1454].y =  -0.4279052608; grid_points[1454].z =   3.7672864320; grid_points[1454].w_fixed =    0.3659008737; grid_points[1454].w_total = 0.2578929044; grid_points[1454].i_atom = 1;
    grid_points[1455].x =   1.5959473017; grid_points[1455].y =  -0.4991974941; grid_points[1455].z =   5.7911389945; grid_points[1455].w_fixed =    0.4430347113; grid_points[1455].w_total = 0.1354214943; grid_points[1455].i_atom = 1;
    grid_points[1456].x =   0.9146585961; grid_points[1456].y =  -0.9146585961; grid_points[1456].z =   6.1109403775; grid_points[1456].w_fixed =    0.4275576805; grid_points[1456].w_total = 0.0649370027; grid_points[1456].i_atom = 1;
    grid_points[1457].x =   1.9157486847; grid_points[1457].y =  -0.9146585961; grid_points[1457].z =   5.1098502889; grid_points[1457].w_fixed =    0.4275576805; grid_points[1457].w_total = 0.3174678581; grid_points[1457].i_atom = 1;
    grid_points[1458].x =   1.9157486847; grid_points[1458].y =  -0.9146585961; grid_points[1458].z =   3.2805330966; grid_points[1458].w_fixed =    0.4275576805; grid_points[1458].w_total = 0.1674035783; grid_points[1458].i_atom = 1;
    grid_points[1459].x =   2.0299139882; grid_points[1459].y =  -1.1057770654; grid_points[1459].z =   4.1951916927; grid_points[1459].w_fixed =    0.4319923443; grid_points[1459].w_total = 0.3838311380; grid_points[1459].i_atom = 1;
    grid_points[1460].x =   1.1057770654; grid_points[1460].y =   0.0000000000; grid_points[1460].z =   6.2251056810; grid_points[1460].w_fixed =    0.4319923443; grid_points[1460].w_total = 0.0483760601; grid_points[1460].i_atom = 1;
    grid_points[1461].x =   2.0299139882; grid_points[1461].y =   0.0000000000; grid_points[1461].z =   5.3009687581; grid_points[1461].w_fixed =    0.4319923443; grid_points[1461].w_total = 0.2694163787; grid_points[1461].i_atom = 1;
    grid_points[1462].x =   2.5684830011; grid_points[1462].y =  -1.2263012757; grid_points[1462].z =   5.4214929684; grid_points[1462].w_fixed =    1.0767762850; grid_points[1462].w_total = 0.2278386018; grid_points[1462].i_atom = 1;
    grid_points[1463].x =   2.7215467321; grid_points[1463].y =   0.0000000000; grid_points[1463].z =   5.6777293649; grid_points[1463].w_fixed =    1.0879446980; grid_points[1463].w_total = 0.1296664136; grid_points[1463].i_atom = 1;
    grid_points[1464].x =   2.9827802928; grid_points[1464].y =  -2.9827802928; grid_points[1464].z =   4.1951916927; grid_points[1464].w_fixed =    6.7598180997; grid_points[1464].w_total = 0.0643528246; grid_points[1464].i_atom = 1;
    grid_points[1465].x =   3.0991506199; grid_points[1465].y =   0.0000000000; grid_points[1465].z =   4.1951916927; grid_points[1465].w_fixed =    0.4295975400; grid_points[1465].w_total = 0.2078856501; grid_points[1465].i_atom = 1;
    grid_points[1466].x =   2.9910647516; grid_points[1466].y =  -0.5737012360; grid_points[1466].z =   4.7688929288; grid_points[1466].w_fixed =    0.9214976165; grid_points[1466].w_total = 0.4792995334; grid_points[1466].i_atom = 1;
    grid_points[1467].x =   2.9910647516; grid_points[1467].y =  -0.5737012360; grid_points[1467].z =   3.6214904567; grid_points[1467].w_fixed =    0.9214976165; grid_points[1467].w_total = 0.1977834187; grid_points[1467].i_atom = 1;
    grid_points[1468].x =  -7.5660149734; grid_points[1468].y =  -2.5220049911; grid_points[1468].z =   6.7171966839; grid_points[1468].w_fixed =   55.2891360098; grid_points[1468].w_total = 0.0000000000; grid_points[1468].i_atom = 1;
    grid_points[1469].x =   1.7666549613; grid_points[1469].y =  -1.7666549613; grid_points[1469].z =  -1.1047731912; grid_points[1469].w_fixed =   17.3848280115; grid_points[1469].w_total = 0.0000000000; grid_points[1469].i_atom = 1;
    grid_points[1470].x =   0.6692841749; grid_points[1470].y =  -2.1397188201; grid_points[1470].z =   6.3349105128; grid_points[1470].w_fixed =    1.1157541830; grid_points[1470].w_total = 0.0256935379; grid_points[1470].i_atom = 1;
    grid_points[1471].x =   2.4354299108; grid_points[1471].y =  -2.4354299108; grid_points[1471].z =   6.6306216035; grid_points[1471].w_fixed =    6.3162875542; grid_points[1471].w_total = 0.0011258350; grid_points[1471].i_atom = 1;
    grid_points[1472].x =   1.2718617903; grid_points[1472].y =  -1.2718617903; grid_points[1472].z =   8.0107770637; grid_points[1472].w_fixed =    6.0406797362; grid_points[1472].w_total = 0.0000010065; grid_points[1472].i_atom = 1;
    grid_points[1473].x =   0.4279052608; grid_points[1473].y =  -0.4279052608; grid_points[1473].z =   6.4261304847; grid_points[1473].w_fixed =    0.3659008737; grid_points[1473].w_total = 0.0221696713; grid_points[1473].i_atom = 1;
    grid_points[1474].x =   2.1397188201; grid_points[1474].y =  -0.6692841749; grid_points[1474].z =   6.3349105128; grid_points[1474].w_fixed =    1.1157541830; grid_points[1474].w_total = 0.0256935379; grid_points[1474].i_atom = 1;
    grid_points[1475].x =   0.5737012360; grid_points[1475].y =  -0.5737012360; grid_points[1475].z =   7.1862564443; grid_points[1475].w_fixed =    0.9214976165; grid_points[1475].w_total = 0.0005834114; grid_points[1475].i_atom = 1;
    grid_points[1476].x =   1.2263012757; grid_points[1476].y =  -1.2263012757; grid_points[1476].z =   6.7636746939; grid_points[1476].w_fixed =    1.0767762850; grid_points[1476].w_total = 0.0058282274; grid_points[1476].i_atom = 1;
    grid_points[1477].x =   1.4825376722; grid_points[1477].y =   0.0000000000; grid_points[1477].z =   6.9167384249; grid_points[1477].w_fixed =    1.0879446980; grid_points[1477].w_total = 0.0030241246; grid_points[1477].i_atom = 1;
    grid_points[1478].x =   2.5220049911; grid_points[1478].y =   7.5660149734; grid_points[1478].z =   6.7171966839; grid_points[1478].w_fixed =   55.2891360098; grid_points[1478].w_total = 0.0000000000; grid_points[1478].i_atom = 1;
    grid_points[1479].x =   2.9827802928; grid_points[1479].y =   0.0000000000; grid_points[1479].z =   7.1779719856; grid_points[1479].w_fixed =    6.7598180997; grid_points[1479].w_total = 0.0001667045; grid_points[1479].i_atom = 1;
    grid_points[1480].x =  -5.2999648839; grid_points[1480].y =   1.7666549613; grid_points[1480].z =   2.4285367314; grid_points[1480].w_fixed =   17.3848280115; grid_points[1480].w_total = 0.0000000010; grid_points[1480].i_atom = 1;
    grid_points[1481].x =  -3.8155853710; grid_points[1481].y =   1.2718617903; grid_points[1481].z =   2.9233299024; grid_points[1481].w_fixed =    6.0406797362; grid_points[1481].w_total = 0.0032694313; grid_points[1481].i_atom = 1;
    grid_points[1482].x =   5.2999648839; grid_points[1482].y =   1.7666549613; grid_points[1482].z =   2.4285367314; grid_points[1482].w_fixed =   17.3848280115; grid_points[1482].w_total = 0.0000000010; grid_points[1482].i_atom = 1;
    grid_points[1483].x =   3.8155853710; grid_points[1483].y =   1.2718617903; grid_points[1483].z =   2.9233299024; grid_points[1483].w_fixed =    6.0406797362; grid_points[1483].w_total = 0.0032694313; grid_points[1483].i_atom = 1;
    grid_points[1484].x =  -5.2999648839; grid_points[1484].y =   1.7666549613; grid_points[1484].z =   5.9618466541; grid_points[1484].w_fixed =   17.3848280115; grid_points[1484].w_total = 0.0000000058; grid_points[1484].i_atom = 1;
    grid_points[1485].x =  -3.8155853710; grid_points[1485].y =   1.2718617903; grid_points[1485].z =   5.4670534831; grid_points[1485].w_fixed =    6.0406797362; grid_points[1485].w_total = 0.0167747560; grid_points[1485].i_atom = 1;
    grid_points[1486].x =   3.8155853710; grid_points[1486].y =   1.2718617903; grid_points[1486].z =   5.4670534831; grid_points[1486].w_fixed =    6.0406797362; grid_points[1486].w_total = 0.0167747560; grid_points[1486].i_atom = 1;
    grid_points[1487].x =   5.2999648839; grid_points[1487].y =   1.7666549613; grid_points[1487].z =   5.9618466541; grid_points[1487].w_fixed =   17.3848280115; grid_points[1487].w_total = 0.0000000058; grid_points[1487].i_atom = 1;
    grid_points[1488].x =  -5.9146259784; grid_points[1488].y =   0.0000000000; grid_points[1488].z =  -1.7194342857; grid_points[1488].w_fixed =   61.8712659234; grid_points[1488].w_total = 0.0000000000; grid_points[1488].i_atom = 1;
    grid_points[1489].x =   5.9146259784; grid_points[1489].y =   0.0000000000; grid_points[1489].z =  -1.7194342857; grid_points[1489].w_fixed =   61.8712659234; grid_points[1489].w_total = 0.0000000000; grid_points[1489].i_atom = 1;
    grid_points[1490].x =  -3.3828867000; grid_points[1490].y =  -3.3828867000; grid_points[1490].z =   0.8123049927; grid_points[1490].w_fixed =   18.1780159844; grid_points[1490].w_total = 0.0000000000; grid_points[1490].i_atom = 1;
    grid_points[1491].x =  -3.3828867000; grid_points[1491].y =   3.3828867000; grid_points[1491].z =   0.8123049927; grid_points[1491].w_fixed =   18.1780159844; grid_points[1491].w_total = 0.0000000000; grid_points[1491].i_atom = 1;
    grid_points[1492].x =   3.3828867000; grid_points[1492].y =  -3.3828867000; grid_points[1492].z =   0.8123049927; grid_points[1492].w_fixed =   18.1780159844; grid_points[1492].w_total = 0.0000000000; grid_points[1492].i_atom = 1;
    grid_points[1493].x =   3.3828867000; grid_points[1493].y =   3.3828867000; grid_points[1493].z =   0.8123049927; grid_points[1493].w_fixed =   18.1780159844; grid_points[1493].w_total = 0.0000000000; grid_points[1493].i_atom = 1;
    grid_points[1494].x =  -5.9146259784; grid_points[1494].y =  -5.9146259784; grid_points[1494].z =   4.1951916927; grid_points[1494].w_fixed =   61.8712659234; grid_points[1494].w_total = 0.0000000000; grid_points[1494].i_atom = 1;
    grid_points[1495].x =  -5.9146259784; grid_points[1495].y =   5.9146259784; grid_points[1495].z =   4.1951916927; grid_points[1495].w_fixed =   61.8712659234; grid_points[1495].w_total = 0.0000000000; grid_points[1495].i_atom = 1;
    grid_points[1496].x =  -4.1431731363; grid_points[1496].y =  -4.1431731363; grid_points[1496].z =   4.1951916927; grid_points[1496].w_fixed =   19.4544786654; grid_points[1496].w_total = 0.0000002660; grid_points[1496].i_atom = 1;
    grid_points[1497].x =  -4.1431731363; grid_points[1497].y =   4.1431731363; grid_points[1497].z =   4.1951916927; grid_points[1497].w_fixed =   19.4544786654; grid_points[1497].w_total = 0.0000002660; grid_points[1497].i_atom = 1;
    grid_points[1498].x =   5.9146259784; grid_points[1498].y =  -5.9146259784; grid_points[1498].z =   4.1951916927; grid_points[1498].w_fixed =   61.8712659234; grid_points[1498].w_total = 0.0000000000; grid_points[1498].i_atom = 1;
    grid_points[1499].x =   4.1431731363; grid_points[1499].y =  -4.1431731363; grid_points[1499].z =   4.1951916927; grid_points[1499].w_fixed =   19.4544786654; grid_points[1499].w_total = 0.0000002660; grid_points[1499].i_atom = 1;
    grid_points[1500].x =   4.1431731363; grid_points[1500].y =   4.1431731363; grid_points[1500].z =   4.1951916927; grid_points[1500].w_fixed =   19.4544786654; grid_points[1500].w_total = 0.0000002660; grid_points[1500].i_atom = 1;
    grid_points[1501].x =   5.9146259784; grid_points[1501].y =   5.9146259784; grid_points[1501].z =   4.1951916927; grid_points[1501].w_fixed =   61.8712659234; grid_points[1501].w_total = 0.0000000000; grid_points[1501].i_atom = 1;
    grid_points[1502].x =  -4.8292718888; grid_points[1502].y =  -4.8292718888; grid_points[1502].z =  -0.6340801961; grid_points[1502].w_fixed =   57.8117193617; grid_points[1502].w_total = 0.0000000000; grid_points[1502].i_atom = 1;
    grid_points[1503].x =  -3.3828867000; grid_points[1503].y =  -3.3828867000; grid_points[1503].z =   7.5780783928; grid_points[1503].w_fixed =   18.1780159844; grid_points[1503].w_total = 0.0000000005; grid_points[1503].i_atom = 1;
    grid_points[1504].x =  -4.8292718888; grid_points[1504].y =  -4.8292718888; grid_points[1504].z =   9.0244635816; grid_points[1504].w_fixed =   57.8117193617; grid_points[1504].w_total = 0.0000000000; grid_points[1504].i_atom = 1;
    grid_points[1505].x =  -4.8292718888; grid_points[1505].y =   4.8292718888; grid_points[1505].z =  -0.6340801961; grid_points[1505].w_fixed =   57.8117193617; grid_points[1505].w_total = 0.0000000000; grid_points[1505].i_atom = 1;
    grid_points[1506].x =  -3.3828867000; grid_points[1506].y =   3.3828867000; grid_points[1506].z =   7.5780783928; grid_points[1506].w_fixed =   18.1780159844; grid_points[1506].w_total = 0.0000000005; grid_points[1506].i_atom = 1;
    grid_points[1507].x =  -4.8292718888; grid_points[1507].y =   4.8292718888; grid_points[1507].z =   9.0244635816; grid_points[1507].w_fixed =   57.8117193617; grid_points[1507].w_total = 0.0000000000; grid_points[1507].i_atom = 1;
    grid_points[1508].x =   4.8292718888; grid_points[1508].y =  -4.8292718888; grid_points[1508].z =  -0.6340801961; grid_points[1508].w_fixed =   57.8117193617; grid_points[1508].w_total = 0.0000000000; grid_points[1508].i_atom = 1;
    grid_points[1509].x =   3.3828867000; grid_points[1509].y =  -3.3828867000; grid_points[1509].z =   7.5780783928; grid_points[1509].w_fixed =   18.1780159844; grid_points[1509].w_total = 0.0000000005; grid_points[1509].i_atom = 1;
    grid_points[1510].x =   4.8292718888; grid_points[1510].y =  -4.8292718888; grid_points[1510].z =   9.0244635816; grid_points[1510].w_fixed =   57.8117193617; grid_points[1510].w_total = 0.0000000000; grid_points[1510].i_atom = 1;
    grid_points[1511].x =   4.8292718888; grid_points[1511].y =   4.8292718888; grid_points[1511].z =  -0.6340801961; grid_points[1511].w_fixed =   57.8117193617; grid_points[1511].w_total = 0.0000000000; grid_points[1511].i_atom = 1;
    grid_points[1512].x =   3.3828867000; grid_points[1512].y =   3.3828867000; grid_points[1512].z =   7.5780783928; grid_points[1512].w_fixed =   18.1780159844; grid_points[1512].w_total = 0.0000000005; grid_points[1512].i_atom = 1;
    grid_points[1513].x =   4.8292718888; grid_points[1513].y =   4.8292718888; grid_points[1513].z =   9.0244635816; grid_points[1513].w_fixed =   57.8117193617; grid_points[1513].w_total = 0.0000000000; grid_points[1513].i_atom = 1;
    grid_points[1514].x =  -5.2999648839; grid_points[1514].y =  -1.7666549613; grid_points[1514].z =   2.4285367314; grid_points[1514].w_fixed =   17.3848280115; grid_points[1514].w_total = 0.0000000010; grid_points[1514].i_atom = 1;
    grid_points[1515].x =  -5.9146259784; grid_points[1515].y =   0.0000000000; grid_points[1515].z =  10.1098176712; grid_points[1515].w_fixed =   61.8712659234; grid_points[1515].w_total = 0.0000000000; grid_points[1515].i_atom = 1;
    grid_points[1516].x =  -3.8155853710; grid_points[1516].y =  -1.2718617903; grid_points[1516].z =   2.9233299024; grid_points[1516].w_fixed =    6.0406797362; grid_points[1516].w_total = 0.0032694313; grid_points[1516].i_atom = 1;
    grid_points[1517].x =  -4.1431731363; grid_points[1517].y =   0.0000000000; grid_points[1517].z =   0.0520185564; grid_points[1517].w_fixed =   19.4544786654; grid_points[1517].w_total = 0.0000000000; grid_points[1517].i_atom = 1;
    grid_points[1518].x =   5.2999648839; grid_points[1518].y =  -1.7666549613; grid_points[1518].z =   2.4285367314; grid_points[1518].w_fixed =   17.3848280115; grid_points[1518].w_total = 0.0000000010; grid_points[1518].i_atom = 1;
    grid_points[1519].x =   3.8155853710; grid_points[1519].y =  -1.2718617903; grid_points[1519].z =   2.9233299024; grid_points[1519].w_fixed =    6.0406797362; grid_points[1519].w_total = 0.0032694313; grid_points[1519].i_atom = 1;
    grid_points[1520].x =   4.1431731363; grid_points[1520].y =   0.0000000000; grid_points[1520].z =   0.0520185564; grid_points[1520].w_fixed =   19.4544786654; grid_points[1520].w_total = 0.0000000000; grid_points[1520].i_atom = 1;
    grid_points[1521].x =   5.9146259784; grid_points[1521].y =   0.0000000000; grid_points[1521].z =  10.1098176712; grid_points[1521].w_fixed =   61.8712659234; grid_points[1521].w_total = 0.0000000000; grid_points[1521].i_atom = 1;
    grid_points[1522].x =  -5.2999648839; grid_points[1522].y =  -1.7666549613; grid_points[1522].z =   5.9618466541; grid_points[1522].w_fixed =   17.3848280115; grid_points[1522].w_total = 0.0000000058; grid_points[1522].i_atom = 1;
    grid_points[1523].x =  -5.8593316407; grid_points[1523].y =   0.0000000000; grid_points[1523].z =   4.1951916927; grid_points[1523].w_fixed =   10.9431442493; grid_points[1523].w_total = 0.0000000054; grid_points[1523].i_atom = 1;
    grid_points[1524].x =  -3.8155853710; grid_points[1524].y =  -1.2718617903; grid_points[1524].z =   5.4670534831; grid_points[1524].w_fixed =    6.0406797362; grid_points[1524].w_total = 0.0167747560; grid_points[1524].i_atom = 1;
    grid_points[1525].x =  -4.2182883437; grid_points[1525].y =   0.0000000000; grid_points[1525].z =   4.1951916927; grid_points[1525].w_fixed =    3.8023976811; grid_points[1525].w_total = 0.0216466137; grid_points[1525].i_atom = 1;
    grid_points[1526].x =   3.8155853710; grid_points[1526].y =  -1.2718617903; grid_points[1526].z =   5.4670534831; grid_points[1526].w_fixed =    6.0406797362; grid_points[1526].w_total = 0.0167747560; grid_points[1526].i_atom = 1;
    grid_points[1527].x =   5.2999648839; grid_points[1527].y =  -1.7666549613; grid_points[1527].z =   5.9618466541; grid_points[1527].w_fixed =   17.3848280115; grid_points[1527].w_total = 0.0000000058; grid_points[1527].i_atom = 1;
    grid_points[1528].x =   4.2182883437; grid_points[1528].y =   0.0000000000; grid_points[1528].z =   4.1951916927; grid_points[1528].w_fixed =    3.8023976811; grid_points[1528].w_total = 0.0216466137; grid_points[1528].i_atom = 1;
    grid_points[1529].x =   5.8593316407; grid_points[1529].y =   0.0000000000; grid_points[1529].z =   4.1951916927; grid_points[1529].w_fixed =   10.9431442493; grid_points[1529].w_total = 0.0000000054; grid_points[1529].i_atom = 1;
    grid_points[1530].x =  -4.1431731363; grid_points[1530].y =   0.0000000000; grid_points[1530].z =   8.3383648291; grid_points[1530].w_fixed =   19.4544786654; grid_points[1530].w_total = 0.0000000000; grid_points[1530].i_atom = 1;
    grid_points[1531].x =   4.1431731363; grid_points[1531].y =   0.0000000000; grid_points[1531].z =   8.3383648291; grid_points[1531].w_fixed =   19.4544786654; grid_points[1531].w_total = 0.0000000000; grid_points[1531].i_atom = 1;
    grid_points[1532].x =  -2.5220049911; grid_points[1532].y =  -7.5660149734; grid_points[1532].z =   1.6731867016; grid_points[1532].w_fixed =   55.2891360098; grid_points[1532].w_total = 0.0000000000; grid_points[1532].i_atom = 1;
    grid_points[1533].x =  -1.5959473017; grid_points[1533].y =   0.4991974941; grid_points[1533].z =   2.5992443910; grid_points[1533].w_fixed =    0.4430347113; grid_points[1533].w_total = 0.0354740456; grid_points[1533].i_atom = 1;
    grid_points[1534].x =  -1.7892954446; grid_points[1534].y =   1.7892954446; grid_points[1534].z =   2.4058962481; grid_points[1534].w_fixed =    1.0990251469; grid_points[1534].w_total = 0.0123810957; grid_points[1534].i_atom = 1;
    grid_points[1535].x =  -2.1397188201; grid_points[1535].y =   0.6692841749; grid_points[1535].z =   2.0554728727; grid_points[1535].w_fixed =    1.1157541830; grid_points[1535].w_total = 0.0031940004; grid_points[1535].i_atom = 1;
    grid_points[1536].x =  -2.4354299108; grid_points[1536].y =   2.4354299108; grid_points[1536].z =   1.7597617820; grid_points[1536].w_fixed =    6.3162875542; grid_points[1536].w_total = 0.0001328651; grid_points[1536].i_atom = 1;
    grid_points[1537].x =  -1.2718617903; grid_points[1537].y =   1.2718617903; grid_points[1537].z =   0.3796063218; grid_points[1537].w_fixed =    6.0406797362; grid_points[1537].w_total = 0.0000000021; grid_points[1537].i_atom = 1;
    grid_points[1538].x =  -1.3345778006; grid_points[1538].y =   1.3345778006; grid_points[1538].z =   2.8606138921; grid_points[1538].w_fixed =    0.4363920800; grid_points[1538].w_total = 0.0730945598; grid_points[1538].i_atom = 1;
    grid_points[1539].x =  -2.5684830011; grid_points[1539].y =   1.2263012757; grid_points[1539].z =   2.9688904170; grid_points[1539].w_fixed =    1.0767762850; grid_points[1539].w_total = 0.0592143299; grid_points[1539].i_atom = 1;
    grid_points[1540].x =  -1.7666549613; grid_points[1540].y =   1.7666549613; grid_points[1540].z =   9.4951565767; grid_points[1540].w_fixed =   17.3848280115; grid_points[1540].w_total = 0.0000000000; grid_points[1540].i_atom = 1;
    grid_points[1541].x =  -2.5220049911; grid_points[1541].y =   2.5220049911; grid_points[1541].z =  11.7612066662; grid_points[1541].w_fixed =   55.2891360098; grid_points[1541].w_total = 0.0000000000; grid_points[1541].i_atom = 1;
    grid_points[1542].x =  -0.3227069453; grid_points[1542].y =   0.3227069453; grid_points[1542].z =   2.5127177700; grid_points[1542].w_fixed =    0.1518581277; grid_points[1542].w_total = 0.0160452730; grid_points[1542].i_atom = 1;
    grid_points[1543].x =   0.0000000000; grid_points[1543].y =   0.8339274406; grid_points[1543].z =   2.6643216559; grid_points[1543].w_fixed =    0.1792877615; grid_points[1543].w_total = 0.0310706777; grid_points[1543].i_atom = 1;
    grid_points[1544].x =  -0.4279052608; grid_points[1544].y =   0.4279052608; grid_points[1544].z =   1.9642529008; grid_points[1544].w_fixed =    0.3659008737; grid_points[1544].w_total = 0.0015586987; grid_points[1544].i_atom = 1;
    grid_points[1545].x =  -0.4991974941; grid_points[1545].y =   1.5959473017; grid_points[1545].z =   2.5992443910; grid_points[1545].w_fixed =    0.4430347113; grid_points[1545].w_total = 0.0354740456; grid_points[1545].i_atom = 1;
    grid_points[1546].x =  -0.9146585961; grid_points[1546].y =   0.9146585961; grid_points[1546].z =   2.2794430080; grid_points[1546].w_fixed =    0.4275576805; grid_points[1546].w_total = 0.0104656786; grid_points[1546].i_atom = 1;
    grid_points[1547].x =   0.0000000000; grid_points[1547].y =   1.1057770654; grid_points[1547].z =   2.1652777045; grid_points[1547].w_fixed =    0.4319923443; grid_points[1547].w_total = 0.0061142765; grid_points[1547].i_atom = 1;
    grid_points[1548].x =  -0.5737012360; grid_points[1548].y =   0.5737012360; grid_points[1548].z =   1.2041269412; grid_points[1548].w_fixed =    0.9214976165; grid_points[1548].w_total = 0.0000026829; grid_points[1548].i_atom = 1;
    grid_points[1549].x =  -0.6692841749; grid_points[1549].y =   2.1397188201; grid_points[1549].z =   2.0554728727; grid_points[1549].w_fixed =    1.1157541830; grid_points[1549].w_total = 0.0031940004; grid_points[1549].i_atom = 1;
    grid_points[1550].x =  -1.2263012757; grid_points[1550].y =   1.2263012757; grid_points[1550].z =   1.6267086916; grid_points[1550].w_fixed =    1.0767762850; grid_points[1550].w_total = 0.0002888020; grid_points[1550].i_atom = 1;
    grid_points[1551].x =   0.0000000000; grid_points[1551].y =   1.4825376722; grid_points[1551].z =   1.4736449606; grid_points[1551].w_fixed =    1.0879446980; grid_points[1551].w_total = 0.0000849175; grid_points[1551].i_atom = 1;
    grid_points[1552].x =   0.0000000000; grid_points[1552].y =   2.7215467321; grid_points[1552].z =   2.7126540206; grid_points[1552].w_fixed =    1.0879446980; grid_points[1552].w_total = 0.0291805697; grid_points[1552].i_atom = 1;
    grid_points[1553].x =  -0.2452657651; grid_points[1553].y =   0.2452657651; grid_points[1553].z =   2.9164672649; grid_points[1553].w_fixed =    0.0651606456; grid_points[1553].w_total = 0.0288023961; grid_points[1553].i_atom = 1;
    grid_points[1554].x =  -0.5242619359; grid_points[1554].y =   0.5242619359; grid_points[1554].z =   3.0971272145; grid_points[1554].w_fixed =    0.0761406613; grid_points[1554].w_total = 0.0448544141; grid_points[1554].i_atom = 1;
    grid_points[1555].x =   0.0000000000; grid_points[1555].y =   0.6338067859; grid_points[1555].z =   3.0316902052; grid_points[1555].w_fixed =    0.0769303986; grid_points[1555].w_total = 0.0412765861; grid_points[1555].i_atom = 1;
    grid_points[1556].x =  -1.2035918363; grid_points[1556].y =   0.3764723484; grid_points[1556].z =   2.9915998565; grid_points[1556].w_fixed =    0.1838706234; grid_points[1556].w_total = 0.0694127716; grid_points[1556].i_atom = 1;
    grid_points[1557].x =  -0.3764723484; grid_points[1557].y =   1.2035918363; grid_points[1557].z =   2.9915998565; grid_points[1557].w_fixed =    0.1838706234; grid_points[1557].w_total = 0.0694127716; grid_points[1557].i_atom = 1;
    grid_points[1558].x =  -0.6897944676; grid_points[1558].y =   0.6897944676; grid_points[1558].z =   2.7504200046; grid_points[1558].w_fixed =    0.1774472638; grid_points[1558].w_total = 0.0390890156; grid_points[1558].i_atom = 1;
    grid_points[1559].x =   0.0000000000; grid_points[1559].y =   2.0299139882; grid_points[1559].z =   3.0894146273; grid_points[1559].w_fixed =    0.4319923443; grid_points[1559].w_total = 0.1180180438; grid_points[1559].i_atom = 1;
    grid_points[1560].x =  -1.2263012757; grid_points[1560].y =   2.5684830011; grid_points[1560].z =   2.9688904170; grid_points[1560].w_fixed =    1.0767762850; grid_points[1560].w_total = 0.0592143299; grid_points[1560].i_atom = 1;
    grid_points[1561].x =   0.0000000000; grid_points[1561].y =   2.9827802928; grid_points[1561].z =   1.2124113999; grid_points[1561].w_fixed =    6.7598180997; grid_points[1561].w_total = 0.0000096089; grid_points[1561].i_atom = 1;
    grid_points[1562].x =   7.5660149734; grid_points[1562].y =   2.5220049911; grid_points[1562].z =   1.6731867016; grid_points[1562].w_fixed =   55.2891360098; grid_points[1562].w_total = 0.0000000000; grid_points[1562].i_atom = 1;
    grid_points[1563].x =  -2.5220049911; grid_points[1563].y =   2.5220049911; grid_points[1563].z =  -3.3708232807; grid_points[1563].w_fixed =   55.2891360098; grid_points[1563].w_total = 0.0000000000; grid_points[1563].i_atom = 1;
    grid_points[1564].x =  -1.2787244278; grid_points[1564].y =   0.2452657651; grid_points[1564].z =   4.4404574578; grid_points[1564].w_fixed =    0.0651606456; grid_points[1564].w_total = 0.0648319178; grid_points[1564].i_atom = 1;
    grid_points[1565].x =  -1.2787244278; grid_points[1565].y =   0.2452657651; grid_points[1565].z =   3.9499259276; grid_points[1565].w_fixed =    0.0651606456; grid_points[1565].w_total = 0.0634126219; grid_points[1565].i_atom = 1;
    grid_points[1566].x =  -1.6824739228; grid_points[1566].y =   0.3227069453; grid_points[1566].z =   4.5178986380; grid_points[1566].w_fixed =    0.1518581277; grid_points[1566].w_total = 0.1485902764; grid_points[1566].i_atom = 1;
    grid_points[1567].x =  -1.6824739228; grid_points[1567].y =   0.3227069453; grid_points[1567].z =   3.8724847475; grid_points[1567].w_fixed =    0.1518581277; grid_points[1567].w_total = 0.1376185851; grid_points[1567].i_atom = 1;
    grid_points[1568].x =  -1.4447716881; grid_points[1568].y =   0.6897944676; grid_points[1568].z =   4.8849861603; grid_points[1568].w_fixed =    0.1774472638; grid_points[1568].w_total = 0.1658131153; grid_points[1568].i_atom = 1;
    grid_points[1569].x =  -1.4447716881; grid_points[1569].y =   0.6897944676; grid_points[1569].z =   3.5053972252; grid_points[1569].w_fixed =    0.1774472638; grid_points[1569].w_total = 0.1314016342; grid_points[1569].i_atom = 1;
    grid_points[1570].x =  -1.5308700368; grid_points[1570].y =   0.8339274406; grid_points[1570].z =   4.1951916927; grid_points[1570].w_fixed =    0.1792877615; grid_points[1570].w_total = 0.1738017119; grid_points[1570].i_atom = 1;
    grid_points[1571].x =  -1.3345778006; grid_points[1571].y =   1.3345778006; grid_points[1571].z =   5.5297694934; grid_points[1571].w_fixed =    0.4363920800; grid_points[1571].w_total = 0.2072585535; grid_points[1571].i_atom = 1;
    grid_points[1572].x =  -2.2309387920; grid_points[1572].y =   0.4279052608; grid_points[1572].z =   4.6230969535; grid_points[1572].w_fixed =    0.3659008737; grid_points[1572].w_total = 0.3324753846; grid_points[1572].i_atom = 1;
    grid_points[1573].x =  -2.2309387920; grid_points[1573].y =   0.4279052608; grid_points[1573].z =   3.7672864320; grid_points[1573].w_fixed =    0.3659008737; grid_points[1573].w_total = 0.2578929044; grid_points[1573].i_atom = 1;
    grid_points[1574].x =  -1.5959473017; grid_points[1574].y =   1.5959473017; grid_points[1574].z =   4.6943891868; grid_points[1574].w_fixed =    0.4430347113; grid_points[1574].w_total = 0.3987026542; grid_points[1574].i_atom = 1;
    grid_points[1575].x =  -1.5959473017; grid_points[1575].y =   1.5959473017; grid_points[1575].z =   3.6959941987; grid_points[1575].w_fixed =    0.4430347113; grid_points[1575].w_total = 0.2967672396; grid_points[1575].i_atom = 1;
    grid_points[1576].x =  -1.5959473017; grid_points[1576].y =   0.4991974941; grid_points[1576].z =   5.7911389945; grid_points[1576].w_fixed =    0.4430347113; grid_points[1576].w_total = 0.1354214943; grid_points[1576].i_atom = 1;
    grid_points[1577].x =  -1.9157486847; grid_points[1577].y =   0.9146585961; grid_points[1577].z =   5.1098502889; grid_points[1577].w_fixed =    0.4275576805; grid_points[1577].w_total = 0.3174678581; grid_points[1577].i_atom = 1;
    grid_points[1578].x =  -1.9157486847; grid_points[1578].y =   0.9146585961; grid_points[1578].z =   3.2805330966; grid_points[1578].w_fixed =    0.4275576805; grid_points[1578].w_total = 0.1674035783; grid_points[1578].i_atom = 1;
    grid_points[1579].x =  -2.0299139882; grid_points[1579].y =   1.1057770654; grid_points[1579].z =   4.1951916927; grid_points[1579].w_fixed =    0.4319923443; grid_points[1579].w_total = 0.3838311380; grid_points[1579].i_atom = 1;
    grid_points[1580].x =  -1.7892954446; grid_points[1580].y =   1.7892954446; grid_points[1580].z =   5.9844871374; grid_points[1580].w_fixed =    1.0990251469; grid_points[1580].w_total = 0.0671861197; grid_points[1580].i_atom = 1;
    grid_points[1581].x =  -2.9910647516; grid_points[1581].y =   0.5737012360; grid_points[1581].z =   4.7688929288; grid_points[1581].w_fixed =    0.9214976165; grid_points[1581].w_total = 0.4792995334; grid_points[1581].i_atom = 1;
    grid_points[1582].x =  -2.9910647516; grid_points[1582].y =   0.5737012360; grid_points[1582].z =   3.6214904567; grid_points[1582].w_fixed =    0.9214976165; grid_points[1582].w_total = 0.1977834187; grid_points[1582].i_atom = 1;
    grid_points[1583].x =  -2.1397188201; grid_points[1583].y =   2.1397188201; grid_points[1583].z =   4.8644758676; grid_points[1583].w_fixed =    1.1157541830; grid_points[1583].w_total = 0.5622597816; grid_points[1583].i_atom = 1;
    grid_points[1584].x =  -2.1397188201; grid_points[1584].y =   2.1397188201; grid_points[1584].z =   3.5259075179; grid_points[1584].w_fixed =    1.1157541830; grid_points[1584].w_total = 0.2172616206; grid_points[1584].i_atom = 1;
    grid_points[1585].x =  -2.5684830011; grid_points[1585].y =   1.2263012757; grid_points[1585].z =   5.4214929684; grid_points[1585].w_fixed =    1.0767762850; grid_points[1585].w_total = 0.2278386018; grid_points[1585].i_atom = 1;
    grid_points[1586].x =  -1.4825376722; grid_points[1586].y =   2.7215467321; grid_points[1586].z =   4.1951916927; grid_points[1586].w_fixed =    1.0879446980; grid_points[1586].w_total = 0.5523500480; grid_points[1586].i_atom = 1;
    grid_points[1587].x =  -2.7215467321; grid_points[1587].y =   1.4825376722; grid_points[1587].z =   4.1951916927; grid_points[1587].w_fixed =    1.0879446980; grid_points[1587].w_total = 0.5523500480; grid_points[1587].i_atom = 1;
    grid_points[1588].x =  -2.9827802928; grid_points[1588].y =   2.9827802928; grid_points[1588].z =   4.1951916927; grid_points[1588].w_fixed =    6.7598180997; grid_points[1588].w_total = 0.0643528246; grid_points[1588].i_atom = 1;
    grid_points[1589].x =   0.0000000000; grid_points[1589].y =   0.0009212695; grid_points[1589].z =   4.1951916927; grid_points[1589].w_fixed =    0.0000000034; grid_points[1589].w_total = 0.0000000034; grid_points[1589].i_atom = 1;
    grid_points[1590].x =   0.0000000000; grid_points[1590].y =   0.0039529982; grid_points[1590].z =   4.1951916927; grid_points[1590].w_fixed =    0.0000001386; grid_points[1590].w_total = 0.0000001386; grid_points[1590].i_atom = 1;
    grid_points[1591].x =   0.0000000000; grid_points[1591].y =   0.0095652797; grid_points[1591].z =   4.1951916927; grid_points[1591].w_fixed =    0.0000013577; grid_points[1591].w_total = 0.0000013577; grid_points[1591].i_atom = 1;
    grid_points[1592].x =   0.0000000000; grid_points[1592].y =   0.0183381693; grid_points[1592].z =   4.1951916927; grid_points[1592].w_fixed =    0.0000074515; grid_points[1592].w_total = 0.0000074515; grid_points[1592].i_atom = 1;
    grid_points[1593].x =   0.0000000000; grid_points[1593].y =   0.0309915062; grid_points[1593].z =   4.1951916927; grid_points[1593].w_fixed =    0.0000299246; grid_points[1593].w_total = 0.0000299246; grid_points[1593].i_atom = 1;
    grid_points[1594].x =   0.0000000000; grid_points[1594].y =   0.0484242284; grid_points[1594].z =   4.1951916927; grid_points[1594].w_fixed =    0.0000990913; grid_points[1594].w_total = 0.0000990913; grid_points[1594].i_atom = 1;
    grid_points[1595].x =   0.0000000000; grid_points[1595].y =   0.0717667204; grid_points[1595].z =   4.1951916927; grid_points[1595].w_fixed =    0.0002885050; grid_points[1595].w_total = 0.0002885050; grid_points[1595].i_atom = 1;
    grid_points[1596].x =   0.0000000000; grid_points[1596].y =   0.1024512602; grid_points[1596].z =   4.1951916927; grid_points[1596].w_fixed =    0.0007678018; grid_points[1596].w_total = 0.0007678018; grid_points[1596].i_atom = 1;
    grid_points[1597].x =   0.0000000000; grid_points[1597].y =   0.1423079366; grid_points[1597].z =   4.1951916927; grid_points[1597].w_fixed =    0.0005474789; grid_points[1597].w_total = 0.0005474789; grid_points[1597].i_atom = 1;
    grid_points[1598].x =   0.0000000000; grid_points[1598].y =   0.1006269070; grid_points[1598].z =   4.2958185997; grid_points[1598].w_fixed =    0.0004379831; grid_points[1598].w_total = 0.0004379831; grid_points[1598].i_atom = 1;
    grid_points[1599].x =   0.0000000000; grid_points[1599].y =   0.1006269070; grid_points[1599].z =   4.0945647857; grid_points[1599].w_fixed =    0.0004379831; grid_points[1599].w_total = 0.0004379831; grid_points[1599].i_atom = 1;
    grid_points[1600].x =  -0.1006269070; grid_points[1600].y =   0.1006269070; grid_points[1600].z =   4.1951916927; grid_points[1600].w_fixed =    0.0004379831; grid_points[1600].w_total = 0.0004379831; grid_points[1600].i_atom = 1;
    grid_points[1601].x =  -0.0821615255; grid_points[1601].y =   0.0821615255; grid_points[1601].z =   4.2773532183; grid_points[1601].w_fixed =    0.0003695483; grid_points[1601].w_total = 0.0003695483; grid_points[1601].i_atom = 1;
    grid_points[1602].x =  -0.0821615255; grid_points[1602].y =   0.0821615255; grid_points[1602].z =   4.1130301672; grid_points[1602].w_fixed =    0.0003695483; grid_points[1602].w_total = 0.0003695482; grid_points[1602].i_atom = 1;
    grid_points[1603].x =   0.0000000000; grid_points[1603].y =   0.1936969137; grid_points[1603].z =   4.1951916927; grid_points[1603].w_fixed =    0.0013046077; grid_points[1603].w_total = 0.0013046076; grid_points[1603].i_atom = 1;
    grid_points[1604].x =   0.0000000000; grid_points[1604].y =   0.1369644012; grid_points[1604].z =   4.3321560939; grid_points[1604].w_fixed =    0.0010436861; grid_points[1604].w_total = 0.0010436860; grid_points[1604].i_atom = 1;
    grid_points[1605].x =   0.0000000000; grid_points[1605].y =   0.1369644012; grid_points[1605].z =   4.0582272915; grid_points[1605].w_fixed =    0.0010436861; grid_points[1605].w_total = 0.0010436847; grid_points[1605].i_atom = 1;
    grid_points[1606].x =  -0.1369644012; grid_points[1606].y =   0.1369644012; grid_points[1606].z =   4.1951916927; grid_points[1606].w_fixed =    0.0010436861; grid_points[1606].w_total = 0.0010436861; grid_points[1606].i_atom = 1;
    grid_points[1607].x =  -0.1118309653; grid_points[1607].y =   0.1118309653; grid_points[1607].z =   4.3070226580; grid_points[1607].w_fixed =    0.0008806102; grid_points[1607].w_total = 0.0008806101; grid_points[1607].i_atom = 1;
    grid_points[1608].x =  -0.1118309653; grid_points[1608].y =   0.1118309653; grid_points[1608].z =   4.0833607275; grid_points[1608].w_fixed =    0.0008806102; grid_points[1608].w_total = 0.0008806095; grid_points[1608].i_atom = 1;
    grid_points[1609].x =   0.0000000000; grid_points[1609].y =   0.2596933691; grid_points[1609].z =   4.1951916927; grid_points[1609].w_fixed =    0.0030086947; grid_points[1609].w_total = 0.0030086941; grid_points[1609].i_atom = 1;
    grid_points[1610].x =   0.0000000000; grid_points[1610].y =   0.1836309423; grid_points[1610].z =   4.3788226351; grid_points[1610].w_fixed =    0.0024069557; grid_points[1610].w_total = 0.0024069533; grid_points[1610].i_atom = 1;
    grid_points[1611].x =   0.0000000000; grid_points[1611].y =   0.1836309423; grid_points[1611].z =   4.0115607504; grid_points[1611].w_fixed =    0.0024069557; grid_points[1611].w_total = 0.0024069249; grid_points[1611].i_atom = 1;
    grid_points[1612].x =  -0.1836309423; grid_points[1612].y =   0.1836309423; grid_points[1612].z =   4.1951916927; grid_points[1612].w_fixed =    0.0024069557; grid_points[1612].w_total = 0.0024069553; grid_points[1612].i_atom = 1;
    grid_points[1613].x =  -0.1499340366; grid_points[1613].y =   0.1499340366; grid_points[1613].z =   4.3451257293; grid_points[1613].w_fixed =    0.0020308689; grid_points[1613].w_total = 0.0020308678; grid_points[1613].i_atom = 1;
    grid_points[1614].x =  -0.1499340366; grid_points[1614].y =   0.1499340366; grid_points[1614].z =   4.0452576562; grid_points[1614].w_fixed =    0.0020308689; grid_points[1614].w_total = 0.0020308551; grid_points[1614].i_atom = 1;
    grid_points[1615].x =   0.0000000000; grid_points[1615].y =   0.3443500689; grid_points[1615].z =   4.1951916927; grid_points[1615].w_fixed =    0.0067871679; grid_points[1615].w_total = 0.0067871580; grid_points[1615].i_atom = 1;
    grid_points[1616].x =   0.0000000000; grid_points[1616].y =   0.2434922688; grid_points[1616].z =   4.4386839615; grid_points[1616].w_fixed =    0.0054297343; grid_points[1616].w_total = 0.0054296870; grid_points[1616].i_atom = 1;
    grid_points[1617].x =   0.0000000000; grid_points[1617].y =   0.2434922688; grid_points[1617].z =   3.9516994239; grid_points[1617].w_fixed =    0.0054297343; grid_points[1617].w_total = 0.0054291838; grid_points[1617].i_atom = 1;
    grid_points[1618].x =  -0.2434922688; grid_points[1618].y =   0.2434922688; grid_points[1618].z =   4.1951916927; grid_points[1618].w_fixed =    0.0054297343; grid_points[1618].w_total = 0.0054297264; grid_points[1618].i_atom = 1;
    grid_points[1619].x =  -0.1988106050; grid_points[1619].y =   0.1988106050; grid_points[1619].z =   4.3940022977; grid_points[1619].w_fixed =    0.0045813383; grid_points[1619].w_total = 0.0045813179; grid_points[1619].i_atom = 1;
    grid_points[1620].x =  -0.1988106050; grid_points[1620].y =   0.1988106050; grid_points[1620].z =   3.9963810878; grid_points[1620].w_fixed =    0.0045813383; grid_points[1620].w_total = 0.0045810919; grid_points[1620].i_atom = 1;
    grid_points[1621].x =   0.0000000000; grid_points[1621].y =   0.4530765180; grid_points[1621].z =   4.1951916927; grid_points[1621].w_fixed =    0.0040293310; grid_points[1621].w_total = 0.0040292893; grid_points[1621].i_atom = 1;
    grid_points[1622].x =   0.0000000000; grid_points[1622].y =   0.3203734782; grid_points[1622].z =   4.5155651710; grid_points[1622].w_fixed =    0.0071632550; grid_points[1622].w_total = 0.0071627669; grid_points[1622].i_atom = 1;
    grid_points[1623].x =   0.0000000000; grid_points[1623].y =   0.3203734782; grid_points[1623].z =   3.8748182145; grid_points[1623].w_fixed =    0.0071632550; grid_points[1623].w_total = 0.0071582150; grid_points[1623].i_atom = 1;
    grid_points[1624].x =  -0.3203734782; grid_points[1624].y =   0.3203734782; grid_points[1624].z =   4.1951916927; grid_points[1624].w_fixed =    0.0071632550; grid_points[1624].w_total = 0.0071631809; grid_points[1624].i_atom = 1;
    grid_points[1625].x =  -0.2615838496; grid_points[1625].y =   0.2615838496; grid_points[1625].z =   4.4567755423; grid_points[1625].w_fixed =    0.0066932539; grid_points[1625].w_total = 0.0066930222; grid_points[1625].i_atom = 1;
    grid_points[1626].x =  -0.2615838496; grid_points[1626].y =   0.2615838496; grid_points[1626].z =   3.9336078431; grid_points[1626].w_fixed =    0.0066932539; grid_points[1626].w_total = 0.0066907381; grid_points[1626].i_atom = 1;
    grid_points[1627].x =  -0.1366077101; grid_points[1627].y =   0.1366077101; grid_points[1627].z =   4.6050148231; grid_points[1627].w_fixed =    0.0064011973; grid_points[1627].w_total = 0.0064000836; grid_points[1627].i_atom = 1;
    grid_points[1628].x =  -0.1366077101; grid_points[1628].y =   0.1366077101; grid_points[1628].z =   3.7853685624; grid_points[1628].w_fixed =    0.0064011973; grid_points[1628].w_total = 0.0063905489; grid_points[1628].i_atom = 1;
    grid_points[1629].x =  -0.1366077101; grid_points[1629].y =   0.4098231304; grid_points[1629].z =   4.3317994029; grid_points[1629].w_fixed =    0.0064011973; grid_points[1629].w_total = 0.0064011513; grid_points[1629].i_atom = 1;
    grid_points[1630].x =  -0.1366077101; grid_points[1630].y =   0.4098231304; grid_points[1630].z =   4.0585839826; grid_points[1630].w_fixed =    0.0064011973; grid_points[1630].w_total = 0.0064006801; grid_points[1630].i_atom = 1;
    grid_points[1631].x =  -0.4098231304; grid_points[1631].y =   0.1366077101; grid_points[1631].z =   4.3317994029; grid_points[1631].w_fixed =    0.0064011973; grid_points[1631].w_total = 0.0064011513; grid_points[1631].i_atom = 1;
    grid_points[1632].x =  -0.4098231304; grid_points[1632].y =   0.1366077101; grid_points[1632].z =   4.0585839826; grid_points[1632].w_fixed =    0.0064011973; grid_points[1632].w_total = 0.0064006801; grid_points[1632].i_atom = 1;
    grid_points[1633].x =   0.0000000000; grid_points[1633].y =   0.5931967983; grid_points[1633].z =   4.1951916927; grid_points[1633].w_fixed =    0.0089219324; grid_points[1633].w_total = 0.0089213452; grid_points[1633].i_atom = 1;
    grid_points[1634].x =   0.0000000000; grid_points[1634].y =   0.4194534787; grid_points[1634].z =   4.6146451714; grid_points[1634].w_fixed =    0.0158612132; grid_points[1634].w_total = 0.0158535719; grid_points[1634].i_atom = 1;
    grid_points[1635].x =   0.0000000000; grid_points[1635].y =   0.4194534787; grid_points[1635].z =   3.7757382141; grid_points[1635].w_fixed =    0.0158612132; grid_points[1635].w_total = 0.0157937824; grid_points[1635].i_atom = 1;
    grid_points[1636].x =  -0.4194534787; grid_points[1636].y =   0.4194534787; grid_points[1636].z =   4.1951916927; grid_points[1636].w_fixed =    0.0158612132; grid_points[1636].w_total = 0.0158601692; grid_points[1636].i_atom = 1;
    grid_points[1637].x =  -0.3424823312; grid_points[1637].y =   0.3424823312; grid_points[1637].z =   4.5376740239; grid_points[1637].w_fixed =    0.0148205147; grid_points[1637].w_total = 0.0148168993; grid_points[1637].i_atom = 1;
    grid_points[1638].x =  -0.3424823312; grid_points[1638].y =   0.3424823312; grid_points[1638].z =   3.8527093615; grid_points[1638].w_fixed =    0.0148205147; grid_points[1638].w_total = 0.0147863672; grid_points[1638].i_atom = 1;
    grid_points[1639].x =  -0.1788555643; grid_points[1639].y =   0.1788555643; grid_points[1639].z =   4.7317583855; grid_points[1639].w_fixed =    0.0141738295; grid_points[1639].w_total = 0.0141564043; grid_points[1639].i_atom = 1;
    grid_points[1640].x =  -0.1788555643; grid_points[1640].y =   0.1788555643; grid_points[1640].z =   3.6586249999; grid_points[1640].w_fixed =    0.0141738295; grid_points[1640].w_total = 0.0140351874; grid_points[1640].i_atom = 1;
    grid_points[1641].x =  -0.1788555643; grid_points[1641].y =   0.5365666928; grid_points[1641].z =   4.3740472570; grid_points[1641].w_fixed =    0.0141738295; grid_points[1641].w_total = 0.0141731261; grid_points[1641].i_atom = 1;
    grid_points[1642].x =  -0.1788555643; grid_points[1642].y =   0.5365666928; grid_points[1642].z =   4.0163361285; grid_points[1642].w_fixed =    0.0141738295; grid_points[1642].w_total = 0.0141666532; grid_points[1642].i_atom = 1;
    grid_points[1643].x =  -0.5365666928; grid_points[1643].y =   0.1788555643; grid_points[1643].z =   4.3740472570; grid_points[1643].w_fixed =    0.0141738295; grid_points[1643].w_total = 0.0141731261; grid_points[1643].i_atom = 1;
    grid_points[1644].x =  -0.5365666928; grid_points[1644].y =   0.1788555643; grid_points[1644].z =   4.0163361285; grid_points[1644].w_fixed =    0.0141738295; grid_points[1644].w_total = 0.0141666532; grid_points[1644].i_atom = 1;
    grid_points[1645].x =   0.0000000000; grid_points[1645].y =   0.7747876550; grid_points[1645].z =   4.1951916927; grid_points[1645].w_fixed =    0.0059666325; grid_points[1645].w_total = 0.0059644033; grid_points[1645].i_atom = 1;
    grid_points[1646].x =  -0.4473238612; grid_points[1646].y =   0.4473238612; grid_points[1646].z =   4.6425155539; grid_points[1646].w_fixed =    0.0152642382; grid_points[1646].w_total = 0.0152405851; grid_points[1646].i_atom = 1;
    grid_points[1647].x =  -0.4473238612; grid_points[1647].y =   0.4473238612; grid_points[1647].z =   3.7478678316; grid_points[1647].w_fixed =    0.0152642382; grid_points[1647].w_total = 0.0150785295; grid_points[1647].i_atom = 1;
    grid_points[1648].x =  -0.1434253090; grid_points[1648].y =   0.1434253090; grid_points[1648].z =   4.9429578806; grid_points[1648].w_fixed =    0.0127985780; grid_points[1648].w_total = 0.0126709482; grid_points[1648].i_atom = 1;
    grid_points[1649].x =  -0.1434253090; grid_points[1649].y =   0.1434253090; grid_points[1649].z =   3.4474255048; grid_points[1649].w_fixed =    0.0127985780; grid_points[1649].w_total = 0.0120391551; grid_points[1649].i_atom = 1;
    grid_points[1650].x =  -0.1434253090; grid_points[1650].y =   0.7477661879; grid_points[1650].z =   4.3386170017; grid_points[1650].w_fixed =    0.0127985780; grid_points[1650].w_total = 0.0127960888; grid_points[1650].i_atom = 1;
    grid_points[1651].x =  -0.1434253090; grid_points[1651].y =   0.7477661879; grid_points[1651].z =   4.0517663837; grid_points[1651].w_fixed =    0.0127985780; grid_points[1651].w_total = 0.0127814266; grid_points[1651].i_atom = 1;
    grid_points[1652].x =  -0.7477661879; grid_points[1652].y =   0.1434253090; grid_points[1652].z =   4.3386170017; grid_points[1652].w_fixed =    0.0127985780; grid_points[1652].w_total = 0.0127960888; grid_points[1652].i_atom = 1;
    grid_points[1653].x =  -0.7477661879; grid_points[1653].y =   0.1434253090; grid_points[1653].z =   4.0517663837; grid_points[1653].w_fixed =    0.0127985780; grid_points[1653].w_total = 0.0127814266; grid_points[1653].i_atom = 1;
    grid_points[1654].x =  -0.5349297050; grid_points[1654].y =   0.5349297050; grid_points[1654].z =   4.3625127365; grid_points[1654].w_fixed =    0.0154965859; grid_points[1654].w_total = 0.0154933439; grid_points[1654].i_atom = 1;
    grid_points[1655].x =  -0.5349297050; grid_points[1655].y =   0.5349297050; grid_points[1655].z =   4.0278706490; grid_points[1655].w_fixed =    0.0154965859; grid_points[1655].w_total = 0.0154711868; grid_points[1655].i_atom = 1;
    grid_points[1656].x =  -0.5349297050; grid_points[1656].y =   0.1673210437; grid_points[1656].z =   4.7301213978; grid_points[1656].w_fixed =    0.0154965859; grid_points[1656].w_total = 0.0154531108; grid_points[1656].i_atom = 1;
    grid_points[1657].x =  -0.5349297050; grid_points[1657].y =   0.1673210437; grid_points[1657].z =   3.6602619877; grid_points[1657].w_fixed =    0.0154965859; grid_points[1657].w_total = 0.0151815178; grid_points[1657].i_atom = 1;
    grid_points[1658].x =  -0.1673210437; grid_points[1658].y =   0.5349297050; grid_points[1658].z =   4.7301213978; grid_points[1658].w_fixed =    0.0154965859; grid_points[1658].w_total = 0.0154531108; grid_points[1658].i_atom = 1;
    grid_points[1659].x =  -0.1673210437; grid_points[1659].y =   0.5349297050; grid_points[1659].z =   3.6602619877; grid_points[1659].w_fixed =    0.0154965859; grid_points[1659].w_total = 0.0151815178; grid_points[1659].i_atom = 1;
    grid_points[1660].x =  -0.3065753189; grid_points[1660].y =   0.3065753189; grid_points[1660].z =   4.8373124430; grid_points[1660].w_fixed =    0.0149552262; grid_points[1660].w_total = 0.0148734574; grid_points[1660].i_atom = 1;
    grid_points[1661].x =  -0.3065753189; grid_points[1661].y =   0.3065753189; grid_points[1661].z =   3.5530709425; grid_points[1661].w_fixed =    0.0149552262; grid_points[1661].w_total = 0.0144183652; grid_points[1661].i_atom = 1;
    grid_points[1662].x =  -0.3065753189; grid_points[1662].y =   0.6421207503; grid_points[1662].z =   4.5017670117; grid_points[1662].w_fixed =    0.0149552262; grid_points[1662].w_total = 0.0149471131; grid_points[1662].i_atom = 1;
    grid_points[1663].x =  -0.3065753189; grid_points[1663].y =   0.6421207503; grid_points[1663].z =   3.8886163738; grid_points[1663].w_fixed =    0.0149552262; grid_points[1663].w_total = 0.0148835083; grid_points[1663].i_atom = 1;
    grid_points[1664].x =  -0.6421207503; grid_points[1664].y =   0.3065753189; grid_points[1664].z =   4.5017670117; grid_points[1664].w_fixed =    0.0149552262; grid_points[1664].w_total = 0.0149471131; grid_points[1664].i_atom = 1;
    grid_points[1665].x =  -0.6421207503; grid_points[1665].y =   0.3065753189; grid_points[1665].z =   3.8886163738; grid_points[1665].w_fixed =    0.0149552262; grid_points[1665].w_total = 0.0148835083; grid_points[1665].i_atom = 1;
    grid_points[1666].x =  -0.3706344180; grid_points[1666].y =   0.6803866830; grid_points[1666].z =   4.1951916927; grid_points[1666].w_fixed =    0.0151103430; grid_points[1666].w_total = 0.0151046974; grid_points[1666].i_atom = 1;
    grid_points[1667].x =  -0.6803866830; grid_points[1667].y =   0.3706344180; grid_points[1667].z =   4.1951916927; grid_points[1667].w_fixed =    0.0151103430; grid_points[1667].w_total = 0.0151046974; grid_points[1667].i_atom = 1;
    grid_points[1668].x =   0.0000000000; grid_points[1668].y =   0.3706344180; grid_points[1668].z =   4.8755783758; grid_points[1668].w_fixed =    0.0151103430; grid_points[1668].w_total = 0.0150070014; grid_points[1668].i_atom = 1;
    grid_points[1669].x =   0.0000000000; grid_points[1669].y =   0.3706344180; grid_points[1669].z =   3.5148050097; grid_points[1669].w_fixed =    0.0151103430; grid_points[1669].w_total = 0.0144554995; grid_points[1669].i_atom = 1;
    grid_points[1670].x =   0.0000000000; grid_points[1670].y =   0.6803866830; grid_points[1670].z =   4.5658261108; grid_points[1670].w_fixed =    0.0151103430; grid_points[1670].w_total = 0.0150969464; grid_points[1670].i_atom = 1;
    grid_points[1671].x =   0.0000000000; grid_points[1671].y =   0.6803866830; grid_points[1671].z =   3.8245572747; grid_points[1671].w_fixed =    0.0151103430; grid_points[1671].w_total = 0.0149977013; grid_points[1671].i_atom = 1;
    grid_points[1672].x =   0.0000000000; grid_points[1672].y =   1.0119675493; grid_points[1672].z =   4.1951916927; grid_points[1672].w_fixed =    0.0133541260; grid_points[1672].w_total = 0.0133291503; grid_points[1672].i_atom = 1;
    grid_points[1673].x =  -0.5842597370; grid_points[1673].y =   0.5842597370; grid_points[1673].z =   4.7794514298; grid_points[1673].w_fixed =    0.0341634180; grid_points[1673].w_total = 0.0338660424; grid_points[1673].i_atom = 1;
    grid_points[1674].x =  -0.5842597370; grid_points[1674].y =   0.5842597370; grid_points[1674].z =   3.6109319557; grid_points[1674].w_fixed =    0.0341634180; grid_points[1674].w_total = 0.0323366953; grid_points[1674].i_atom = 1;
    grid_points[1675].x =  -0.1873310158; grid_points[1675].y =   0.1873310158; grid_points[1675].z =   5.1718658973; grid_points[1675].w_fixed =    0.0286449390; grid_points[1675].w_total = 0.0271394907; grid_points[1675].i_atom = 1;
    grid_points[1676].x =  -0.1873310158; grid_points[1676].y =   0.1873310158; grid_points[1676].z =   3.2185174881; grid_points[1676].w_fixed =    0.0286449390; grid_points[1676].w_total = 0.0223716243; grid_points[1676].i_atom = 1;
    grid_points[1677].x =  -0.1873310158; grid_points[1677].y =   0.9766742046; grid_points[1677].z =   4.3825227086; grid_points[1677].w_fixed =    0.0286449390; grid_points[1677].w_total = 0.0286147757; grid_points[1677].i_atom = 1;
    grid_points[1678].x =  -0.1873310158; grid_points[1678].y =   0.9766742046; grid_points[1678].z =   4.0078606769; grid_points[1678].w_fixed =    0.0286449390; grid_points[1678].w_total = 0.0284587127; grid_points[1678].i_atom = 1;
    grid_points[1679].x =  -0.9766742046; grid_points[1679].y =   0.1873310158; grid_points[1679].z =   4.3825227086; grid_points[1679].w_fixed =    0.0286449390; grid_points[1679].w_total = 0.0286147757; grid_points[1679].i_atom = 1;
    grid_points[1680].x =  -0.9766742046; grid_points[1680].y =   0.1873310158; grid_points[1680].z =   4.0078606769; grid_points[1680].w_fixed =    0.0286449390; grid_points[1680].w_total = 0.0284587127; grid_points[1680].i_atom = 1;
    grid_points[1681].x =  -0.6986836963; grid_points[1681].y =   0.6986836963; grid_points[1681].z =   4.4137334641; grid_points[1681].w_fixed =    0.0346834435; grid_points[1681].w_total = 0.0346436110; grid_points[1681].i_atom = 1;
    grid_points[1682].x =  -0.6986836963; grid_points[1682].y =   0.6986836963; grid_points[1682].z =   3.9766499214; grid_points[1682].w_fixed =    0.0346834435; grid_points[1682].w_total = 0.0344092150; grid_points[1682].i_atom = 1;
    grid_points[1683].x =  -0.6986836963; grid_points[1683].y =   0.2185417714; grid_points[1683].z =   4.8938753891; grid_points[1683].w_fixed =    0.0346834435; grid_points[1683].w_total = 0.0341432423; grid_points[1683].i_atom = 1;
    grid_points[1684].x =  -0.6986836963; grid_points[1684].y =   0.2185417714; grid_points[1684].z =   3.4965079964; grid_points[1684].w_fixed =    0.0346834435; grid_points[1684].w_total = 0.0317117114; grid_points[1684].i_atom = 1;
    grid_points[1685].x =  -0.2185417714; grid_points[1685].y =   0.6986836963; grid_points[1685].z =   4.8938753891; grid_points[1685].w_fixed =    0.0346834435; grid_points[1685].w_total = 0.0341432423; grid_points[1685].i_atom = 1;
    grid_points[1686].x =  -0.2185417714; grid_points[1686].y =   0.6986836963; grid_points[1686].z =   3.4965079964; grid_points[1686].w_fixed =    0.0346834435; grid_points[1686].w_total = 0.0317117114; grid_points[1686].i_atom = 1;
    grid_points[1687].x =  -0.4004249063; grid_points[1687].y =   0.4004249063; grid_points[1687].z =   5.0338800196; grid_points[1687].w_fixed =    0.0334718077; grid_points[1687].w_total = 0.0324778009; grid_points[1687].i_atom = 1;
    grid_points[1688].x =  -0.4004249063; grid_points[1688].y =   0.4004249063; grid_points[1688].z =   3.3565033658; grid_points[1688].w_fixed =    0.0334718077; grid_points[1688].w_total = 0.0287086336; grid_points[1688].i_atom = 1;
    grid_points[1689].x =  -0.4004249063; grid_points[1689].y =   0.8386883269; grid_points[1689].z =   4.5956165991; grid_points[1689].w_fixed =    0.0334718077; grid_points[1689].w_total = 0.0333691081; grid_points[1689].i_atom = 1;
    grid_points[1690].x =  -0.4004249063; grid_points[1690].y =   0.8386883269; grid_points[1690].z =   3.7947667864; grid_points[1690].w_fixed =    0.0334718077; grid_points[1690].w_total = 0.0327270168; grid_points[1690].i_atom = 1;
    grid_points[1691].x =  -0.8386883269; grid_points[1691].y =   0.4004249063; grid_points[1691].z =   4.5956165991; grid_points[1691].w_fixed =    0.0334718077; grid_points[1691].w_total = 0.0333691081; grid_points[1691].i_atom = 1;
    grid_points[1692].x =  -0.8386883269; grid_points[1692].y =   0.4004249063; grid_points[1692].z =   3.7947667864; grid_points[1692].w_fixed =    0.0334718077; grid_points[1692].w_total = 0.0327270168; grid_points[1692].i_atom = 1;
    grid_points[1693].x =  -0.4840939338; grid_points[1693].y =   0.8886683207; grid_points[1693].z =   4.1951916927; grid_points[1693].w_fixed =    0.0338189800; grid_points[1693].w_total = 0.0337557378; grid_points[1693].i_atom = 1;
    grid_points[1694].x =  -0.8886683207; grid_points[1694].y =   0.4840939338; grid_points[1694].z =   4.1951916927; grid_points[1694].w_fixed =    0.0338189800; grid_points[1694].w_total = 0.0337557378; grid_points[1694].i_atom = 1;
    grid_points[1695].x =   0.0000000000; grid_points[1695].y =   0.4840939338; grid_points[1695].z =   5.0838600134; grid_points[1695].w_fixed =    0.0338189800; grid_points[1695].w_total = 0.0325750639; grid_points[1695].i_atom = 1;
    grid_points[1696].x =   0.0000000000; grid_points[1696].y =   0.4840939338; grid_points[1696].z =   3.3065233720; grid_points[1696].w_fixed =    0.0338189800; grid_points[1696].w_total = 0.0281497735; grid_points[1696].i_atom = 1;
    grid_points[1697].x =   0.0000000000; grid_points[1697].y =   0.8886683207; grid_points[1697].z =   4.6792856265; grid_points[1697].w_fixed =    0.0338189800; grid_points[1697].w_total = 0.0336495489; grid_points[1697].i_atom = 1;
    grid_points[1698].x =   0.0000000000; grid_points[1698].y =   0.8886683207; grid_points[1698].z =   3.7110977590; grid_points[1698].w_fixed =    0.0338189800; grid_points[1698].w_total = 0.0326753125; grid_points[1698].i_atom = 1;
    grid_points[1699].x =   0.0000000000; grid_points[1699].y =   1.3249327354; grid_points[1699].z =   4.1951916927; grid_points[1699].w_fixed =    0.0303775643; grid_points[1699].w_total = 0.0301307226; grid_points[1699].i_atom = 1;
    grid_points[1700].x =  -0.7649502715; grid_points[1700].y =   0.7649502715; grid_points[1700].z =   4.9601419642; grid_points[1700].w_fixed =    0.0777139157; grid_points[1700].w_total = 0.0744471256; grid_points[1700].i_atom = 1;
    grid_points[1701].x =  -0.7649502715; grid_points[1701].y =   0.7649502715; grid_points[1701].z =   3.4302414213; grid_points[1701].w_fixed =    0.0777139157; grid_points[1701].w_total = 0.0632755179; grid_points[1701].i_atom = 1;
    grid_points[1702].x =  -0.2452657651; grid_points[1702].y =   0.2452657651; grid_points[1702].z =   5.4739161206; grid_points[1702].w_fixed =    0.0651606456; grid_points[1702].w_total = 0.0509705226; grid_points[1702].i_atom = 1;
    grid_points[1703].x =  -0.2452657651; grid_points[1703].y =   1.2787244278; grid_points[1703].z =   4.4404574578; grid_points[1703].w_fixed =    0.0651606456; grid_points[1703].w_total = 0.0648319178; grid_points[1703].i_atom = 1;
    grid_points[1704].x =  -0.2452657651; grid_points[1704].y =   1.2787244278; grid_points[1704].z =   3.9499259276; grid_points[1704].w_fixed =    0.0651606456; grid_points[1704].w_total = 0.0634126219; grid_points[1704].i_atom = 1;
    grid_points[1705].x =  -0.9147614482; grid_points[1705].y =   0.9147614482; grid_points[1705].z =   4.4813205782; grid_points[1705].w_fixed =    0.0788968540; grid_points[1705].w_total = 0.0784571533; grid_points[1705].i_atom = 1;
    grid_points[1706].x =  -0.9147614482; grid_points[1706].y =   0.9147614482; grid_points[1706].z =   3.9090628073; grid_points[1706].w_fixed =    0.0788968540; grid_points[1706].w_total = 0.0763502659; grid_points[1706].i_atom = 1;
    grid_points[1707].x =  -0.9147614482; grid_points[1707].y =   0.2861288854; grid_points[1707].z =   5.1099531410; grid_points[1707].w_fixed =    0.0788968540; grid_points[1707].w_total = 0.0731312984; grid_points[1707].i_atom = 1;
    grid_points[1708].x =  -0.9147614482; grid_points[1708].y =   0.2861288854; grid_points[1708].z =   3.2804302445; grid_points[1708].w_fixed =    0.0788968540; grid_points[1708].w_total = 0.0571063680; grid_points[1708].i_atom = 1;
    grid_points[1709].x =  -0.2861288854; grid_points[1709].y =   0.9147614482; grid_points[1709].z =   5.1099531410; grid_points[1709].w_fixed =    0.0788968540; grid_points[1709].w_total = 0.0731312984; grid_points[1709].i_atom = 1;
    grid_points[1710].x =  -0.2861288854; grid_points[1710].y =   0.9147614482; grid_points[1710].z =   3.2804302445; grid_points[1710].w_fixed =    0.0788968540; grid_points[1710].w_total = 0.0571063680; grid_points[1710].i_atom = 1;
    grid_points[1711].x =  -0.5242619359; grid_points[1711].y =   0.5242619359; grid_points[1711].z =   5.2932561710; grid_points[1711].w_fixed =    0.0761406613; grid_points[1711].w_total = 0.0660727372; grid_points[1711].i_atom = 1;
    grid_points[1712].x =  -0.5242619359; grid_points[1712].y =   1.0980644783; grid_points[1712].z =   4.7194536287; grid_points[1712].w_fixed =    0.0761406613; grid_points[1712].w_total = 0.0749832017; grid_points[1712].i_atom = 1;
    grid_points[1713].x =  -0.5242619359; grid_points[1713].y =   1.0980644783; grid_points[1713].z =   3.6709297568; grid_points[1713].w_fixed =    0.0761406613; grid_points[1713].w_total = 0.0696636336; grid_points[1713].i_atom = 1;
    grid_points[1714].x =  -1.0980644783; grid_points[1714].y =   0.5242619359; grid_points[1714].z =   4.7194536287; grid_points[1714].w_fixed =    0.0761406613; grid_points[1714].w_total = 0.0749832017; grid_points[1714].i_atom = 1;
    grid_points[1715].x =  -1.0980644783; grid_points[1715].y =   0.5242619359; grid_points[1715].z =   3.6709297568; grid_points[1715].w_fixed =    0.0761406613; grid_points[1715].w_total = 0.0696636336; grid_points[1715].i_atom = 1;
    grid_points[1716].x =  -0.6338067859; grid_points[1716].y =   1.1635014876; grid_points[1716].z =   4.1951916927; grid_points[1716].w_fixed =    0.0769303986; grid_points[1716].w_total = 0.0763061913; grid_points[1716].i_atom = 1;
    grid_points[1717].x =  -1.1635014876; grid_points[1717].y =   0.6338067859; grid_points[1717].z =   4.1951916927; grid_points[1717].w_fixed =    0.0769303986; grid_points[1717].w_total = 0.0763061913; grid_points[1717].i_atom = 1;
    grid_points[1718].x =   0.0000000000; grid_points[1718].y =   0.6338067859; grid_points[1718].z =   5.3586931803; grid_points[1718].w_fixed =    0.0769303986; grid_points[1718].w_total = 0.0646229177; grid_points[1718].i_atom = 1;
    grid_points[1719].x =   0.0000000000; grid_points[1719].y =   1.1635014876; grid_points[1719].z =   4.8289984786; grid_points[1719].w_fixed =    0.0769303986; grid_points[1719].w_total = 0.0750352153; grid_points[1719].i_atom = 1;
    grid_points[1720].x =   0.0000000000; grid_points[1720].y =   1.1635014876; grid_points[1720].z =   3.5613849069; grid_points[1720].w_fixed =    0.0769303986; grid_points[1720].w_total = 0.0673622975; grid_points[1720].i_atom = 1;
    grid_points[1721].x =   0.0000000000; grid_points[1721].y =   1.7432722237; grid_points[1721].z =   4.1951916927; grid_points[1721].w_fixed =    0.0707954930; grid_points[1721].w_total = 0.0686117539; grid_points[1721].i_atom = 1;
    grid_points[1722].x =  -1.0064786876; grid_points[1722].y =   1.0064786876; grid_points[1722].z =   5.2016703804; grid_points[1722].w_fixed =    0.1811137632; grid_points[1722].w_total = 0.1504983152; grid_points[1722].i_atom = 1;
    grid_points[1723].x =  -1.0064786876; grid_points[1723].y =   1.0064786876; grid_points[1723].z =   3.1887130051; grid_points[1723].w_fixed =    0.1811137632; grid_points[1723].w_total = 0.0948100401; grid_points[1723].i_atom = 1;
    grid_points[1724].x =  -0.3227069453; grid_points[1724].y =   0.3227069453; grid_points[1724].z =   5.8776656155; grid_points[1724].w_fixed =    0.1518581277; grid_points[1724].w_total = 0.0603864433; grid_points[1724].i_atom = 1;
    grid_points[1725].x =  -0.3227069453; grid_points[1725].y =   1.6824739228; grid_points[1725].z =   4.5178986380; grid_points[1725].w_fixed =    0.1518581277; grid_points[1725].w_total = 0.1485902764; grid_points[1725].i_atom = 1;
    grid_points[1726].x =  -0.3227069453; grid_points[1726].y =   1.6824739228; grid_points[1726].z =   3.8724847475; grid_points[1726].w_fixed =    0.1518581277; grid_points[1726].w_total = 0.1376185851; grid_points[1726].i_atom = 1;
    grid_points[1727].x =  -1.2035918363; grid_points[1727].y =   1.2035918363; grid_points[1727].z =   4.5716640411; grid_points[1727].w_fixed =    0.1838706234; grid_points[1727].w_total = 0.1794730010; grid_points[1727].i_atom = 1;
    grid_points[1728].x =  -1.2035918363; grid_points[1728].y =   1.2035918363; grid_points[1728].z =   3.8187193444; grid_points[1728].w_fixed =    0.1838706234; grid_points[1728].w_total = 0.1635857757; grid_points[1728].i_atom = 1;
    grid_points[1729].x =  -1.2035918363; grid_points[1729].y =   0.3764723484; grid_points[1729].z =   5.3987835290; grid_points[1729].w_fixed =    0.1838706234; grid_points[1729].w_total = 0.1337193515; grid_points[1729].i_atom = 1;
    grid_points[1730].x =  -0.3764723484; grid_points[1730].y =   1.2035918363; grid_points[1730].z =   5.3987835290; grid_points[1730].w_fixed =    0.1838706234; grid_points[1730].w_total = 0.1337193515; grid_points[1730].i_atom = 1;
    grid_points[1731].x =  -0.6897944676; grid_points[1731].y =   0.6897944676; grid_points[1731].z =   5.6399633809; grid_points[1731].w_fixed =    0.1774472638; grid_points[1731].w_total = 0.1007712847; grid_points[1731].i_atom = 1;
    grid_points[1732].x =  -0.6897944676; grid_points[1732].y =   1.4447716881; grid_points[1732].z =   4.8849861603; grid_points[1732].w_fixed =    0.1774472638; grid_points[1732].w_total = 0.1658131153; grid_points[1732].i_atom = 1;
    grid_points[1733].x =  -0.6897944676; grid_points[1733].y =   1.4447716881; grid_points[1733].z =   3.5053972252; grid_points[1733].w_fixed =    0.1774472638; grid_points[1733].w_total = 0.1314016342; grid_points[1733].i_atom = 1;
    grid_points[1734].x =  -0.8339274406; grid_points[1734].y =   1.5308700368; grid_points[1734].z =   4.1951916927; grid_points[1734].w_fixed =    0.1792877615; grid_points[1734].w_total = 0.1738017119; grid_points[1734].i_atom = 1;
    grid_points[1735].x =   0.0000000000; grid_points[1735].y =   0.8339274406; grid_points[1735].z =   5.7260617296; grid_points[1735].w_fixed =    0.1792877615; grid_points[1735].w_total = 0.0907093839; grid_points[1735].i_atom = 1;
    grid_points[1736].x =   0.0000000000; grid_points[1736].y =   1.5308700368; grid_points[1736].z =   5.0291191333; grid_points[1736].w_fixed =    0.1792877615; grid_points[1736].w_total = 0.1606309925; grid_points[1736].i_atom = 1;
    grid_points[1737].x =   0.0000000000; grid_points[1737].y =   1.5308700368; grid_points[1737].z =   3.3612642521; grid_points[1737].w_fixed =    0.1792877615; grid_points[1737].w_total = 0.1158377808; grid_points[1737].i_atom = 1;
    grid_points[1738].x =   0.0000000000; grid_points[1738].y =   2.3115565574; grid_points[1738].z =   4.1951916927; grid_points[1738].w_fixed =    0.1705811414; grid_points[1738].w_total = 0.1510163491; grid_points[1738].i_atom = 1;
    grid_points[1739].x =  -0.4279052608; grid_points[1739].y =   2.2309387920; grid_points[1739].z =   4.6230969535; grid_points[1739].w_fixed =    0.3659008737; grid_points[1739].w_total = 0.3324753846; grid_points[1739].i_atom = 1;
    grid_points[1740].x =  -0.4279052608; grid_points[1740].y =   2.2309387920; grid_points[1740].z =   3.7672864320; grid_points[1740].w_fixed =    0.3659008737; grid_points[1740].w_total = 0.2578929044; grid_points[1740].i_atom = 1;
    grid_points[1741].x =  -0.4991974941; grid_points[1741].y =   1.5959473017; grid_points[1741].z =   5.7911389945; grid_points[1741].w_fixed =    0.4430347113; grid_points[1741].w_total = 0.1354214943; grid_points[1741].i_atom = 1;
    grid_points[1742].x =  -0.9146585961; grid_points[1742].y =   0.9146585961; grid_points[1742].z =   6.1109403775; grid_points[1742].w_fixed =    0.4275576805; grid_points[1742].w_total = 0.0649370027; grid_points[1742].i_atom = 1;
    grid_points[1743].x =  -0.9146585961; grid_points[1743].y =   1.9157486847; grid_points[1743].z =   5.1098502889; grid_points[1743].w_fixed =    0.4275576805; grid_points[1743].w_total = 0.3174678581; grid_points[1743].i_atom = 1;
    grid_points[1744].x =  -0.9146585961; grid_points[1744].y =   1.9157486847; grid_points[1744].z =   3.2805330966; grid_points[1744].w_fixed =    0.4275576805; grid_points[1744].w_total = 0.1674035783; grid_points[1744].i_atom = 1;
    grid_points[1745].x =  -1.1057770654; grid_points[1745].y =   2.0299139882; grid_points[1745].z =   4.1951916927; grid_points[1745].w_fixed =    0.4319923443; grid_points[1745].w_total = 0.3838311380; grid_points[1745].i_atom = 1;
    grid_points[1746].x =   0.0000000000; grid_points[1746].y =   1.1057770654; grid_points[1746].z =   6.2251056810; grid_points[1746].w_fixed =    0.4319923443; grid_points[1746].w_total = 0.0483760601; grid_points[1746].i_atom = 1;
    grid_points[1747].x =   0.0000000000; grid_points[1747].y =   2.0299139882; grid_points[1747].z =   5.3009687581; grid_points[1747].w_fixed =    0.4319923443; grid_points[1747].w_total = 0.2694163787; grid_points[1747].i_atom = 1;
    grid_points[1748].x =  -1.2263012757; grid_points[1748].y =   2.5684830011; grid_points[1748].z =   5.4214929684; grid_points[1748].w_fixed =    1.0767762850; grid_points[1748].w_total = 0.2278386018; grid_points[1748].i_atom = 1;
    grid_points[1749].x =   0.0000000000; grid_points[1749].y =   2.7215467321; grid_points[1749].z =   5.6777293649; grid_points[1749].w_fixed =    1.0879446980; grid_points[1749].w_total = 0.1296664136; grid_points[1749].i_atom = 1;
    grid_points[1750].x =   0.0000000000; grid_points[1750].y =   3.0991506199; grid_points[1750].z =   4.1951916927; grid_points[1750].w_fixed =    0.4295975400; grid_points[1750].w_total = 0.2078856501; grid_points[1750].i_atom = 1;
    grid_points[1751].x =  -0.5737012360; grid_points[1751].y =   2.9910647516; grid_points[1751].z =   4.7688929288; grid_points[1751].w_fixed =    0.9214976165; grid_points[1751].w_total = 0.4792995334; grid_points[1751].i_atom = 1;
    grid_points[1752].x =  -0.5737012360; grid_points[1752].y =   2.9910647516; grid_points[1752].z =   3.6214904567; grid_points[1752].w_fixed =    0.9214976165; grid_points[1752].w_total = 0.1977834187; grid_points[1752].i_atom = 1;
    grid_points[1753].x =  -2.5220049911; grid_points[1753].y =  -7.5660149734; grid_points[1753].z =   6.7171966839; grid_points[1753].w_fixed =   55.2891360098; grid_points[1753].w_total = 0.0000000000; grid_points[1753].i_atom = 1;
    grid_points[1754].x =  -1.7666549613; grid_points[1754].y =   1.7666549613; grid_points[1754].z =  -1.1047731912; grid_points[1754].w_fixed =   17.3848280115; grid_points[1754].w_total = 0.0000000000; grid_points[1754].i_atom = 1;
    grid_points[1755].x =  -2.1397188201; grid_points[1755].y =   0.6692841749; grid_points[1755].z =   6.3349105128; grid_points[1755].w_fixed =    1.1157541830; grid_points[1755].w_total = 0.0256935379; grid_points[1755].i_atom = 1;
    grid_points[1756].x =  -2.4354299108; grid_points[1756].y =   2.4354299108; grid_points[1756].z =   6.6306216035; grid_points[1756].w_fixed =    6.3162875542; grid_points[1756].w_total = 0.0011258350; grid_points[1756].i_atom = 1;
    grid_points[1757].x =  -1.2718617903; grid_points[1757].y =   1.2718617903; grid_points[1757].z =   8.0107770637; grid_points[1757].w_fixed =    6.0406797362; grid_points[1757].w_total = 0.0000010065; grid_points[1757].i_atom = 1;
    grid_points[1758].x =  -0.4279052608; grid_points[1758].y =   0.4279052608; grid_points[1758].z =   6.4261304847; grid_points[1758].w_fixed =    0.3659008737; grid_points[1758].w_total = 0.0221696713; grid_points[1758].i_atom = 1;
    grid_points[1759].x =  -0.6692841749; grid_points[1759].y =   2.1397188201; grid_points[1759].z =   6.3349105128; grid_points[1759].w_fixed =    1.1157541830; grid_points[1759].w_total = 0.0256935379; grid_points[1759].i_atom = 1;
    grid_points[1760].x =  -0.5737012360; grid_points[1760].y =   0.5737012360; grid_points[1760].z =   7.1862564443; grid_points[1760].w_fixed =    0.9214976165; grid_points[1760].w_total = 0.0005834114; grid_points[1760].i_atom = 1;
    grid_points[1761].x =  -1.2263012757; grid_points[1761].y =   1.2263012757; grid_points[1761].z =   6.7636746939; grid_points[1761].w_fixed =    1.0767762850; grid_points[1761].w_total = 0.0058282274; grid_points[1761].i_atom = 1;
    grid_points[1762].x =   0.0000000000; grid_points[1762].y =   1.4825376722; grid_points[1762].z =   6.9167384249; grid_points[1762].w_fixed =    1.0879446980; grid_points[1762].w_total = 0.0030241246; grid_points[1762].i_atom = 1;
    grid_points[1763].x =   0.0000000000; grid_points[1763].y =   2.9827802928; grid_points[1763].z =   7.1779719856; grid_points[1763].w_fixed =    6.7598180997; grid_points[1763].w_total = 0.0001667045; grid_points[1763].i_atom = 1;
    grid_points[1764].x =   7.5660149734; grid_points[1764].y =   2.5220049911; grid_points[1764].z =   6.7171966839; grid_points[1764].w_fixed =   55.2891360098; grid_points[1764].w_total = 0.0000000000; grid_points[1764].i_atom = 1;
    grid_points[1765].x =  -1.7666549613; grid_points[1765].y =  -5.2999648839; grid_points[1765].z =   2.4285367314; grid_points[1765].w_fixed =   17.3848280115; grid_points[1765].w_total = 0.0000000010; grid_points[1765].i_atom = 1;
    grid_points[1766].x =  -1.2718617903; grid_points[1766].y =  -3.8155853710; grid_points[1766].z =   2.9233299024; grid_points[1766].w_fixed =    6.0406797362; grid_points[1766].w_total = 0.0032694313; grid_points[1766].i_atom = 1;
    grid_points[1767].x =  -1.7666549613; grid_points[1767].y =   5.2999648839; grid_points[1767].z =   2.4285367314; grid_points[1767].w_fixed =   17.3848280115; grid_points[1767].w_total = 0.0000000010; grid_points[1767].i_atom = 1;
    grid_points[1768].x =  -1.2718617903; grid_points[1768].y =   3.8155853710; grid_points[1768].z =   2.9233299024; grid_points[1768].w_fixed =    6.0406797362; grid_points[1768].w_total = 0.0032694313; grid_points[1768].i_atom = 1;
    grid_points[1769].x =   0.0000000000; grid_points[1769].y =  -5.9146259784; grid_points[1769].z =  10.1098176712; grid_points[1769].w_fixed =   61.8712659234; grid_points[1769].w_total = 0.0000000000; grid_points[1769].i_atom = 1;
    grid_points[1770].x =   0.0000000000; grid_points[1770].y =   4.1431731363; grid_points[1770].z =   0.0520185564; grid_points[1770].w_fixed =   19.4544786654; grid_points[1770].w_total = 0.0000000000; grid_points[1770].i_atom = 1;
    grid_points[1771].x =   0.0000000000; grid_points[1771].y =   5.9146259784; grid_points[1771].z =  10.1098176712; grid_points[1771].w_fixed =   61.8712659234; grid_points[1771].w_total = 0.0000000000; grid_points[1771].i_atom = 1;
    grid_points[1772].x =  -1.7666549613; grid_points[1772].y =  -5.2999648839; grid_points[1772].z =   5.9618466541; grid_points[1772].w_fixed =   17.3848280115; grid_points[1772].w_total = 0.0000000058; grid_points[1772].i_atom = 1;
    grid_points[1773].x =  -1.2718617903; grid_points[1773].y =  -3.8155853710; grid_points[1773].z =   5.4670534831; grid_points[1773].w_fixed =    6.0406797362; grid_points[1773].w_total = 0.0167747560; grid_points[1773].i_atom = 1;
    grid_points[1774].x =  -1.2718617903; grid_points[1774].y =   3.8155853710; grid_points[1774].z =   5.4670534831; grid_points[1774].w_fixed =    6.0406797362; grid_points[1774].w_total = 0.0167747560; grid_points[1774].i_atom = 1;
    grid_points[1775].x =  -1.7666549613; grid_points[1775].y =   5.2999648839; grid_points[1775].z =   5.9618466541; grid_points[1775].w_fixed =   17.3848280115; grid_points[1775].w_total = 0.0000000058; grid_points[1775].i_atom = 1;
    grid_points[1776].x =   0.0000000000; grid_points[1776].y =   4.2182883437; grid_points[1776].z =   4.1951916927; grid_points[1776].w_fixed =    3.8023976811; grid_points[1776].w_total = 0.0216466137; grid_points[1776].i_atom = 1;
    grid_points[1777].x =   0.0000000000; grid_points[1777].y =   5.8593316407; grid_points[1777].z =   4.1951916927; grid_points[1777].w_fixed =   10.9431442493; grid_points[1777].w_total = 0.0000000054; grid_points[1777].i_atom = 1;
    grid_points[1778].x =   0.0000000000; grid_points[1778].y =  -4.1431731363; grid_points[1778].z =   8.3383648291; grid_points[1778].w_fixed =   19.4544786654; grid_points[1778].w_total = 0.0000000000; grid_points[1778].i_atom = 1;
    grid_points[1779].x =   0.0000000000; grid_points[1779].y =   5.9146259784; grid_points[1779].z =  -1.7194342857; grid_points[1779].w_fixed =   61.8712659234; grid_points[1779].w_total = 0.0000000000; grid_points[1779].i_atom = 1;
    grid_points[1780].x =   0.0000000000; grid_points[1780].y =   4.1431731363; grid_points[1780].z =   8.3383648291; grid_points[1780].w_fixed =   19.4544786654; grid_points[1780].w_total = 0.0000000000; grid_points[1780].i_atom = 1;
    grid_points[1781].x =  -1.7892954446; grid_points[1781].y =  -1.7892954446; grid_points[1781].z =   2.4058962481; grid_points[1781].w_fixed =    1.0990251469; grid_points[1781].w_total = 0.0123810957; grid_points[1781].i_atom = 1;
    grid_points[1782].x =  -2.4354299108; grid_points[1782].y =  -2.4354299108; grid_points[1782].z =   1.7597617820; grid_points[1782].w_fixed =    6.3162875542; grid_points[1782].w_total = 0.0001328651; grid_points[1782].i_atom = 1;
    grid_points[1783].x =  -1.2718617903; grid_points[1783].y =  -1.2718617903; grid_points[1783].z =   0.3796063218; grid_points[1783].w_fixed =    6.0406797362; grid_points[1783].w_total = 0.0000000021; grid_points[1783].i_atom = 1;
    grid_points[1784].x =  -1.3345778006; grid_points[1784].y =  -1.3345778006; grid_points[1784].z =   2.8606138921; grid_points[1784].w_fixed =    0.4363920800; grid_points[1784].w_total = 0.0730945598; grid_points[1784].i_atom = 1;
    grid_points[1785].x =  -1.7666549613; grid_points[1785].y =  -1.7666549613; grid_points[1785].z =   9.4951565767; grid_points[1785].w_fixed =   17.3848280115; grid_points[1785].w_total = 0.0000000000; grid_points[1785].i_atom = 1;
    grid_points[1786].x =  -2.5220049911; grid_points[1786].y =  -2.5220049911; grid_points[1786].z =  11.7612066662; grid_points[1786].w_fixed =   55.2891360098; grid_points[1786].w_total = 0.0000000000; grid_points[1786].i_atom = 1;
    grid_points[1787].x =  -1.5959473017; grid_points[1787].y =  -0.4991974941; grid_points[1787].z =   2.5992443910; grid_points[1787].w_fixed =    0.4430347113; grid_points[1787].w_total = 0.0354740456; grid_points[1787].i_atom = 1;
    grid_points[1788].x =  -2.1397188201; grid_points[1788].y =  -0.6692841749; grid_points[1788].z =   2.0554728727; grid_points[1788].w_fixed =    1.1157541830; grid_points[1788].w_total = 0.0031940004; grid_points[1788].i_atom = 1;
    grid_points[1789].x =  -1.4825376722; grid_points[1789].y =   0.0000000000; grid_points[1789].z =   1.4736449606; grid_points[1789].w_fixed =    1.0879446980; grid_points[1789].w_total = 0.0000849175; grid_points[1789].i_atom = 1;
    grid_points[1790].x =  -2.7215467321; grid_points[1790].y =   0.0000000000; grid_points[1790].z =   2.7126540206; grid_points[1790].w_fixed =    1.0879446980; grid_points[1790].w_total = 0.0291805697; grid_points[1790].i_atom = 1;
    grid_points[1791].x =  -2.9827802928; grid_points[1791].y =   0.0000000000; grid_points[1791].z =   1.2124113999; grid_points[1791].w_fixed =    6.7598180997; grid_points[1791].w_total = 0.0000096089; grid_points[1791].i_atom = 1;
    grid_points[1792].x =  -2.0299139882; grid_points[1792].y =   0.0000000000; grid_points[1792].z =   3.0894146273; grid_points[1792].w_fixed =    0.4319923443; grid_points[1792].w_total = 0.1180180438; grid_points[1792].i_atom = 1;
    grid_points[1793].x =  -2.5684830011; grid_points[1793].y =  -1.2263012757; grid_points[1793].z =   2.9688904170; grid_points[1793].w_fixed =    1.0767762850; grid_points[1793].w_total = 0.0592143299; grid_points[1793].i_atom = 1;
    grid_points[1794].x =  -2.5220049911; grid_points[1794].y =   7.5660149734; grid_points[1794].z =   1.6731867016; grid_points[1794].w_fixed =   55.2891360098; grid_points[1794].w_total = 0.0000000000; grid_points[1794].i_atom = 1;
    grid_points[1795].x =  -0.4991974941; grid_points[1795].y =  -1.5959473017; grid_points[1795].z =   2.5992443910; grid_points[1795].w_fixed =    0.4430347113; grid_points[1795].w_total = 0.0354740456; grid_points[1795].i_atom = 1;
    grid_points[1796].x =  -0.6692841749; grid_points[1796].y =  -2.1397188201; grid_points[1796].z =   2.0554728727; grid_points[1796].w_fixed =    1.1157541830; grid_points[1796].w_total = 0.0031940004; grid_points[1796].i_atom = 1;
    grid_points[1797].x =   0.0000000000; grid_points[1797].y =  -2.0299139882; grid_points[1797].z =   3.0894146273; grid_points[1797].w_fixed =    0.4319923443; grid_points[1797].w_total = 0.1180180438; grid_points[1797].i_atom = 1;
    grid_points[1798].x =  -1.2263012757; grid_points[1798].y =  -2.5684830011; grid_points[1798].z =   2.9688904170; grid_points[1798].w_fixed =    1.0767762850; grid_points[1798].w_total = 0.0592143299; grid_points[1798].i_atom = 1;
    grid_points[1799].x =   0.0000000000; grid_points[1799].y =   0.0000000000; grid_points[1799].z =   2.4519194691; grid_points[1799].w_fixed =    0.0707954930; grid_points[1799].w_total = 0.0059395634; grid_points[1799].i_atom = 1;
    grid_points[1800].x =  -0.3227069453; grid_points[1800].y =  -0.3227069453; grid_points[1800].z =   2.5127177700; grid_points[1800].w_fixed =    0.1518581277; grid_points[1800].w_total = 0.0160452730; grid_points[1800].i_atom = 1;
    grid_points[1801].x =  -0.8339274406; grid_points[1801].y =   0.0000000000; grid_points[1801].z =   2.6643216559; grid_points[1801].w_fixed =    0.1792877615; grid_points[1801].w_total = 0.0310706777; grid_points[1801].i_atom = 1;
    grid_points[1802].x =   0.0000000000; grid_points[1802].y =  -0.8339274406; grid_points[1802].z =   2.6643216559; grid_points[1802].w_fixed =    0.1792877615; grid_points[1802].w_total = 0.0310706777; grid_points[1802].i_atom = 1;
    grid_points[1803].x =   0.0000000000; grid_points[1803].y =   0.0000000000; grid_points[1803].z =   1.8836351354; grid_points[1803].w_fixed =    0.1705811414; grid_points[1803].w_total = 0.0003951358; grid_points[1803].i_atom = 1;
    grid_points[1804].x =  -0.4279052608; grid_points[1804].y =  -0.4279052608; grid_points[1804].z =   1.9642529008; grid_points[1804].w_fixed =    0.3659008737; grid_points[1804].w_total = 0.0015586987; grid_points[1804].i_atom = 1;
    grid_points[1805].x =  -0.9146585961; grid_points[1805].y =  -0.9146585961; grid_points[1805].z =   2.2794430080; grid_points[1805].w_fixed =    0.4275576805; grid_points[1805].w_total = 0.0104656786; grid_points[1805].i_atom = 1;
    grid_points[1806].x =  -1.1057770654; grid_points[1806].y =   0.0000000000; grid_points[1806].z =   2.1652777045; grid_points[1806].w_fixed =    0.4319923443; grid_points[1806].w_total = 0.0061142765; grid_points[1806].i_atom = 1;
    grid_points[1807].x =   0.0000000000; grid_points[1807].y =  -1.1057770654; grid_points[1807].z =   2.1652777045; grid_points[1807].w_fixed =    0.4319923443; grid_points[1807].w_total = 0.0061142765; grid_points[1807].i_atom = 1;
    grid_points[1808].x =   0.0000000000; grid_points[1808].y =   0.0000000000; grid_points[1808].z =   1.0960410729; grid_points[1808].w_fixed =    0.4295975400; grid_points[1808].w_total = 0.0000001272; grid_points[1808].i_atom = 1;
    grid_points[1809].x =  -0.5737012360; grid_points[1809].y =  -0.5737012360; grid_points[1809].z =   1.2041269412; grid_points[1809].w_fixed =    0.9214976165; grid_points[1809].w_total = 0.0000026829; grid_points[1809].i_atom = 1;
    grid_points[1810].x =  -1.2263012757; grid_points[1810].y =  -1.2263012757; grid_points[1810].z =   1.6267086916; grid_points[1810].w_fixed =    1.0767762850; grid_points[1810].w_total = 0.0002888020; grid_points[1810].i_atom = 1;
    grid_points[1811].x =   0.0000000000; grid_points[1811].y =   0.0000000000; grid_points[1811].z =   2.8702589573; grid_points[1811].w_fixed =    0.0303775643; grid_points[1811].w_total = 0.0122786614; grid_points[1811].i_atom = 1;
    grid_points[1812].x =  -0.2452657651; grid_points[1812].y =  -0.2452657651; grid_points[1812].z =   2.9164672649; grid_points[1812].w_fixed =    0.0651606456; grid_points[1812].w_total = 0.0288023961; grid_points[1812].i_atom = 1;
    grid_points[1813].x =  -0.5242619359; grid_points[1813].y =  -0.5242619359; grid_points[1813].z =   3.0971272145; grid_points[1813].w_fixed =    0.0761406613; grid_points[1813].w_total = 0.0448544141; grid_points[1813].i_atom = 1;
    grid_points[1814].x =  -0.6338067859; grid_points[1814].y =   0.0000000000; grid_points[1814].z =   3.0316902052; grid_points[1814].w_fixed =    0.0769303986; grid_points[1814].w_total = 0.0412765861; grid_points[1814].i_atom = 1;
    grid_points[1815].x =   0.0000000000; grid_points[1815].y =  -0.6338067859; grid_points[1815].z =   3.0316902052; grid_points[1815].w_fixed =    0.0769303986; grid_points[1815].w_total = 0.0412765861; grid_points[1815].i_atom = 1;
    grid_points[1816].x =  -1.2035918363; grid_points[1816].y =  -0.3764723484; grid_points[1816].z =   2.9915998565; grid_points[1816].w_fixed =    0.1838706234; grid_points[1816].w_total = 0.0694127716; grid_points[1816].i_atom = 1;
    grid_points[1817].x =  -0.3764723484; grid_points[1817].y =  -1.2035918363; grid_points[1817].z =   2.9915998565; grid_points[1817].w_fixed =    0.1838706234; grid_points[1817].w_total = 0.0694127716; grid_points[1817].i_atom = 1;
    grid_points[1818].x =  -0.6897944676; grid_points[1818].y =  -0.6897944676; grid_points[1818].z =   2.7504200046; grid_points[1818].w_fixed =    0.1774472638; grid_points[1818].w_total = 0.0390890156; grid_points[1818].i_atom = 1;
    grid_points[1819].x =   0.0000000000; grid_points[1819].y =   0.0000000000; grid_points[1819].z =  10.0545233334; grid_points[1819].w_fixed =   10.9431442493; grid_points[1819].w_total = 0.0000000000; grid_points[1819].i_atom = 1;
    grid_points[1820].x =   0.0000000000; grid_points[1820].y =   0.0000000000; grid_points[1820].z =  12.5597359678; grid_points[1820].w_fixed =   34.8025870819; grid_points[1820].w_total = 0.0000000000; grid_points[1820].i_atom = 1;
    grid_points[1821].x =   7.5660149734; grid_points[1821].y =  -2.5220049911; grid_points[1821].z =   1.6731867016; grid_points[1821].w_fixed =   55.2891360098; grid_points[1821].w_total = 0.0000000000; grid_points[1821].i_atom = 1;
    grid_points[1822].x =  -2.5220049911; grid_points[1822].y =  -2.5220049911; grid_points[1822].z =  -3.3708232807; grid_points[1822].w_fixed =   55.2891360098; grid_points[1822].w_total = 0.0000000000; grid_points[1822].i_atom = 1;
    grid_points[1823].x =  -1.3345778006; grid_points[1823].y =  -1.3345778006; grid_points[1823].z =   5.5297694934; grid_points[1823].w_fixed =    0.4363920800; grid_points[1823].w_total = 0.2072585535; grid_points[1823].i_atom = 1;
    grid_points[1824].x =  -1.5959473017; grid_points[1824].y =  -1.5959473017; grid_points[1824].z =   4.6943891868; grid_points[1824].w_fixed =    0.4430347113; grid_points[1824].w_total = 0.3987026542; grid_points[1824].i_atom = 1;
    grid_points[1825].x =  -1.5959473017; grid_points[1825].y =  -1.5959473017; grid_points[1825].z =   3.6959941987; grid_points[1825].w_fixed =    0.4430347113; grid_points[1825].w_total = 0.2967672396; grid_points[1825].i_atom = 1;
    grid_points[1826].x =  -1.7892954446; grid_points[1826].y =  -1.7892954446; grid_points[1826].z =   5.9844871374; grid_points[1826].w_fixed =    1.0990251469; grid_points[1826].w_total = 0.0671861197; grid_points[1826].i_atom = 1;
    grid_points[1827].x =  -2.1397188201; grid_points[1827].y =  -2.1397188201; grid_points[1827].z =   4.8644758676; grid_points[1827].w_fixed =    1.1157541830; grid_points[1827].w_total = 0.5622597816; grid_points[1827].i_atom = 1;
    grid_points[1828].x =  -2.1397188201; grid_points[1828].y =  -2.1397188201; grid_points[1828].z =   3.5259075179; grid_points[1828].w_fixed =    1.1157541830; grid_points[1828].w_total = 0.2172616206; grid_points[1828].i_atom = 1;
    grid_points[1829].x =  -1.4825376722; grid_points[1829].y =  -2.7215467321; grid_points[1829].z =   4.1951916927; grid_points[1829].w_fixed =    1.0879446980; grid_points[1829].w_total = 0.5523500480; grid_points[1829].i_atom = 1;
    grid_points[1830].x =  -2.7215467321; grid_points[1830].y =  -1.4825376722; grid_points[1830].z =   4.1951916927; grid_points[1830].w_fixed =    1.0879446980; grid_points[1830].w_total = 0.5523500480; grid_points[1830].i_atom = 1;
    grid_points[1831].x =  -2.9827802928; grid_points[1831].y =  -2.9827802928; grid_points[1831].z =   4.1951916927; grid_points[1831].w_fixed =    6.7598180997; grid_points[1831].w_total = 0.0643528246; grid_points[1831].i_atom = 1;
    grid_points[1832].x =  -1.3249327354; grid_points[1832].y =   0.0000000000; grid_points[1832].z =   4.1951916927; grid_points[1832].w_fixed =    0.0303775643; grid_points[1832].w_total = 0.0301307226; grid_points[1832].i_atom = 1;
    grid_points[1833].x =  -1.2787244278; grid_points[1833].y =  -0.2452657651; grid_points[1833].z =   4.4404574578; grid_points[1833].w_fixed =    0.0651606456; grid_points[1833].w_total = 0.0648319178; grid_points[1833].i_atom = 1;
    grid_points[1834].x =  -1.2787244278; grid_points[1834].y =  -0.2452657651; grid_points[1834].z =   3.9499259276; grid_points[1834].w_fixed =    0.0651606456; grid_points[1834].w_total = 0.0634126219; grid_points[1834].i_atom = 1;
    grid_points[1835].x =  -1.7432722237; grid_points[1835].y =   0.0000000000; grid_points[1835].z =   4.1951916927; grid_points[1835].w_fixed =    0.0707954930; grid_points[1835].w_total = 0.0686117539; grid_points[1835].i_atom = 1;
    grid_points[1836].x =  -1.6824739228; grid_points[1836].y =  -0.3227069453; grid_points[1836].z =   4.5178986380; grid_points[1836].w_fixed =    0.1518581277; grid_points[1836].w_total = 0.1485902764; grid_points[1836].i_atom = 1;
    grid_points[1837].x =  -1.6824739228; grid_points[1837].y =  -0.3227069453; grid_points[1837].z =   3.8724847475; grid_points[1837].w_fixed =    0.1518581277; grid_points[1837].w_total = 0.1376185851; grid_points[1837].i_atom = 1;
    grid_points[1838].x =  -1.4447716881; grid_points[1838].y =  -0.6897944676; grid_points[1838].z =   4.8849861603; grid_points[1838].w_fixed =    0.1774472638; grid_points[1838].w_total = 0.1658131153; grid_points[1838].i_atom = 1;
    grid_points[1839].x =  -1.4447716881; grid_points[1839].y =  -0.6897944676; grid_points[1839].z =   3.5053972252; grid_points[1839].w_fixed =    0.1774472638; grid_points[1839].w_total = 0.1314016342; grid_points[1839].i_atom = 1;
    grid_points[1840].x =  -1.5308700368; grid_points[1840].y =  -0.8339274406; grid_points[1840].z =   4.1951916927; grid_points[1840].w_fixed =    0.1792877615; grid_points[1840].w_total = 0.1738017119; grid_points[1840].i_atom = 1;
    grid_points[1841].x =  -1.5308700368; grid_points[1841].y =   0.0000000000; grid_points[1841].z =   5.0291191333; grid_points[1841].w_fixed =    0.1792877615; grid_points[1841].w_total = 0.1606309925; grid_points[1841].i_atom = 1;
    grid_points[1842].x =  -1.5308700368; grid_points[1842].y =   0.0000000000; grid_points[1842].z =   3.3612642521; grid_points[1842].w_fixed =    0.1792877615; grid_points[1842].w_total = 0.1158377808; grid_points[1842].i_atom = 1;
    grid_points[1843].x =  -2.3115565574; grid_points[1843].y =   0.0000000000; grid_points[1843].z =   4.1951916927; grid_points[1843].w_fixed =    0.1705811414; grid_points[1843].w_total = 0.1510163491; grid_points[1843].i_atom = 1;
    grid_points[1844].x =  -2.2309387920; grid_points[1844].y =  -0.4279052608; grid_points[1844].z =   4.6230969535; grid_points[1844].w_fixed =    0.3659008737; grid_points[1844].w_total = 0.3324753846; grid_points[1844].i_atom = 1;
    grid_points[1845].x =  -2.2309387920; grid_points[1845].y =  -0.4279052608; grid_points[1845].z =   3.7672864320; grid_points[1845].w_fixed =    0.3659008737; grid_points[1845].w_total = 0.2578929044; grid_points[1845].i_atom = 1;
    grid_points[1846].x =  -1.5959473017; grid_points[1846].y =  -0.4991974941; grid_points[1846].z =   5.7911389945; grid_points[1846].w_fixed =    0.4430347113; grid_points[1846].w_total = 0.1354214943; grid_points[1846].i_atom = 1;
    grid_points[1847].x =  -1.9157486847; grid_points[1847].y =  -0.9146585961; grid_points[1847].z =   5.1098502889; grid_points[1847].w_fixed =    0.4275576805; grid_points[1847].w_total = 0.3174678581; grid_points[1847].i_atom = 1;
    grid_points[1848].x =  -1.9157486847; grid_points[1848].y =  -0.9146585961; grid_points[1848].z =   3.2805330966; grid_points[1848].w_fixed =    0.4275576805; grid_points[1848].w_total = 0.1674035783; grid_points[1848].i_atom = 1;
    grid_points[1849].x =  -2.0299139882; grid_points[1849].y =  -1.1057770654; grid_points[1849].z =   4.1951916927; grid_points[1849].w_fixed =    0.4319923443; grid_points[1849].w_total = 0.3838311380; grid_points[1849].i_atom = 1;
    grid_points[1850].x =  -2.0299139882; grid_points[1850].y =   0.0000000000; grid_points[1850].z =   5.3009687581; grid_points[1850].w_fixed =    0.4319923443; grid_points[1850].w_total = 0.2694163787; grid_points[1850].i_atom = 1;
    grid_points[1851].x =  -3.0991506199; grid_points[1851].y =   0.0000000000; grid_points[1851].z =   4.1951916927; grid_points[1851].w_fixed =    0.4295975400; grid_points[1851].w_total = 0.2078856501; grid_points[1851].i_atom = 1;
    grid_points[1852].x =  -2.9910647516; grid_points[1852].y =  -0.5737012360; grid_points[1852].z =   4.7688929288; grid_points[1852].w_fixed =    0.9214976165; grid_points[1852].w_total = 0.4792995334; grid_points[1852].i_atom = 1;
    grid_points[1853].x =  -2.9910647516; grid_points[1853].y =  -0.5737012360; grid_points[1853].z =   3.6214904567; grid_points[1853].w_fixed =    0.9214976165; grid_points[1853].w_total = 0.1977834187; grid_points[1853].i_atom = 1;
    grid_points[1854].x =  -2.5684830011; grid_points[1854].y =  -1.2263012757; grid_points[1854].z =   5.4214929684; grid_points[1854].w_fixed =    1.0767762850; grid_points[1854].w_total = 0.2278386018; grid_points[1854].i_atom = 1;
    grid_points[1855].x =  -2.7215467321; grid_points[1855].y =   0.0000000000; grid_points[1855].z =   5.6777293649; grid_points[1855].w_fixed =    1.0879446980; grid_points[1855].w_total = 0.1296664136; grid_points[1855].i_atom = 1;
    grid_points[1856].x =   0.0000000000; grid_points[1856].y =  -1.3249327354; grid_points[1856].z =   4.1951916927; grid_points[1856].w_fixed =    0.0303775643; grid_points[1856].w_total = 0.0301307226; grid_points[1856].i_atom = 1;
    grid_points[1857].x =  -0.2452657651; grid_points[1857].y =  -1.2787244278; grid_points[1857].z =   4.4404574578; grid_points[1857].w_fixed =    0.0651606456; grid_points[1857].w_total = 0.0648319178; grid_points[1857].i_atom = 1;
    grid_points[1858].x =  -0.2452657651; grid_points[1858].y =  -1.2787244278; grid_points[1858].z =   3.9499259276; grid_points[1858].w_fixed =    0.0651606456; grid_points[1858].w_total = 0.0634126219; grid_points[1858].i_atom = 1;
    grid_points[1859].x =   0.0000000000; grid_points[1859].y =  -1.7432722237; grid_points[1859].z =   4.1951916927; grid_points[1859].w_fixed =    0.0707954930; grid_points[1859].w_total = 0.0686117539; grid_points[1859].i_atom = 1;
    grid_points[1860].x =  -0.3227069453; grid_points[1860].y =  -1.6824739228; grid_points[1860].z =   4.5178986380; grid_points[1860].w_fixed =    0.1518581277; grid_points[1860].w_total = 0.1485902764; grid_points[1860].i_atom = 1;
    grid_points[1861].x =  -0.3227069453; grid_points[1861].y =  -1.6824739228; grid_points[1861].z =   3.8724847475; grid_points[1861].w_fixed =    0.1518581277; grid_points[1861].w_total = 0.1376185851; grid_points[1861].i_atom = 1;
    grid_points[1862].x =  -0.6897944676; grid_points[1862].y =  -1.4447716881; grid_points[1862].z =   4.8849861603; grid_points[1862].w_fixed =    0.1774472638; grid_points[1862].w_total = 0.1658131153; grid_points[1862].i_atom = 1;
    grid_points[1863].x =  -0.6897944676; grid_points[1863].y =  -1.4447716881; grid_points[1863].z =   3.5053972252; grid_points[1863].w_fixed =    0.1774472638; grid_points[1863].w_total = 0.1314016342; grid_points[1863].i_atom = 1;
    grid_points[1864].x =  -0.8339274406; grid_points[1864].y =  -1.5308700368; grid_points[1864].z =   4.1951916927; grid_points[1864].w_fixed =    0.1792877615; grid_points[1864].w_total = 0.1738017119; grid_points[1864].i_atom = 1;
    grid_points[1865].x =   0.0000000000; grid_points[1865].y =  -1.5308700368; grid_points[1865].z =   5.0291191333; grid_points[1865].w_fixed =    0.1792877615; grid_points[1865].w_total = 0.1606309925; grid_points[1865].i_atom = 1;
    grid_points[1866].x =   0.0000000000; grid_points[1866].y =  -1.5308700368; grid_points[1866].z =   3.3612642521; grid_points[1866].w_fixed =    0.1792877615; grid_points[1866].w_total = 0.1158377808; grid_points[1866].i_atom = 1;
    grid_points[1867].x =   0.0000000000; grid_points[1867].y =  -2.3115565574; grid_points[1867].z =   4.1951916927; grid_points[1867].w_fixed =    0.1705811414; grid_points[1867].w_total = 0.1510163491; grid_points[1867].i_atom = 1;
    grid_points[1868].x =  -0.4279052608; grid_points[1868].y =  -2.2309387920; grid_points[1868].z =   4.6230969535; grid_points[1868].w_fixed =    0.3659008737; grid_points[1868].w_total = 0.3324753846; grid_points[1868].i_atom = 1;
    grid_points[1869].x =  -0.4279052608; grid_points[1869].y =  -2.2309387920; grid_points[1869].z =   3.7672864320; grid_points[1869].w_fixed =    0.3659008737; grid_points[1869].w_total = 0.2578929044; grid_points[1869].i_atom = 1;
    grid_points[1870].x =  -0.4991974941; grid_points[1870].y =  -1.5959473017; grid_points[1870].z =   5.7911389945; grid_points[1870].w_fixed =    0.4430347113; grid_points[1870].w_total = 0.1354214943; grid_points[1870].i_atom = 1;
    grid_points[1871].x =  -0.9146585961; grid_points[1871].y =  -1.9157486847; grid_points[1871].z =   5.1098502889; grid_points[1871].w_fixed =    0.4275576805; grid_points[1871].w_total = 0.3174678581; grid_points[1871].i_atom = 1;
    grid_points[1872].x =  -0.9146585961; grid_points[1872].y =  -1.9157486847; grid_points[1872].z =   3.2805330966; grid_points[1872].w_fixed =    0.4275576805; grid_points[1872].w_total = 0.1674035783; grid_points[1872].i_atom = 1;
    grid_points[1873].x =  -1.1057770654; grid_points[1873].y =  -2.0299139882; grid_points[1873].z =   4.1951916927; grid_points[1873].w_fixed =    0.4319923443; grid_points[1873].w_total = 0.3838311380; grid_points[1873].i_atom = 1;
    grid_points[1874].x =   0.0000000000; grid_points[1874].y =  -2.0299139882; grid_points[1874].z =   5.3009687581; grid_points[1874].w_fixed =    0.4319923443; grid_points[1874].w_total = 0.2694163787; grid_points[1874].i_atom = 1;
    grid_points[1875].x =   0.0000000000; grid_points[1875].y =  -3.0991506199; grid_points[1875].z =   4.1951916927; grid_points[1875].w_fixed =    0.4295975400; grid_points[1875].w_total = 0.2078856501; grid_points[1875].i_atom = 1;
    grid_points[1876].x =  -0.5737012360; grid_points[1876].y =  -2.9910647516; grid_points[1876].z =   4.7688929288; grid_points[1876].w_fixed =    0.9214976165; grid_points[1876].w_total = 0.4792995334; grid_points[1876].i_atom = 1;
    grid_points[1877].x =  -0.5737012360; grid_points[1877].y =  -2.9910647516; grid_points[1877].z =   3.6214904567; grid_points[1877].w_fixed =    0.9214976165; grid_points[1877].w_total = 0.1977834187; grid_points[1877].i_atom = 1;
    grid_points[1878].x =  -1.2263012757; grid_points[1878].y =  -2.5684830011; grid_points[1878].z =   5.4214929684; grid_points[1878].w_fixed =    1.0767762850; grid_points[1878].w_total = 0.2278386018; grid_points[1878].i_atom = 1;
    grid_points[1879].x =   0.0000000000; grid_points[1879].y =  -2.7215467321; grid_points[1879].z =   5.6777293649; grid_points[1879].w_fixed =    1.0879446980; grid_points[1879].w_total = 0.1296664136; grid_points[1879].i_atom = 1;
    grid_points[1880].x =  -0.0009212695; grid_points[1880].y =   0.0000000000; grid_points[1880].z =   4.1951916927; grid_points[1880].w_fixed =    0.0000000034; grid_points[1880].w_total = 0.0000000034; grid_points[1880].i_atom = 1;
    grid_points[1881].x =   0.0000000000; grid_points[1881].y =  -0.0009212695; grid_points[1881].z =   4.1951916927; grid_points[1881].w_fixed =    0.0000000034; grid_points[1881].w_total = 0.0000000034; grid_points[1881].i_atom = 1;
    grid_points[1882].x =   0.0000000000; grid_points[1882].y =   0.0000000000; grid_points[1882].z =   4.1961129622; grid_points[1882].w_fixed =    0.0000000034; grid_points[1882].w_total = 0.0000000034; grid_points[1882].i_atom = 1;
    grid_points[1883].x =   0.0000000000; grid_points[1883].y =   0.0000000000; grid_points[1883].z =   4.1942704232; grid_points[1883].w_fixed =    0.0000000034; grid_points[1883].w_total = 0.0000000034; grid_points[1883].i_atom = 1;
    grid_points[1884].x =  -0.0039529982; grid_points[1884].y =   0.0000000000; grid_points[1884].z =   4.1951916927; grid_points[1884].w_fixed =    0.0000001386; grid_points[1884].w_total = 0.0000001386; grid_points[1884].i_atom = 1;
    grid_points[1885].x =   0.0000000000; grid_points[1885].y =  -0.0039529982; grid_points[1885].z =   4.1951916927; grid_points[1885].w_fixed =    0.0000001386; grid_points[1885].w_total = 0.0000001386; grid_points[1885].i_atom = 1;
    grid_points[1886].x =   0.0000000000; grid_points[1886].y =   0.0000000000; grid_points[1886].z =   4.1991446910; grid_points[1886].w_fixed =    0.0000001386; grid_points[1886].w_total = 0.0000001386; grid_points[1886].i_atom = 1;
    grid_points[1887].x =   0.0000000000; grid_points[1887].y =   0.0000000000; grid_points[1887].z =   4.1912386945; grid_points[1887].w_fixed =    0.0000001386; grid_points[1887].w_total = 0.0000001386; grid_points[1887].i_atom = 1;
    grid_points[1888].x =  -0.0095652797; grid_points[1888].y =   0.0000000000; grid_points[1888].z =   4.1951916927; grid_points[1888].w_fixed =    0.0000013577; grid_points[1888].w_total = 0.0000013577; grid_points[1888].i_atom = 1;
    grid_points[1889].x =   0.0000000000; grid_points[1889].y =  -0.0095652797; grid_points[1889].z =   4.1951916927; grid_points[1889].w_fixed =    0.0000013577; grid_points[1889].w_total = 0.0000013577; grid_points[1889].i_atom = 1;
    grid_points[1890].x =   0.0000000000; grid_points[1890].y =   0.0000000000; grid_points[1890].z =   4.2047569724; grid_points[1890].w_fixed =    0.0000013577; grid_points[1890].w_total = 0.0000013577; grid_points[1890].i_atom = 1;
    grid_points[1891].x =   0.0000000000; grid_points[1891].y =   0.0000000000; grid_points[1891].z =   4.1856264131; grid_points[1891].w_fixed =    0.0000013577; grid_points[1891].w_total = 0.0000013577; grid_points[1891].i_atom = 1;
    grid_points[1892].x =  -0.0183381693; grid_points[1892].y =   0.0000000000; grid_points[1892].z =   4.1951916927; grid_points[1892].w_fixed =    0.0000074515; grid_points[1892].w_total = 0.0000074515; grid_points[1892].i_atom = 1;
    grid_points[1893].x =   0.0000000000; grid_points[1893].y =  -0.0183381693; grid_points[1893].z =   4.1951916927; grid_points[1893].w_fixed =    0.0000074515; grid_points[1893].w_total = 0.0000074515; grid_points[1893].i_atom = 1;
    grid_points[1894].x =   0.0000000000; grid_points[1894].y =   0.0000000000; grid_points[1894].z =   4.2135298621; grid_points[1894].w_fixed =    0.0000074515; grid_points[1894].w_total = 0.0000074515; grid_points[1894].i_atom = 1;
    grid_points[1895].x =   0.0000000000; grid_points[1895].y =   0.0000000000; grid_points[1895].z =   4.1768535234; grid_points[1895].w_fixed =    0.0000074515; grid_points[1895].w_total = 0.0000074515; grid_points[1895].i_atom = 1;
    grid_points[1896].x =  -0.0309915062; grid_points[1896].y =   0.0000000000; grid_points[1896].z =   4.1951916927; grid_points[1896].w_fixed =    0.0000299246; grid_points[1896].w_total = 0.0000299246; grid_points[1896].i_atom = 1;
    grid_points[1897].x =   0.0000000000; grid_points[1897].y =  -0.0309915062; grid_points[1897].z =   4.1951916927; grid_points[1897].w_fixed =    0.0000299246; grid_points[1897].w_total = 0.0000299246; grid_points[1897].i_atom = 1;
    grid_points[1898].x =   0.0000000000; grid_points[1898].y =   0.0000000000; grid_points[1898].z =   4.2261831989; grid_points[1898].w_fixed =    0.0000299246; grid_points[1898].w_total = 0.0000299246; grid_points[1898].i_atom = 1;
    grid_points[1899].x =   0.0000000000; grid_points[1899].y =   0.0000000000; grid_points[1899].z =   4.1642001865; grid_points[1899].w_fixed =    0.0000299246; grid_points[1899].w_total = 0.0000299246; grid_points[1899].i_atom = 1;
    grid_points[1900].x =  -0.0484242284; grid_points[1900].y =   0.0000000000; grid_points[1900].z =   4.1951916927; grid_points[1900].w_fixed =    0.0000990913; grid_points[1900].w_total = 0.0000990913; grid_points[1900].i_atom = 1;
    grid_points[1901].x =   0.0000000000; grid_points[1901].y =  -0.0484242284; grid_points[1901].z =   4.1951916927; grid_points[1901].w_fixed =    0.0000990913; grid_points[1901].w_total = 0.0000990913; grid_points[1901].i_atom = 1;
    grid_points[1902].x =   0.0000000000; grid_points[1902].y =   0.0000000000; grid_points[1902].z =   4.2436159212; grid_points[1902].w_fixed =    0.0000990913; grid_points[1902].w_total = 0.0000990913; grid_points[1902].i_atom = 1;
    grid_points[1903].x =   0.0000000000; grid_points[1903].y =   0.0000000000; grid_points[1903].z =   4.1467674643; grid_points[1903].w_fixed =    0.0000990913; grid_points[1903].w_total = 0.0000990913; grid_points[1903].i_atom = 1;
    grid_points[1904].x =  -0.0717667204; grid_points[1904].y =   0.0000000000; grid_points[1904].z =   4.1951916927; grid_points[1904].w_fixed =    0.0002885050; grid_points[1904].w_total = 0.0002885050; grid_points[1904].i_atom = 1;
    grid_points[1905].x =   0.0000000000; grid_points[1905].y =  -0.0717667204; grid_points[1905].z =   4.1951916927; grid_points[1905].w_fixed =    0.0002885050; grid_points[1905].w_total = 0.0002885050; grid_points[1905].i_atom = 1;
    grid_points[1906].x =   0.0000000000; grid_points[1906].y =   0.0000000000; grid_points[1906].z =   4.2669584131; grid_points[1906].w_fixed =    0.0002885050; grid_points[1906].w_total = 0.0002885050; grid_points[1906].i_atom = 1;
    grid_points[1907].x =   0.0000000000; grid_points[1907].y =   0.0000000000; grid_points[1907].z =   4.1234249723; grid_points[1907].w_fixed =    0.0002885050; grid_points[1907].w_total = 0.0002885050; grid_points[1907].i_atom = 1;
    grid_points[1908].x =  -0.1024512602; grid_points[1908].y =   0.0000000000; grid_points[1908].z =   4.1951916927; grid_points[1908].w_fixed =    0.0007678018; grid_points[1908].w_total = 0.0007678018; grid_points[1908].i_atom = 1;
    grid_points[1909].x =   0.0000000000; grid_points[1909].y =  -0.1024512602; grid_points[1909].z =   4.1951916927; grid_points[1909].w_fixed =    0.0007678018; grid_points[1909].w_total = 0.0007678018; grid_points[1909].i_atom = 1;
    grid_points[1910].x =   0.0000000000; grid_points[1910].y =   0.0000000000; grid_points[1910].z =   4.2976429529; grid_points[1910].w_fixed =    0.0007678018; grid_points[1910].w_total = 0.0007678018; grid_points[1910].i_atom = 1;
    grid_points[1911].x =   0.0000000000; grid_points[1911].y =   0.0000000000; grid_points[1911].z =   4.0927404326; grid_points[1911].w_fixed =    0.0007678018; grid_points[1911].w_total = 0.0007678017; grid_points[1911].i_atom = 1;
    grid_points[1912].x =  -0.1423079366; grid_points[1912].y =   0.0000000000; grid_points[1912].z =   4.1951916927; grid_points[1912].w_fixed =    0.0005474789; grid_points[1912].w_total = 0.0005474789; grid_points[1912].i_atom = 1;
    grid_points[1913].x =   0.0000000000; grid_points[1913].y =  -0.1423079366; grid_points[1913].z =   4.1951916927; grid_points[1913].w_fixed =    0.0005474789; grid_points[1913].w_total = 0.0005474789; grid_points[1913].i_atom = 1;
    grid_points[1914].x =   0.0000000000; grid_points[1914].y =   0.0000000000; grid_points[1914].z =   4.3374996294; grid_points[1914].w_fixed =    0.0005474789; grid_points[1914].w_total = 0.0005474789; grid_points[1914].i_atom = 1;
    grid_points[1915].x =   0.0000000000; grid_points[1915].y =   0.0000000000; grid_points[1915].z =   4.0528837561; grid_points[1915].w_fixed =    0.0005474789; grid_points[1915].w_total = 0.0005474786; grid_points[1915].i_atom = 1;
    grid_points[1916].x =   0.0000000000; grid_points[1916].y =  -0.1006269070; grid_points[1916].z =   4.2958185997; grid_points[1916].w_fixed =    0.0004379831; grid_points[1916].w_total = 0.0004379831; grid_points[1916].i_atom = 1;
    grid_points[1917].x =   0.0000000000; grid_points[1917].y =  -0.1006269070; grid_points[1917].z =   4.0945647857; grid_points[1917].w_fixed =    0.0004379831; grid_points[1917].w_total = 0.0004379831; grid_points[1917].i_atom = 1;
    grid_points[1918].x =  -0.1006269070; grid_points[1918].y =   0.0000000000; grid_points[1918].z =   4.2958185997; grid_points[1918].w_fixed =    0.0004379831; grid_points[1918].w_total = 0.0004379831; grid_points[1918].i_atom = 1;
    grid_points[1919].x =  -0.1006269070; grid_points[1919].y =   0.0000000000; grid_points[1919].z =   4.0945647857; grid_points[1919].w_fixed =    0.0004379831; grid_points[1919].w_total = 0.0004379831; grid_points[1919].i_atom = 1;
    grid_points[1920].x =  -0.1006269070; grid_points[1920].y =  -0.1006269070; grid_points[1920].z =   4.1951916927; grid_points[1920].w_fixed =    0.0004379831; grid_points[1920].w_total = 0.0004379831; grid_points[1920].i_atom = 1;
    grid_points[1921].x =  -0.0821615255; grid_points[1921].y =  -0.0821615255; grid_points[1921].z =   4.2773532183; grid_points[1921].w_fixed =    0.0003695483; grid_points[1921].w_total = 0.0003695483; grid_points[1921].i_atom = 1;
    grid_points[1922].x =  -0.0821615255; grid_points[1922].y =  -0.0821615255; grid_points[1922].z =   4.1130301672; grid_points[1922].w_fixed =    0.0003695483; grid_points[1922].w_total = 0.0003695482; grid_points[1922].i_atom = 1;
    grid_points[1923].x =  -0.1936969137; grid_points[1923].y =   0.0000000000; grid_points[1923].z =   4.1951916927; grid_points[1923].w_fixed =    0.0013046077; grid_points[1923].w_total = 0.0013046076; grid_points[1923].i_atom = 1;
    grid_points[1924].x =   0.0000000000; grid_points[1924].y =  -0.1936969137; grid_points[1924].z =   4.1951916927; grid_points[1924].w_fixed =    0.0013046077; grid_points[1924].w_total = 0.0013046076; grid_points[1924].i_atom = 1;
    grid_points[1925].x =   0.0000000000; grid_points[1925].y =   0.0000000000; grid_points[1925].z =   4.3888886065; grid_points[1925].w_fixed =    0.0013046077; grid_points[1925].w_total = 0.0013046072; grid_points[1925].i_atom = 1;
    grid_points[1926].x =   0.0000000000; grid_points[1926].y =   0.0000000000; grid_points[1926].z =   4.0014947790; grid_points[1926].w_fixed =    0.0013046077; grid_points[1926].w_total = 0.0013046011; grid_points[1926].i_atom = 1;
    grid_points[1927].x =   0.0000000000; grid_points[1927].y =  -0.1369644012; grid_points[1927].z =   4.3321560939; grid_points[1927].w_fixed =    0.0010436861; grid_points[1927].w_total = 0.0010436860; grid_points[1927].i_atom = 1;
    grid_points[1928].x =   0.0000000000; grid_points[1928].y =  -0.1369644012; grid_points[1928].z =   4.0582272915; grid_points[1928].w_fixed =    0.0010436861; grid_points[1928].w_total = 0.0010436847; grid_points[1928].i_atom = 1;
    grid_points[1929].x =  -0.1369644012; grid_points[1929].y =   0.0000000000; grid_points[1929].z =   4.3321560939; grid_points[1929].w_fixed =    0.0010436861; grid_points[1929].w_total = 0.0010436860; grid_points[1929].i_atom = 1;
    grid_points[1930].x =  -0.1369644012; grid_points[1930].y =   0.0000000000; grid_points[1930].z =   4.0582272915; grid_points[1930].w_fixed =    0.0010436861; grid_points[1930].w_total = 0.0010436847; grid_points[1930].i_atom = 1;
    grid_points[1931].x =  -0.1369644012; grid_points[1931].y =  -0.1369644012; grid_points[1931].z =   4.1951916927; grid_points[1931].w_fixed =    0.0010436861; grid_points[1931].w_total = 0.0010436861; grid_points[1931].i_atom = 1;
    grid_points[1932].x =  -0.1118309653; grid_points[1932].y =  -0.1118309653; grid_points[1932].z =   4.3070226580; grid_points[1932].w_fixed =    0.0008806102; grid_points[1932].w_total = 0.0008806101; grid_points[1932].i_atom = 1;
    grid_points[1933].x =  -0.1118309653; grid_points[1933].y =  -0.1118309653; grid_points[1933].z =   4.0833607275; grid_points[1933].w_fixed =    0.0008806102; grid_points[1933].w_total = 0.0008806095; grid_points[1933].i_atom = 1;
    grid_points[1934].x =  -0.2596933691; grid_points[1934].y =   0.0000000000; grid_points[1934].z =   4.1951916927; grid_points[1934].w_fixed =    0.0030086947; grid_points[1934].w_total = 0.0030086941; grid_points[1934].i_atom = 1;
    grid_points[1935].x =   0.0000000000; grid_points[1935].y =  -0.2596933691; grid_points[1935].z =   4.1951916927; grid_points[1935].w_fixed =    0.0030086947; grid_points[1935].w_total = 0.0030086941; grid_points[1935].i_atom = 1;
    grid_points[1936].x =   0.0000000000; grid_points[1936].y =   0.0000000000; grid_points[1936].z =   4.4548850619; grid_points[1936].w_fixed =    0.0030086947; grid_points[1936].w_total = 0.0030086831; grid_points[1936].i_atom = 1;
    grid_points[1937].x =   0.0000000000; grid_points[1937].y =   0.0000000000; grid_points[1937].z =   3.9354983236; grid_points[1937].w_fixed =    0.0030086947; grid_points[1937].w_total = 0.0030085569; grid_points[1937].i_atom = 1;
    grid_points[1938].x =   0.0000000000; grid_points[1938].y =  -0.1836309423; grid_points[1938].z =   4.3788226351; grid_points[1938].w_fixed =    0.0024069557; grid_points[1938].w_total = 0.0024069533; grid_points[1938].i_atom = 1;
    grid_points[1939].x =   0.0000000000; grid_points[1939].y =  -0.1836309423; grid_points[1939].z =   4.0115607504; grid_points[1939].w_fixed =    0.0024069557; grid_points[1939].w_total = 0.0024069249; grid_points[1939].i_atom = 1;
    grid_points[1940].x =  -0.1836309423; grid_points[1940].y =   0.0000000000; grid_points[1940].z =   4.3788226351; grid_points[1940].w_fixed =    0.0024069557; grid_points[1940].w_total = 0.0024069533; grid_points[1940].i_atom = 1;
    grid_points[1941].x =  -0.1836309423; grid_points[1941].y =   0.0000000000; grid_points[1941].z =   4.0115607504; grid_points[1941].w_fixed =    0.0024069557; grid_points[1941].w_total = 0.0024069249; grid_points[1941].i_atom = 1;
    grid_points[1942].x =  -0.1836309423; grid_points[1942].y =  -0.1836309423; grid_points[1942].z =   4.1951916927; grid_points[1942].w_fixed =    0.0024069557; grid_points[1942].w_total = 0.0024069553; grid_points[1942].i_atom = 1;
    grid_points[1943].x =  -0.1499340366; grid_points[1943].y =  -0.1499340366; grid_points[1943].z =   4.3451257293; grid_points[1943].w_fixed =    0.0020308689; grid_points[1943].w_total = 0.0020308678; grid_points[1943].i_atom = 1;
    grid_points[1944].x =  -0.1499340366; grid_points[1944].y =  -0.1499340366; grid_points[1944].z =   4.0452576562; grid_points[1944].w_fixed =    0.0020308689; grid_points[1944].w_total = 0.0020308551; grid_points[1944].i_atom = 1;
    grid_points[1945].x =  -0.3443500689; grid_points[1945].y =   0.0000000000; grid_points[1945].z =   4.1951916927; grid_points[1945].w_fixed =    0.0067871679; grid_points[1945].w_total = 0.0067871580; grid_points[1945].i_atom = 1;
    grid_points[1946].x =   0.0000000000; grid_points[1946].y =  -0.3443500689; grid_points[1946].z =   4.1951916927; grid_points[1946].w_fixed =    0.0067871679; grid_points[1946].w_total = 0.0067871580; grid_points[1946].i_atom = 1;
    grid_points[1947].x =   0.0000000000; grid_points[1947].y =   0.0000000000; grid_points[1947].z =   4.5395417616; grid_points[1947].w_fixed =    0.0067871679; grid_points[1947].w_total = 0.0067869400; grid_points[1947].i_atom = 1;
    grid_points[1948].x =   0.0000000000; grid_points[1948].y =   0.0000000000; grid_points[1948].z =   3.8508416239; grid_points[1948].w_fixed =    0.0067871679; grid_points[1948].w_total = 0.0067847341; grid_points[1948].i_atom = 1;
    grid_points[1949].x =   0.0000000000; grid_points[1949].y =  -0.2434922688; grid_points[1949].z =   4.4386839615; grid_points[1949].w_fixed =    0.0054297343; grid_points[1949].w_total = 0.0054296870; grid_points[1949].i_atom = 1;
    grid_points[1950].x =   0.0000000000; grid_points[1950].y =  -0.2434922688; grid_points[1950].z =   3.9516994239; grid_points[1950].w_fixed =    0.0054297343; grid_points[1950].w_total = 0.0054291838; grid_points[1950].i_atom = 1;
    grid_points[1951].x =  -0.2434922688; grid_points[1951].y =   0.0000000000; grid_points[1951].z =   4.4386839615; grid_points[1951].w_fixed =    0.0054297343; grid_points[1951].w_total = 0.0054296870; grid_points[1951].i_atom = 1;
    grid_points[1952].x =  -0.2434922688; grid_points[1952].y =   0.0000000000; grid_points[1952].z =   3.9516994239; grid_points[1952].w_fixed =    0.0054297343; grid_points[1952].w_total = 0.0054291838; grid_points[1952].i_atom = 1;
    grid_points[1953].x =  -0.2434922688; grid_points[1953].y =  -0.2434922688; grid_points[1953].z =   4.1951916927; grid_points[1953].w_fixed =    0.0054297343; grid_points[1953].w_total = 0.0054297264; grid_points[1953].i_atom = 1;
    grid_points[1954].x =  -0.1988106050; grid_points[1954].y =  -0.1988106050; grid_points[1954].z =   4.3940022977; grid_points[1954].w_fixed =    0.0045813383; grid_points[1954].w_total = 0.0045813179; grid_points[1954].i_atom = 1;
    grid_points[1955].x =  -0.1988106050; grid_points[1955].y =  -0.1988106050; grid_points[1955].z =   3.9963810878; grid_points[1955].w_fixed =    0.0045813383; grid_points[1955].w_total = 0.0045810919; grid_points[1955].i_atom = 1;
    grid_points[1956].x =  -0.4530765180; grid_points[1956].y =   0.0000000000; grid_points[1956].z =   4.1951916927; grid_points[1956].w_fixed =    0.0040293310; grid_points[1956].w_total = 0.0040292893; grid_points[1956].i_atom = 1;
    grid_points[1957].x =   0.0000000000; grid_points[1957].y =  -0.4530765180; grid_points[1957].z =   4.1951916927; grid_points[1957].w_fixed =    0.0040293310; grid_points[1957].w_total = 0.0040292893; grid_points[1957].i_atom = 1;
    grid_points[1958].x =   0.0000000000; grid_points[1958].y =   0.0000000000; grid_points[1958].z =   4.6482682107; grid_points[1958].w_fixed =    0.0040293310; grid_points[1958].w_total = 0.0040282660; grid_points[1958].i_atom = 1;
    grid_points[1959].x =   0.0000000000; grid_points[1959].y =   0.0000000000; grid_points[1959].z =   3.7421151748; grid_points[1959].w_fixed =    0.0040293310; grid_points[1959].w_total = 0.0040195193; grid_points[1959].i_atom = 1;
    grid_points[1960].x =   0.0000000000; grid_points[1960].y =  -0.3203734782; grid_points[1960].z =   4.5155651710; grid_points[1960].w_fixed =    0.0071632550; grid_points[1960].w_total = 0.0071627669; grid_points[1960].i_atom = 1;
    grid_points[1961].x =   0.0000000000; grid_points[1961].y =  -0.3203734782; grid_points[1961].z =   3.8748182145; grid_points[1961].w_fixed =    0.0071632550; grid_points[1961].w_total = 0.0071582150; grid_points[1961].i_atom = 1;
    grid_points[1962].x =  -0.3203734782; grid_points[1962].y =   0.0000000000; grid_points[1962].z =   4.5155651710; grid_points[1962].w_fixed =    0.0071632550; grid_points[1962].w_total = 0.0071627669; grid_points[1962].i_atom = 1;
    grid_points[1963].x =  -0.3203734782; grid_points[1963].y =   0.0000000000; grid_points[1963].z =   3.8748182145; grid_points[1963].w_fixed =    0.0071632550; grid_points[1963].w_total = 0.0071582150; grid_points[1963].i_atom = 1;
    grid_points[1964].x =  -0.3203734782; grid_points[1964].y =  -0.3203734782; grid_points[1964].z =   4.1951916927; grid_points[1964].w_fixed =    0.0071632550; grid_points[1964].w_total = 0.0071631809; grid_points[1964].i_atom = 1;
    grid_points[1965].x =  -0.2615838496; grid_points[1965].y =  -0.2615838496; grid_points[1965].z =   4.4567755423; grid_points[1965].w_fixed =    0.0066932539; grid_points[1965].w_total = 0.0066930222; grid_points[1965].i_atom = 1;
    grid_points[1966].x =  -0.2615838496; grid_points[1966].y =  -0.2615838496; grid_points[1966].z =   3.9336078431; grid_points[1966].w_fixed =    0.0066932539; grid_points[1966].w_total = 0.0066907381; grid_points[1966].i_atom = 1;
    grid_points[1967].x =  -0.1366077101; grid_points[1967].y =  -0.1366077101; grid_points[1967].z =   4.6050148231; grid_points[1967].w_fixed =    0.0064011973; grid_points[1967].w_total = 0.0064000836; grid_points[1967].i_atom = 1;
    grid_points[1968].x =  -0.1366077101; grid_points[1968].y =  -0.1366077101; grid_points[1968].z =   3.7853685624; grid_points[1968].w_fixed =    0.0064011973; grid_points[1968].w_total = 0.0063905489; grid_points[1968].i_atom = 1;
    grid_points[1969].x =  -0.1366077101; grid_points[1969].y =  -0.4098231304; grid_points[1969].z =   4.3317994029; grid_points[1969].w_fixed =    0.0064011973; grid_points[1969].w_total = 0.0064011513; grid_points[1969].i_atom = 1;
    grid_points[1970].x =  -0.1366077101; grid_points[1970].y =  -0.4098231304; grid_points[1970].z =   4.0585839826; grid_points[1970].w_fixed =    0.0064011973; grid_points[1970].w_total = 0.0064006801; grid_points[1970].i_atom = 1;
    grid_points[1971].x =  -0.4098231304; grid_points[1971].y =  -0.1366077101; grid_points[1971].z =   4.3317994029; grid_points[1971].w_fixed =    0.0064011973; grid_points[1971].w_total = 0.0064011513; grid_points[1971].i_atom = 1;
    grid_points[1972].x =  -0.4098231304; grid_points[1972].y =  -0.1366077101; grid_points[1972].z =   4.0585839826; grid_points[1972].w_fixed =    0.0064011973; grid_points[1972].w_total = 0.0064006801; grid_points[1972].i_atom = 1;
    grid_points[1973].x =  -0.5931967983; grid_points[1973].y =   0.0000000000; grid_points[1973].z =   4.1951916927; grid_points[1973].w_fixed =    0.0089219324; grid_points[1973].w_total = 0.0089213452; grid_points[1973].i_atom = 1;
    grid_points[1974].x =   0.0000000000; grid_points[1974].y =  -0.5931967983; grid_points[1974].z =   4.1951916927; grid_points[1974].w_fixed =    0.0089219324; grid_points[1974].w_total = 0.0089213452; grid_points[1974].i_atom = 1;
    grid_points[1975].x =   0.0000000000; grid_points[1975].y =   0.0000000000; grid_points[1975].z =   4.7883884911; grid_points[1975].w_fixed =    0.0089219324; grid_points[1975].w_total = 0.0089053102; grid_points[1975].i_atom = 1;
    grid_points[1976].x =   0.0000000000; grid_points[1976].y =   0.0000000000; grid_points[1976].z =   3.6019948944; grid_points[1976].w_fixed =    0.0089219324; grid_points[1976].w_total = 0.0087961505; grid_points[1976].i_atom = 1;
    grid_points[1977].x =   0.0000000000; grid_points[1977].y =  -0.4194534787; grid_points[1977].z =   4.6146451714; grid_points[1977].w_fixed =    0.0158612132; grid_points[1977].w_total = 0.0158535719; grid_points[1977].i_atom = 1;
    grid_points[1978].x =   0.0000000000; grid_points[1978].y =  -0.4194534787; grid_points[1978].z =   3.7757382141; grid_points[1978].w_fixed =    0.0158612132; grid_points[1978].w_total = 0.0157937824; grid_points[1978].i_atom = 1;
    grid_points[1979].x =  -0.4194534787; grid_points[1979].y =   0.0000000000; grid_points[1979].z =   4.6146451714; grid_points[1979].w_fixed =    0.0158612132; grid_points[1979].w_total = 0.0158535719; grid_points[1979].i_atom = 1;
    grid_points[1980].x =  -0.4194534787; grid_points[1980].y =   0.0000000000; grid_points[1980].z =   3.7757382141; grid_points[1980].w_fixed =    0.0158612132; grid_points[1980].w_total = 0.0157937824; grid_points[1980].i_atom = 1;
    grid_points[1981].x =  -0.4194534787; grid_points[1981].y =  -0.4194534787; grid_points[1981].z =   4.1951916927; grid_points[1981].w_fixed =    0.0158612132; grid_points[1981].w_total = 0.0158601692; grid_points[1981].i_atom = 1;
    grid_points[1982].x =  -0.3424823312; grid_points[1982].y =  -0.3424823312; grid_points[1982].z =   4.5376740239; grid_points[1982].w_fixed =    0.0148205147; grid_points[1982].w_total = 0.0148168993; grid_points[1982].i_atom = 1;
    grid_points[1983].x =  -0.3424823312; grid_points[1983].y =  -0.3424823312; grid_points[1983].z =   3.8527093615; grid_points[1983].w_fixed =    0.0148205147; grid_points[1983].w_total = 0.0147863672; grid_points[1983].i_atom = 1;
    grid_points[1984].x =  -0.1788555643; grid_points[1984].y =  -0.1788555643; grid_points[1984].z =   4.7317583855; grid_points[1984].w_fixed =    0.0141738295; grid_points[1984].w_total = 0.0141564043; grid_points[1984].i_atom = 1;
    grid_points[1985].x =  -0.1788555643; grid_points[1985].y =  -0.1788555643; grid_points[1985].z =   3.6586249999; grid_points[1985].w_fixed =    0.0141738295; grid_points[1985].w_total = 0.0140351874; grid_points[1985].i_atom = 1;
    grid_points[1986].x =  -0.1788555643; grid_points[1986].y =  -0.5365666928; grid_points[1986].z =   4.3740472570; grid_points[1986].w_fixed =    0.0141738295; grid_points[1986].w_total = 0.0141731261; grid_points[1986].i_atom = 1;
    grid_points[1987].x =  -0.1788555643; grid_points[1987].y =  -0.5365666928; grid_points[1987].z =   4.0163361285; grid_points[1987].w_fixed =    0.0141738295; grid_points[1987].w_total = 0.0141666532; grid_points[1987].i_atom = 1;
    grid_points[1988].x =  -0.5365666928; grid_points[1988].y =  -0.1788555643; grid_points[1988].z =   4.3740472570; grid_points[1988].w_fixed =    0.0141738295; grid_points[1988].w_total = 0.0141731261; grid_points[1988].i_atom = 1;
    grid_points[1989].x =  -0.5365666928; grid_points[1989].y =  -0.1788555643; grid_points[1989].z =   4.0163361285; grid_points[1989].w_fixed =    0.0141738295; grid_points[1989].w_total = 0.0141666532; grid_points[1989].i_atom = 1;
    grid_points[1990].x =  -0.7747876550; grid_points[1990].y =   0.0000000000; grid_points[1990].z =   4.1951916927; grid_points[1990].w_fixed =    0.0059666325; grid_points[1990].w_total = 0.0059644033; grid_points[1990].i_atom = 1;
    grid_points[1991].x =   0.0000000000; grid_points[1991].y =  -0.7747876550; grid_points[1991].z =   4.1951916927; grid_points[1991].w_fixed =    0.0059666325; grid_points[1991].w_total = 0.0059644033; grid_points[1991].i_atom = 1;
    grid_points[1992].x =   0.0000000000; grid_points[1992].y =   0.0000000000; grid_points[1992].z =   4.9699793477; grid_points[1992].w_fixed =    0.0059666325; grid_points[1992].w_total = 0.0058978180; grid_points[1992].i_atom = 1;
    grid_points[1993].x =   0.0000000000; grid_points[1993].y =   0.0000000000; grid_points[1993].z =   3.4204040378; grid_points[1993].w_fixed =    0.0059666325; grid_points[1993].w_total = 0.0055674494; grid_points[1993].i_atom = 1;
    grid_points[1994].x =  -0.4473238612; grid_points[1994].y =  -0.4473238612; grid_points[1994].z =   4.6425155539; grid_points[1994].w_fixed =    0.0152642382; grid_points[1994].w_total = 0.0152405851; grid_points[1994].i_atom = 1;
    grid_points[1995].x =  -0.4473238612; grid_points[1995].y =  -0.4473238612; grid_points[1995].z =   3.7478678316; grid_points[1995].w_fixed =    0.0152642382; grid_points[1995].w_total = 0.0150785295; grid_points[1995].i_atom = 1;
    grid_points[1996].x =  -0.1434253090; grid_points[1996].y =  -0.1434253090; grid_points[1996].z =   4.9429578806; grid_points[1996].w_fixed =    0.0127985780; grid_points[1996].w_total = 0.0126709482; grid_points[1996].i_atom = 1;
    grid_points[1997].x =  -0.1434253090; grid_points[1997].y =  -0.1434253090; grid_points[1997].z =   3.4474255048; grid_points[1997].w_fixed =    0.0127985780; grid_points[1997].w_total = 0.0120391551; grid_points[1997].i_atom = 1;
    grid_points[1998].x =  -0.1434253090; grid_points[1998].y =  -0.7477661879; grid_points[1998].z =   4.3386170017; grid_points[1998].w_fixed =    0.0127985780; grid_points[1998].w_total = 0.0127960888; grid_points[1998].i_atom = 1;
    grid_points[1999].x =  -0.1434253090; grid_points[1999].y =  -0.7477661879; grid_points[1999].z =   4.0517663837; grid_points[1999].w_fixed =    0.0127985780; grid_points[1999].w_total = 0.0127814266; grid_points[1999].i_atom = 1;
    grid_points[2000].x =  -0.7477661879; grid_points[2000].y =  -0.1434253090; grid_points[2000].z =   4.3386170017; grid_points[2000].w_fixed =    0.0127985780; grid_points[2000].w_total = 0.0127960888; grid_points[2000].i_atom = 1;
    grid_points[2001].x =  -0.7477661879; grid_points[2001].y =  -0.1434253090; grid_points[2001].z =   4.0517663837; grid_points[2001].w_fixed =    0.0127985780; grid_points[2001].w_total = 0.0127814266; grid_points[2001].i_atom = 1;
    grid_points[2002].x =  -0.5349297050; grid_points[2002].y =  -0.5349297050; grid_points[2002].z =   4.3625127365; grid_points[2002].w_fixed =    0.0154965859; grid_points[2002].w_total = 0.0154933439; grid_points[2002].i_atom = 1;
    grid_points[2003].x =  -0.5349297050; grid_points[2003].y =  -0.5349297050; grid_points[2003].z =   4.0278706490; grid_points[2003].w_fixed =    0.0154965859; grid_points[2003].w_total = 0.0154711868; grid_points[2003].i_atom = 1;
    grid_points[2004].x =  -0.5349297050; grid_points[2004].y =  -0.1673210437; grid_points[2004].z =   4.7301213978; grid_points[2004].w_fixed =    0.0154965859; grid_points[2004].w_total = 0.0154531108; grid_points[2004].i_atom = 1;
    grid_points[2005].x =  -0.5349297050; grid_points[2005].y =  -0.1673210437; grid_points[2005].z =   3.6602619877; grid_points[2005].w_fixed =    0.0154965859; grid_points[2005].w_total = 0.0151815178; grid_points[2005].i_atom = 1;
    grid_points[2006].x =  -0.1673210437; grid_points[2006].y =  -0.5349297050; grid_points[2006].z =   4.7301213978; grid_points[2006].w_fixed =    0.0154965859; grid_points[2006].w_total = 0.0154531108; grid_points[2006].i_atom = 1;
    grid_points[2007].x =  -0.1673210437; grid_points[2007].y =  -0.5349297050; grid_points[2007].z =   3.6602619877; grid_points[2007].w_fixed =    0.0154965859; grid_points[2007].w_total = 0.0151815178; grid_points[2007].i_atom = 1;
    grid_points[2008].x =  -0.3065753189; grid_points[2008].y =  -0.3065753189; grid_points[2008].z =   4.8373124430; grid_points[2008].w_fixed =    0.0149552262; grid_points[2008].w_total = 0.0148734574; grid_points[2008].i_atom = 1;
    grid_points[2009].x =  -0.3065753189; grid_points[2009].y =  -0.3065753189; grid_points[2009].z =   3.5530709425; grid_points[2009].w_fixed =    0.0149552262; grid_points[2009].w_total = 0.0144183652; grid_points[2009].i_atom = 1;
    grid_points[2010].x =  -0.3065753189; grid_points[2010].y =  -0.6421207503; grid_points[2010].z =   4.5017670117; grid_points[2010].w_fixed =    0.0149552262; grid_points[2010].w_total = 0.0149471131; grid_points[2010].i_atom = 1;
    grid_points[2011].x =  -0.3065753189; grid_points[2011].y =  -0.6421207503; grid_points[2011].z =   3.8886163738; grid_points[2011].w_fixed =    0.0149552262; grid_points[2011].w_total = 0.0148835083; grid_points[2011].i_atom = 1;
    grid_points[2012].x =  -0.6421207503; grid_points[2012].y =  -0.3065753189; grid_points[2012].z =   4.5017670117; grid_points[2012].w_fixed =    0.0149552262; grid_points[2012].w_total = 0.0149471131; grid_points[2012].i_atom = 1;
    grid_points[2013].x =  -0.6421207503; grid_points[2013].y =  -0.3065753189; grid_points[2013].z =   3.8886163738; grid_points[2013].w_fixed =    0.0149552262; grid_points[2013].w_total = 0.0148835083; grid_points[2013].i_atom = 1;
    grid_points[2014].x =  -0.3706344180; grid_points[2014].y =  -0.6803866830; grid_points[2014].z =   4.1951916927; grid_points[2014].w_fixed =    0.0151103430; grid_points[2014].w_total = 0.0151046974; grid_points[2014].i_atom = 1;
    grid_points[2015].x =  -0.6803866830; grid_points[2015].y =  -0.3706344180; grid_points[2015].z =   4.1951916927; grid_points[2015].w_fixed =    0.0151103430; grid_points[2015].w_total = 0.0151046974; grid_points[2015].i_atom = 1;
    grid_points[2016].x =  -0.3706344180; grid_points[2016].y =   0.0000000000; grid_points[2016].z =   4.8755783758; grid_points[2016].w_fixed =    0.0151103430; grid_points[2016].w_total = 0.0150070014; grid_points[2016].i_atom = 1;
    grid_points[2017].x =  -0.3706344180; grid_points[2017].y =   0.0000000000; grid_points[2017].z =   3.5148050097; grid_points[2017].w_fixed =    0.0151103430; grid_points[2017].w_total = 0.0144554995; grid_points[2017].i_atom = 1;
    grid_points[2018].x =  -0.6803866830; grid_points[2018].y =   0.0000000000; grid_points[2018].z =   4.5658261108; grid_points[2018].w_fixed =    0.0151103430; grid_points[2018].w_total = 0.0150969464; grid_points[2018].i_atom = 1;
    grid_points[2019].x =  -0.6803866830; grid_points[2019].y =   0.0000000000; grid_points[2019].z =   3.8245572747; grid_points[2019].w_fixed =    0.0151103430; grid_points[2019].w_total = 0.0149977013; grid_points[2019].i_atom = 1;
    grid_points[2020].x =   0.0000000000; grid_points[2020].y =  -0.3706344180; grid_points[2020].z =   4.8755783758; grid_points[2020].w_fixed =    0.0151103430; grid_points[2020].w_total = 0.0150070014; grid_points[2020].i_atom = 1;
    grid_points[2021].x =   0.0000000000; grid_points[2021].y =  -0.3706344180; grid_points[2021].z =   3.5148050097; grid_points[2021].w_fixed =    0.0151103430; grid_points[2021].w_total = 0.0144554995; grid_points[2021].i_atom = 1;
    grid_points[2022].x =   0.0000000000; grid_points[2022].y =  -0.6803866830; grid_points[2022].z =   4.5658261108; grid_points[2022].w_fixed =    0.0151103430; grid_points[2022].w_total = 0.0150969464; grid_points[2022].i_atom = 1;
    grid_points[2023].x =   0.0000000000; grid_points[2023].y =  -0.6803866830; grid_points[2023].z =   3.8245572747; grid_points[2023].w_fixed =    0.0151103430; grid_points[2023].w_total = 0.0149977013; grid_points[2023].i_atom = 1;
    grid_points[2024].x =  -1.0119675493; grid_points[2024].y =   0.0000000000; grid_points[2024].z =   4.1951916927; grid_points[2024].w_fixed =    0.0133541260; grid_points[2024].w_total = 0.0133291503; grid_points[2024].i_atom = 1;
    grid_points[2025].x =   0.0000000000; grid_points[2025].y =  -1.0119675493; grid_points[2025].z =   4.1951916927; grid_points[2025].w_fixed =    0.0133541260; grid_points[2025].w_total = 0.0133291503; grid_points[2025].i_atom = 1;
    grid_points[2026].x =   0.0000000000; grid_points[2026].y =   0.0000000000; grid_points[2026].z =   5.2071592421; grid_points[2026].w_fixed =    0.0133541260; grid_points[2026].w_total = 0.0125498586; grid_points[2026].i_atom = 1;
    grid_points[2027].x =   0.0000000000; grid_points[2027].y =   0.0000000000; grid_points[2027].z =   3.1832241434; grid_points[2027].w_fixed =    0.0133541260; grid_points[2027].w_total = 0.0101223600; grid_points[2027].i_atom = 1;
    grid_points[2028].x =  -0.5842597370; grid_points[2028].y =  -0.5842597370; grid_points[2028].z =   4.7794514298; grid_points[2028].w_fixed =    0.0341634180; grid_points[2028].w_total = 0.0338660424; grid_points[2028].i_atom = 1;
    grid_points[2029].x =  -0.5842597370; grid_points[2029].y =  -0.5842597370; grid_points[2029].z =   3.6109319557; grid_points[2029].w_fixed =    0.0341634180; grid_points[2029].w_total = 0.0323366953; grid_points[2029].i_atom = 1;
    grid_points[2030].x =  -0.1873310158; grid_points[2030].y =  -0.1873310158; grid_points[2030].z =   5.1718658973; grid_points[2030].w_fixed =    0.0286449390; grid_points[2030].w_total = 0.0271394907; grid_points[2030].i_atom = 1;
    grid_points[2031].x =  -0.1873310158; grid_points[2031].y =  -0.1873310158; grid_points[2031].z =   3.2185174881; grid_points[2031].w_fixed =    0.0286449390; grid_points[2031].w_total = 0.0223716243; grid_points[2031].i_atom = 1;
    grid_points[2032].x =  -0.1873310158; grid_points[2032].y =  -0.9766742046; grid_points[2032].z =   4.3825227086; grid_points[2032].w_fixed =    0.0286449390; grid_points[2032].w_total = 0.0286147757; grid_points[2032].i_atom = 1;
    grid_points[2033].x =  -0.1873310158; grid_points[2033].y =  -0.9766742046; grid_points[2033].z =   4.0078606769; grid_points[2033].w_fixed =    0.0286449390; grid_points[2033].w_total = 0.0284587127; grid_points[2033].i_atom = 1;
    grid_points[2034].x =  -0.9766742046; grid_points[2034].y =  -0.1873310158; grid_points[2034].z =   4.3825227086; grid_points[2034].w_fixed =    0.0286449390; grid_points[2034].w_total = 0.0286147757; grid_points[2034].i_atom = 1;
    grid_points[2035].x =  -0.9766742046; grid_points[2035].y =  -0.1873310158; grid_points[2035].z =   4.0078606769; grid_points[2035].w_fixed =    0.0286449390; grid_points[2035].w_total = 0.0284587127; grid_points[2035].i_atom = 1;
    grid_points[2036].x =  -0.6986836963; grid_points[2036].y =  -0.6986836963; grid_points[2036].z =   4.4137334641; grid_points[2036].w_fixed =    0.0346834435; grid_points[2036].w_total = 0.0346436110; grid_points[2036].i_atom = 1;
    grid_points[2037].x =  -0.6986836963; grid_points[2037].y =  -0.6986836963; grid_points[2037].z =   3.9766499214; grid_points[2037].w_fixed =    0.0346834435; grid_points[2037].w_total = 0.0344092150; grid_points[2037].i_atom = 1;
    grid_points[2038].x =  -0.6986836963; grid_points[2038].y =  -0.2185417714; grid_points[2038].z =   4.8938753891; grid_points[2038].w_fixed =    0.0346834435; grid_points[2038].w_total = 0.0341432423; grid_points[2038].i_atom = 1;
    grid_points[2039].x =  -0.6986836963; grid_points[2039].y =  -0.2185417714; grid_points[2039].z =   3.4965079964; grid_points[2039].w_fixed =    0.0346834435; grid_points[2039].w_total = 0.0317117114; grid_points[2039].i_atom = 1;
    grid_points[2040].x =  -0.2185417714; grid_points[2040].y =  -0.6986836963; grid_points[2040].z =   4.8938753891; grid_points[2040].w_fixed =    0.0346834435; grid_points[2040].w_total = 0.0341432423; grid_points[2040].i_atom = 1;
    grid_points[2041].x =  -0.2185417714; grid_points[2041].y =  -0.6986836963; grid_points[2041].z =   3.4965079964; grid_points[2041].w_fixed =    0.0346834435; grid_points[2041].w_total = 0.0317117114; grid_points[2041].i_atom = 1;
    grid_points[2042].x =  -0.4004249063; grid_points[2042].y =  -0.4004249063; grid_points[2042].z =   5.0338800196; grid_points[2042].w_fixed =    0.0334718077; grid_points[2042].w_total = 0.0324778009; grid_points[2042].i_atom = 1;
    grid_points[2043].x =  -0.4004249063; grid_points[2043].y =  -0.4004249063; grid_points[2043].z =   3.3565033658; grid_points[2043].w_fixed =    0.0334718077; grid_points[2043].w_total = 0.0287086336; grid_points[2043].i_atom = 1;
    grid_points[2044].x =  -0.4004249063; grid_points[2044].y =  -0.8386883269; grid_points[2044].z =   4.5956165991; grid_points[2044].w_fixed =    0.0334718077; grid_points[2044].w_total = 0.0333691081; grid_points[2044].i_atom = 1;
    grid_points[2045].x =  -0.4004249063; grid_points[2045].y =  -0.8386883269; grid_points[2045].z =   3.7947667864; grid_points[2045].w_fixed =    0.0334718077; grid_points[2045].w_total = 0.0327270168; grid_points[2045].i_atom = 1;
    grid_points[2046].x =  -0.8386883269; grid_points[2046].y =  -0.4004249063; grid_points[2046].z =   4.5956165991; grid_points[2046].w_fixed =    0.0334718077; grid_points[2046].w_total = 0.0333691081; grid_points[2046].i_atom = 1;
    grid_points[2047].x =  -0.8386883269; grid_points[2047].y =  -0.4004249063; grid_points[2047].z =   3.7947667864; grid_points[2047].w_fixed =    0.0334718077; grid_points[2047].w_total = 0.0327270168; grid_points[2047].i_atom = 1;
    grid_points[2048].x =  -0.4840939338; grid_points[2048].y =  -0.8886683207; grid_points[2048].z =   4.1951916927; grid_points[2048].w_fixed =    0.0338189800; grid_points[2048].w_total = 0.0337557378; grid_points[2048].i_atom = 1;
    grid_points[2049].x =  -0.8886683207; grid_points[2049].y =  -0.4840939338; grid_points[2049].z =   4.1951916927; grid_points[2049].w_fixed =    0.0338189800; grid_points[2049].w_total = 0.0337557378; grid_points[2049].i_atom = 1;
    grid_points[2050].x =  -0.4840939338; grid_points[2050].y =   0.0000000000; grid_points[2050].z =   5.0838600134; grid_points[2050].w_fixed =    0.0338189800; grid_points[2050].w_total = 0.0325750639; grid_points[2050].i_atom = 1;
    grid_points[2051].x =  -0.4840939338; grid_points[2051].y =   0.0000000000; grid_points[2051].z =   3.3065233720; grid_points[2051].w_fixed =    0.0338189800; grid_points[2051].w_total = 0.0281497735; grid_points[2051].i_atom = 1;
    grid_points[2052].x =  -0.8886683207; grid_points[2052].y =   0.0000000000; grid_points[2052].z =   4.6792856265; grid_points[2052].w_fixed =    0.0338189800; grid_points[2052].w_total = 0.0336495489; grid_points[2052].i_atom = 1;
    grid_points[2053].x =  -0.8886683207; grid_points[2053].y =   0.0000000000; grid_points[2053].z =   3.7110977590; grid_points[2053].w_fixed =    0.0338189800; grid_points[2053].w_total = 0.0326753125; grid_points[2053].i_atom = 1;
    grid_points[2054].x =   0.0000000000; grid_points[2054].y =  -0.4840939338; grid_points[2054].z =   5.0838600134; grid_points[2054].w_fixed =    0.0338189800; grid_points[2054].w_total = 0.0325750639; grid_points[2054].i_atom = 1;
    grid_points[2055].x =   0.0000000000; grid_points[2055].y =  -0.4840939338; grid_points[2055].z =   3.3065233720; grid_points[2055].w_fixed =    0.0338189800; grid_points[2055].w_total = 0.0281497735; grid_points[2055].i_atom = 1;
    grid_points[2056].x =   0.0000000000; grid_points[2056].y =  -0.8886683207; grid_points[2056].z =   4.6792856265; grid_points[2056].w_fixed =    0.0338189800; grid_points[2056].w_total = 0.0336495489; grid_points[2056].i_atom = 1;
    grid_points[2057].x =   0.0000000000; grid_points[2057].y =  -0.8886683207; grid_points[2057].z =   3.7110977590; grid_points[2057].w_fixed =    0.0338189800; grid_points[2057].w_total = 0.0326753125; grid_points[2057].i_atom = 1;
    grid_points[2058].x =   0.0000000000; grid_points[2058].y =   0.0000000000; grid_points[2058].z =   5.5201244282; grid_points[2058].w_fixed =    0.0303775643; grid_points[2058].w_total = 0.0229617486; grid_points[2058].i_atom = 1;
    grid_points[2059].x =  -0.7649502715; grid_points[2059].y =  -0.7649502715; grid_points[2059].z =   4.9601419642; grid_points[2059].w_fixed =    0.0777139157; grid_points[2059].w_total = 0.0744471256; grid_points[2059].i_atom = 1;
    grid_points[2060].x =  -0.7649502715; grid_points[2060].y =  -0.7649502715; grid_points[2060].z =   3.4302414213; grid_points[2060].w_fixed =    0.0777139157; grid_points[2060].w_total = 0.0632755179; grid_points[2060].i_atom = 1;
    grid_points[2061].x =  -0.2452657651; grid_points[2061].y =  -0.2452657651; grid_points[2061].z =   5.4739161206; grid_points[2061].w_fixed =    0.0651606456; grid_points[2061].w_total = 0.0509705226; grid_points[2061].i_atom = 1;
    grid_points[2062].x =  -0.9147614482; grid_points[2062].y =  -0.9147614482; grid_points[2062].z =   4.4813205782; grid_points[2062].w_fixed =    0.0788968540; grid_points[2062].w_total = 0.0784571533; grid_points[2062].i_atom = 1;
    grid_points[2063].x =  -0.9147614482; grid_points[2063].y =  -0.9147614482; grid_points[2063].z =   3.9090628073; grid_points[2063].w_fixed =    0.0788968540; grid_points[2063].w_total = 0.0763502659; grid_points[2063].i_atom = 1;
    grid_points[2064].x =  -0.9147614482; grid_points[2064].y =  -0.2861288854; grid_points[2064].z =   5.1099531410; grid_points[2064].w_fixed =    0.0788968540; grid_points[2064].w_total = 0.0731312984; grid_points[2064].i_atom = 1;
    grid_points[2065].x =  -0.9147614482; grid_points[2065].y =  -0.2861288854; grid_points[2065].z =   3.2804302445; grid_points[2065].w_fixed =    0.0788968540; grid_points[2065].w_total = 0.0571063680; grid_points[2065].i_atom = 1;
    grid_points[2066].x =  -0.2861288854; grid_points[2066].y =  -0.9147614482; grid_points[2066].z =   5.1099531410; grid_points[2066].w_fixed =    0.0788968540; grid_points[2066].w_total = 0.0731312984; grid_points[2066].i_atom = 1;
    grid_points[2067].x =  -0.2861288854; grid_points[2067].y =  -0.9147614482; grid_points[2067].z =   3.2804302445; grid_points[2067].w_fixed =    0.0788968540; grid_points[2067].w_total = 0.0571063680; grid_points[2067].i_atom = 1;
    grid_points[2068].x =  -0.5242619359; grid_points[2068].y =  -0.5242619359; grid_points[2068].z =   5.2932561710; grid_points[2068].w_fixed =    0.0761406613; grid_points[2068].w_total = 0.0660727372; grid_points[2068].i_atom = 1;
    grid_points[2069].x =  -0.5242619359; grid_points[2069].y =  -1.0980644783; grid_points[2069].z =   4.7194536287; grid_points[2069].w_fixed =    0.0761406613; grid_points[2069].w_total = 0.0749832017; grid_points[2069].i_atom = 1;
    grid_points[2070].x =  -0.5242619359; grid_points[2070].y =  -1.0980644783; grid_points[2070].z =   3.6709297568; grid_points[2070].w_fixed =    0.0761406613; grid_points[2070].w_total = 0.0696636336; grid_points[2070].i_atom = 1;
    grid_points[2071].x =  -1.0980644783; grid_points[2071].y =  -0.5242619359; grid_points[2071].z =   4.7194536287; grid_points[2071].w_fixed =    0.0761406613; grid_points[2071].w_total = 0.0749832017; grid_points[2071].i_atom = 1;
    grid_points[2072].x =  -1.0980644783; grid_points[2072].y =  -0.5242619359; grid_points[2072].z =   3.6709297568; grid_points[2072].w_fixed =    0.0761406613; grid_points[2072].w_total = 0.0696636336; grid_points[2072].i_atom = 1;
    grid_points[2073].x =  -0.6338067859; grid_points[2073].y =  -1.1635014876; grid_points[2073].z =   4.1951916927; grid_points[2073].w_fixed =    0.0769303986; grid_points[2073].w_total = 0.0763061913; grid_points[2073].i_atom = 1;
    grid_points[2074].x =  -1.1635014876; grid_points[2074].y =  -0.6338067859; grid_points[2074].z =   4.1951916927; grid_points[2074].w_fixed =    0.0769303986; grid_points[2074].w_total = 0.0763061913; grid_points[2074].i_atom = 1;
    grid_points[2075].x =  -0.6338067859; grid_points[2075].y =   0.0000000000; grid_points[2075].z =   5.3586931803; grid_points[2075].w_fixed =    0.0769303986; grid_points[2075].w_total = 0.0646229177; grid_points[2075].i_atom = 1;
    grid_points[2076].x =  -1.1635014876; grid_points[2076].y =   0.0000000000; grid_points[2076].z =   4.8289984786; grid_points[2076].w_fixed =    0.0769303986; grid_points[2076].w_total = 0.0750352153; grid_points[2076].i_atom = 1;
    grid_points[2077].x =  -1.1635014876; grid_points[2077].y =   0.0000000000; grid_points[2077].z =   3.5613849069; grid_points[2077].w_fixed =    0.0769303986; grid_points[2077].w_total = 0.0673622975; grid_points[2077].i_atom = 1;
    grid_points[2078].x =   0.0000000000; grid_points[2078].y =  -0.6338067859; grid_points[2078].z =   5.3586931803; grid_points[2078].w_fixed =    0.0769303986; grid_points[2078].w_total = 0.0646229177; grid_points[2078].i_atom = 1;
    grid_points[2079].x =   0.0000000000; grid_points[2079].y =  -1.1635014876; grid_points[2079].z =   4.8289984786; grid_points[2079].w_fixed =    0.0769303986; grid_points[2079].w_total = 0.0750352153; grid_points[2079].i_atom = 1;
    grid_points[2080].x =   0.0000000000; grid_points[2080].y =  -1.1635014876; grid_points[2080].z =   3.5613849069; grid_points[2080].w_fixed =    0.0769303986; grid_points[2080].w_total = 0.0673622975; grid_points[2080].i_atom = 1;
    grid_points[2081].x =   0.0000000000; grid_points[2081].y =   0.0000000000; grid_points[2081].z =   5.9384639164; grid_points[2081].w_fixed =    0.0707954930; grid_points[2081].w_total = 0.0251826777; grid_points[2081].i_atom = 1;
    grid_points[2082].x =  -1.0064786876; grid_points[2082].y =  -1.0064786876; grid_points[2082].z =   5.2016703804; grid_points[2082].w_fixed =    0.1811137632; grid_points[2082].w_total = 0.1504983152; grid_points[2082].i_atom = 1;
    grid_points[2083].x =  -1.0064786876; grid_points[2083].y =  -1.0064786876; grid_points[2083].z =   3.1887130051; grid_points[2083].w_fixed =    0.1811137632; grid_points[2083].w_total = 0.0948100401; grid_points[2083].i_atom = 1;
    grid_points[2084].x =  -0.3227069453; grid_points[2084].y =  -0.3227069453; grid_points[2084].z =   5.8776656155; grid_points[2084].w_fixed =    0.1518581277; grid_points[2084].w_total = 0.0603864433; grid_points[2084].i_atom = 1;
    grid_points[2085].x =  -1.2035918363; grid_points[2085].y =  -1.2035918363; grid_points[2085].z =   4.5716640411; grid_points[2085].w_fixed =    0.1838706234; grid_points[2085].w_total = 0.1794730010; grid_points[2085].i_atom = 1;
    grid_points[2086].x =  -1.2035918363; grid_points[2086].y =  -1.2035918363; grid_points[2086].z =   3.8187193444; grid_points[2086].w_fixed =    0.1838706234; grid_points[2086].w_total = 0.1635857757; grid_points[2086].i_atom = 1;
    grid_points[2087].x =  -1.2035918363; grid_points[2087].y =  -0.3764723484; grid_points[2087].z =   5.3987835290; grid_points[2087].w_fixed =    0.1838706234; grid_points[2087].w_total = 0.1337193515; grid_points[2087].i_atom = 1;
    grid_points[2088].x =  -0.3764723484; grid_points[2088].y =  -1.2035918363; grid_points[2088].z =   5.3987835290; grid_points[2088].w_fixed =    0.1838706234; grid_points[2088].w_total = 0.1337193515; grid_points[2088].i_atom = 1;
    grid_points[2089].x =  -0.6897944676; grid_points[2089].y =  -0.6897944676; grid_points[2089].z =   5.6399633809; grid_points[2089].w_fixed =    0.1774472638; grid_points[2089].w_total = 0.1007712847; grid_points[2089].i_atom = 1;
    grid_points[2090].x =  -0.8339274406; grid_points[2090].y =   0.0000000000; grid_points[2090].z =   5.7260617296; grid_points[2090].w_fixed =    0.1792877615; grid_points[2090].w_total = 0.0907093839; grid_points[2090].i_atom = 1;
    grid_points[2091].x =   0.0000000000; grid_points[2091].y =  -0.8339274406; grid_points[2091].z =   5.7260617296; grid_points[2091].w_fixed =    0.1792877615; grid_points[2091].w_total = 0.0907093839; grid_points[2091].i_atom = 1;
    grid_points[2092].x =  -0.9146585961; grid_points[2092].y =  -0.9146585961; grid_points[2092].z =   6.1109403775; grid_points[2092].w_fixed =    0.4275576805; grid_points[2092].w_total = 0.0649370027; grid_points[2092].i_atom = 1;
    grid_points[2093].x =  -1.1057770654; grid_points[2093].y =   0.0000000000; grid_points[2093].z =   6.2251056810; grid_points[2093].w_fixed =    0.4319923443; grid_points[2093].w_total = 0.0483760601; grid_points[2093].i_atom = 1;
    grid_points[2094].x =   0.0000000000; grid_points[2094].y =  -1.1057770654; grid_points[2094].z =   6.2251056810; grid_points[2094].w_fixed =    0.4319923443; grid_points[2094].w_total = 0.0483760601; grid_points[2094].i_atom = 1;
    grid_points[2095].x =   0.0000000000; grid_points[2095].y =   8.3645442750; grid_points[2095].z =   4.1951916927; grid_points[2095].w_fixed =   34.8025870819; grid_points[2095].w_total = 0.0000000000; grid_points[2095].i_atom = 1;
    grid_points[2096].x =   8.3645442750; grid_points[2096].y =   0.0000000000; grid_points[2096].z =   4.1951916927; grid_points[2096].w_fixed =   34.8025870819; grid_points[2096].w_total = 0.0000000000; grid_points[2096].i_atom = 1;
    grid_points[2097].x =  -1.7666549613; grid_points[2097].y =  -1.7666549613; grid_points[2097].z =  -1.1047731912; grid_points[2097].w_fixed =   17.3848280115; grid_points[2097].w_total = 0.0000000000; grid_points[2097].i_atom = 1;
    grid_points[2098].x =  -2.4354299108; grid_points[2098].y =  -2.4354299108; grid_points[2098].z =   6.6306216035; grid_points[2098].w_fixed =    6.3162875542; grid_points[2098].w_total = 0.0011258350; grid_points[2098].i_atom = 1;
    grid_points[2099].x =  -1.2718617903; grid_points[2099].y =  -1.2718617903; grid_points[2099].z =   8.0107770637; grid_points[2099].w_fixed =    6.0406797362; grid_points[2099].w_total = 0.0000010065; grid_points[2099].i_atom = 1;
    grid_points[2100].x =  -2.1397188201; grid_points[2100].y =  -0.6692841749; grid_points[2100].z =   6.3349105128; grid_points[2100].w_fixed =    1.1157541830; grid_points[2100].w_total = 0.0256935379; grid_points[2100].i_atom = 1;
    grid_points[2101].x =  -1.4825376722; grid_points[2101].y =   0.0000000000; grid_points[2101].z =   6.9167384249; grid_points[2101].w_fixed =    1.0879446980; grid_points[2101].w_total = 0.0030241246; grid_points[2101].i_atom = 1;
    grid_points[2102].x =  -2.9827802928; grid_points[2102].y =   0.0000000000; grid_points[2102].z =   7.1779719856; grid_points[2102].w_fixed =    6.7598180997; grid_points[2102].w_total = 0.0001667045; grid_points[2102].i_atom = 1;
    grid_points[2103].x =  -2.5220049911; grid_points[2103].y =   7.5660149734; grid_points[2103].z =   6.7171966839; grid_points[2103].w_fixed =   55.2891360098; grid_points[2103].w_total = 0.0000000000; grid_points[2103].i_atom = 1;
    grid_points[2104].x =  -0.6692841749; grid_points[2104].y =  -2.1397188201; grid_points[2104].z =   6.3349105128; grid_points[2104].w_fixed =    1.1157541830; grid_points[2104].w_total = 0.0256935379; grid_points[2104].i_atom = 1;
    grid_points[2105].x =   0.0000000000; grid_points[2105].y =  -1.4825376722; grid_points[2105].z =   6.9167384249; grid_points[2105].w_fixed =    1.0879446980; grid_points[2105].w_total = 0.0030241246; grid_points[2105].i_atom = 1;
    grid_points[2106].x =   0.0000000000; grid_points[2106].y =  -2.9827802928; grid_points[2106].z =   7.1779719856; grid_points[2106].w_fixed =    6.7598180997; grid_points[2106].w_total = 0.0001667045; grid_points[2106].i_atom = 1;
    grid_points[2107].x =   0.0000000000; grid_points[2107].y =   0.0000000000; grid_points[2107].z =   6.5067482501; grid_points[2107].w_fixed =    0.1705811414; grid_points[2107].w_total = 0.0077825024; grid_points[2107].i_atom = 1;
    grid_points[2108].x =  -0.4279052608; grid_points[2108].y =  -0.4279052608; grid_points[2108].z =   6.4261304847; grid_points[2108].w_fixed =    0.3659008737; grid_points[2108].w_total = 0.0221696713; grid_points[2108].i_atom = 1;
    grid_points[2109].x =   0.0000000000; grid_points[2109].y =   0.0000000000; grid_points[2109].z =   7.2943423126; grid_points[2109].w_fixed =    0.4295975400; grid_points[2109].w_total = 0.0001247332; grid_points[2109].i_atom = 1;
    grid_points[2110].x =  -0.5737012360; grid_points[2110].y =  -0.5737012360; grid_points[2110].z =   7.1862564443; grid_points[2110].w_fixed =    0.9214976165; grid_points[2110].w_total = 0.0005834114; grid_points[2110].i_atom = 1;
    grid_points[2111].x =  -1.2263012757; grid_points[2111].y =  -1.2263012757; grid_points[2111].z =   6.7636746939; grid_points[2111].w_fixed =    1.0767762850; grid_points[2111].w_total = 0.0058282274; grid_points[2111].i_atom = 1;
    grid_points[2112].x =   0.0000000000; grid_points[2112].y =   0.0000000000; grid_points[2112].z =   8.4134800364; grid_points[2112].w_fixed =    3.8023976811; grid_points[2112].w_total = 0.0000000004; grid_points[2112].i_atom = 1;
    grid_points[2113].x =   7.5660149734; grid_points[2113].y =  -2.5220049911; grid_points[2113].z =   6.7171966839; grid_points[2113].w_fixed =   55.2891360098; grid_points[2113].w_total = 0.0000000000; grid_points[2113].i_atom = 1;

    std::stable_sort(grid_points, grid_points + n_grid_point, [](const PeriodicBox::GridPoint& a, const PeriodicBox::GridPoint& b) { return a.i_atom < b.i_atom; });

    const double epsilon_xc[n_grid_point] {
        -0.0000269544,
        -0.0004867574,
        -0.0000674222,
        -0.0000675967,
        -0.0000133146,
        -0.0077243539,
        -0.0000501945,
        -0.0000346158,
        -0.0000479764,
        -0.0000289821,
        -0.0000500261,
        -3.6003766093,
        -0.0005650258,
        -0.0000773039,
        -0.0004867574,
        -0.0000674222,
        -0.0000346158,
        -0.0000479764,
        -0.0000289821,
        -0.0000500261,
        -3.6003766093,
        -0.0005650258,
        -0.0078611060,
        -0.0236439453,
        -0.0034082519,
        -3.5973905092,
        -0.0000606923,
       -22.9853780414,
       -22.1567476341,
       -19.1065918217,
       -14.4474167447,
        -9.8824303418,
        -5.8462697009,
        -2.9658756029,
        -1.2404101115,
        -0.4202151511,
        -0.4202156232,
        -0.4197645374,
        -0.1089010630,
        -0.1089021272,
        -0.1090208263,
        -0.0208418730,
        -0.0208441904,
        -0.0212511776,
        -0.0025689291,
        -0.0025735857,
        -0.0030086230,
        -0.0002421192,
        -0.0002490056,
        -0.0005468296,
        -0.0021571094,
        -0.0003174064,
        -0.0003174064,
        -0.0000378659,
        -0.0000447885,
        -0.0002749821,
        -0.0041488022,
        -0.0000791729,
        -0.0000791729,
        -0.0001833732,
        -0.0000401152,
        -0.0005028912,
        -0.0005028912,
        -0.0029318678,
        -0.0000623181,
        -0.0000623181,
        -0.0000260829,
        -0.0000260829,
        -0.0000694916,
        -0.0494977855,
        -0.0000834514,
        -0.0000269193,
        -0.0002043327,
        -0.0000610911,
        -0.0000165798,
        -0.0000269193,
        -0.0002043327,
        -0.0000222680,
        -0.0000120055,
        -1.2613900419,
        -0.0023652036,
        -0.0002096102,
        -0.0000199676,
        -0.0000610911,
        -0.0000272110,
        -0.0000773039,
        -0.5600697020,
        -0.0002096102,
        -0.0273042020,
        -0.0273042020,
        -0.0658109204,
        -3.0893155072,
        -0.4726866296,
     -1894.0692407526,
        -0.0000238220,
        -0.1183084176,
        -0.0188594708,
        -0.0004708701,
        -0.0013931545,
        -0.0013931545,
        -0.0668996625,
        -0.0002449291,
        -0.0011207016,
        -0.0011207016,
        -0.2025001447,
        -0.0023451236,
        -0.0257083764,
        -0.0003642815,
        -0.0000178538,
        -0.0415289527,
        -0.0000295220,
        -0.0000132105,
        -0.0006725059,
        -0.0000334864,
        -0.0000280605,
        -0.0000936675,
        -0.0000215742,
        -0.0000261126,
        -0.0000354026,
        -0.0000866269,
        -0.0000280605,
        -0.0000936675,
        -0.0000215742,
        -0.0000261126,
        -0.0028868455,
        -0.0000367835,
        -0.0033516448,
        -0.0006503505,
        -0.0003805789,
        -0.0003805789,
        -0.0005311651,
        -0.0008400697,
        -0.0001084618,
        -0.0001084618,
        -0.0000736550,
        -0.0003420026,
        -0.0049976082,
        -0.0000526867,
        -0.0025384627,
        -0.0000845860,
        -0.0000845860,
        -0.0003126054,
        -0.0006043240,
        -0.0000319017,
        -0.0006043240,
        -0.0000339708,
       -22.9855652516,
       -22.1574074287,
       -19.1080453836,
       -14.4500598025,
        -9.8861330887,
        -5.8507070548,
        -2.9704664895,
        -1.2445661768,
        -0.4234855095,
        -0.4222128804,
        -0.4222128804,
        -0.4217327482,
        -0.1111468088,
        -0.1101823670,
        -0.1101823670,
        -0.1098380443,
        -0.0222450065,
        -0.0215729696,
        -0.0215729696,
        -0.0213491600,
        -0.0029222335,
        -0.0029222335,
        -0.0027971328,
        -0.0003237755,
        -0.0002529402,
        -0.0002529402,
        -0.0000414908,
        -0.0000414908,
        -0.0000301290,
        -0.0000312218,
        -0.0000312218,
        -0.0003999500,
        -0.0000366961,
        -0.0000561909,
        -0.0000179837,
        -0.0000872246,
        -0.0000320405,
        -0.0000561909,
        -0.0000179837,
        -0.0000115020,
        -0.0000203492,
        -0.0003955841,
        -1.2646862314,
        -0.0001521951,
        -0.0000872246,
        -0.0000320405,
        -0.0000354026,
        -0.5623624376,
        -0.0001521951,
        -0.0000045867,
        -0.0000045867,
        -0.0000002209,
        -0.0000005760,
        -0.0000017149,
        -0.0000078949,
        -0.0000035843,
        -0.0000030744,
        -0.0000086793,
        -0.0000252407,
        -0.0000150705,
        -0.0000258495,
        -0.0000045423,
        -0.0000058395,
        -0.0000055968,
        -0.0000025932,
        -0.0000035843,
        -0.0000086793,
        -0.0000252407,
        -0.0000150705,
        -0.0000045423,
        -0.0000058395,
        -0.0000055968,
        -0.0000078949,
        -0.0000136044,
        -0.0000136044,
        -0.0000105724,
        -0.0000120418,
        -0.0000105724,
        -0.0000120418,
        -0.0000001362,
        -0.0000000567,
        -0.0000266387,
        -0.0000610837,
        -0.0000728889,
        -0.0000125843,
        -0.0000019007,
        -0.0000171451,
        -0.0000019619,
        -0.0000728889,
        -0.0000125843,
        -0.0000019007,
        -0.0000171451,
        -0.0000019619,
        -0.0000001362,
        -0.0000006147,
        -0.0000007795,
        -0.0000006147,
        -0.0000007795,
        -0.0000124835,
        -0.0000084900,
        -0.0000190384,
        -0.0000154343,
        -0.0000063671,
        -0.0000028552,
        -0.0000067133,
        -0.0000037436,
        -0.0000124835,
        -0.0000190384,
        -0.0000063671,
        -0.0000084900,
        -0.0000067133,
        -0.0000037436,
        -0.0000269544,
        -0.0000773039,
        -0.0000346158,
        -0.0000479764,
        -0.0000500261,
        -0.0005650258,
        -0.0000133146,
        -0.0000675967,
        -0.0077243539,
        -0.0000501945,
        -0.0000199676,
        -0.0000610911,
        -0.0004952328,
        -0.0034082519,
        -0.0000606923,
        -0.0025689291,
        -0.0002421192,
        -0.0008189558,
        -0.0002490056,
        -0.0003174064,
        -0.0000378659,
        -0.0006596125,
        -0.0000447885,
        -0.0002749821,
        -0.0000791729,
        -0.0000165798,
        -0.0001833732,
        -0.0000269193,
        -0.0000401152,
        -0.0005028912,
        -0.0000623181,
        -0.0000260829,
        -0.0000260829,
        -0.0000854280,
        -0.0000694916,
        -0.0002043327,
        -0.0000834514,
       -22.9853780414,
       -22.1567476341,
       -19.1065918217,
       -14.4474167447,
        -9.8824303418,
        -5.8462697009,
        -2.9658756029,
        -1.2404101115,
        -0.4202151511,
        -0.4198278901,
        -0.4198278901,
        -0.4202156232,
        -0.4197645374,
        -0.1089010630,
        -0.1092266384,
        -0.1092266384,
        -0.1089021272,
        -0.1090208263,
        -0.0208418730,
        -0.0215279996,
        -0.0215279996,
        -0.0208441904,
        -0.0212511776,
        -0.0032881155,
        -0.0032881155,
        -0.0025735857,
        -0.0030086230,
        -0.0008189558,
        -0.0005468296,
        -0.0021571094,
        -0.0003174064,
        -0.0006596125,
        -0.0041488022,
        -0.0000791729,
        -0.0005028912,
        -0.0029318678,
        -0.0000623181,
        -0.0063573213,
        -0.0000854280,
        -0.0494977855,
        -0.0000346158,
        -0.0000500261,
        -0.0000479764,
        -0.0005650258,
        -0.0004867574,
        -0.0000222680,
        -0.0000120055,
        -0.0000269193,
        -0.0002043327,
        -0.0002096102,
        -1.2613900419,
        -0.0023652036,
        -0.0000272110,
        -0.0000610911,
        -0.0004952328,
        -0.5600697020,
        -0.0000773039,
        -0.0002096102,
        -0.0658109204,
        -0.0415289527,
        -0.0045662702,
        -0.0000238220,
        -0.0004708701,
        -0.0013931545,
        -0.0002449291,
        -0.0011207016,
        -0.0023451236,
        -0.0003642815,
        -0.0000178538,
        -0.0188594708,
        -0.0668996625,
        -0.2025001447,
        -0.0069502465,
        -0.0257083764,
        -0.0012732225,
        -0.0045662702,
        -0.0000866269,
        -0.0000354026,
        -0.0000280605,
        -0.0000936675,
        -0.0000215742,
        -0.0000261126,
        -0.0000132105,
        -0.0000295220,
        -0.0006725059,
        -0.0000334864,
        -0.0000872246,
        -0.0000320405,
        -0.0000367835,
        -0.0003805789,
        -0.0001084618,
        -0.0000736550,
        -0.0000526867,
        -0.0000845860,
        -0.0006043240,
        -0.0000319017,
        -0.0000339708,
        -0.0000561909,
        -0.0002529402,
        -0.0000414908,
        -0.0000179837,
        -0.0000301290,
        -0.0000312218,
        -0.0003999500,
        -0.0000366961,
        -0.0005311651,
        -0.0003420026,
        -0.0025384627,
        -0.0000845860,
        -0.0003126054,
        -0.4222128804,
        -0.4217327482,
        -0.1101823670,
        -0.1098380443,
        -0.0215729696,
        -0.0213491600,
        -0.0029222335,
        -0.0027971328,
        -0.0003237755,
        -0.0002529402,
        -0.0000414908,
        -0.0000312218,
        -0.0000280605,
        -0.0000215742,
        -0.0000936675,
        -0.0000261126,
        -0.0000115020,
        -0.0000203492,
        -0.0000561909,
        -0.0000179837,
        -0.0001521951,
        -0.0003955841,
        -1.2646862314,
        -0.0000320405,
        -0.5623624376,
        -0.0000354026,
        -0.0001521951,
        -0.0000002209,
        -0.0000005760,
        -0.0000017149,
        -0.0000078949,
        -0.0000136044,
        -0.0000045867,
        -0.0000035843,
        -0.0000030744,
        -0.0000086793,
        -0.0000252407,
        -0.0000150705,
        -0.0000045423,
        -0.0000058395,
        -0.0000055968,
        -0.0000045867,
        -0.0000035843,
        -0.0000030744,
        -0.0000086793,
        -0.0000252407,
        -0.0000150705,
        -0.0000017149,
        -0.0000045423,
        -0.0000058395,
        -0.0000055968,
        -0.0000078949,
        -0.0000136044,
        -0.0000005760,
        -0.0000002209,
        -0.0000001362,
        -0.0000266387,
        -0.0000105724,
        -0.0000120418,
        -0.0000728889,
        -0.0000125843,
        -0.0000019007,
        -0.0000171451,
        -0.0000019619,
        -0.0000006147,
        -0.0000007795,
        -0.0000105724,
        -0.0000120418,
        -0.0000728889,
        -0.0000125843,
        -0.0000019007,
        -0.0000266387,
        -0.0000171451,
        -0.0000019619,
        -0.0000006147,
        -0.0000007795,
        -0.0000001362,
        -0.0000124835,
        -0.0000084900,
        -0.0000190384,
        -0.0000154343,
        -0.0000063671,
        -0.0000028552,
        -0.0000067133,
        -0.0000037436,
        -0.0000124835,
        -0.0000190384,
        -0.0000154343,
        -0.0000063671,
        -0.0000028552,
        -0.0000084900,
        -0.0000067133,
        -0.0000037436,
        -0.0000000739,
        -0.0000000739,
        -0.0000108365,
        -0.0000108365,
        -0.0000108365,
        -0.0000049610,
        -0.0000057629,
        -0.0000000733,
        -0.0000001053,
        -0.0000049610,
        -0.0000057629,
        -0.0000000733,
        -0.0000108365,
        -0.0000001053,
        -0.0000108365,
        -0.0000049610,
        -0.0000057629,
        -0.0000000733,
        -0.0000108365,
        -0.0000001053,
        -0.0000049610,
        -0.0000057629,
        -0.0000000733,
        -0.0000108365,
        -0.0000108365,
        -0.0000001053,
        -0.0000000739,
        -0.0000000739,
        -0.0000012324,
        -0.0000011859,
        -0.0000011111,
        -0.0000011320,
        -0.0000011111,
        -0.0000011320,
        -0.0000012324,
        -0.0000011859,
        -0.0000007484,
        -0.0000007421,
        -0.0000007484,
        -0.0000007421,
        -0.0000011111,
        -0.0000011320,
        -0.0000003077,
        -0.0000012285,
        -0.0000011825,
        -0.0000002460,
        -0.0000003077,
        -0.0000012285,
        -0.0000011825,
        -0.0000002460,
        -0.0000011111,
        -0.0000011320,
        -0.0000011111,
        -0.0000011320,
        -0.0000007484,
        -0.0000007421,
        -0.0000003077,
        -0.0000012285,
        -0.0000011825,
        -0.0000002460,
        -0.0000003077,
        -0.0000007484,
        -0.0000012285,
        -0.0000011825,
        -0.0000007421,
        -0.0000002460,
        -0.0000011111,
        -0.0000011320,
        -0.0000012324,
        -0.0000011859,
        -0.0000012324,
        -0.0000011859,
        -0.0000011111,
        -0.0000011320,
        -0.0000011111,
        -0.0000011320,
        -0.0000001028,
        -0.0000000528,
        -0.0000061043,
        -0.0000054852,
        -0.0000001028,
        -0.0000000528,
        -0.0000061043,
        -0.0000054852,
        -0.0000001028,
        -0.0000000528,
        -0.0000061043,
        -0.0000054852,
        -0.0000001028,
        -0.0000000528,
        -0.0000061043,
        -0.0000054852,
        -0.0000136044,
        -0.0000078949,
        -0.0000197836,
        -0.0000252407,
        -0.0000150705,
        -0.0000058395,
        -0.0000055968,
        -0.0000086793,
        -0.0000258495,
        -0.0000045423,
        -0.0000025932,
        -0.0000045867,
        -0.0000136044,
        -0.0000252407,
        -0.0000150705,
        -0.0000058395,
        -0.0000055968,
        -0.0000078949,
        -0.0000086793,
        -0.0000258495,
        -0.0000045423,
        -0.0000025932,
        -0.0000197836,
        -0.0000045867,
        -0.0000000567,
        -0.0000610837,
        -0.0000006147,
        -0.0000007795,
        -0.0000019007,
        -0.0000019619,
        -0.0000125843,
        -0.0000171451,
        -0.0000303908,
        -0.0000105724,
        -0.0000120418,
        -0.0000006147,
        -0.0000007795,
        -0.0000019007,
        -0.0000019619,
        -0.0000125843,
        -0.0000171451,
        -0.0000303908,
        -0.0000610837,
        -0.0000105724,
        -0.0000120418,
        -0.0000000567,
        -0.0000084900,
        -0.0000190384,
        -0.0000063671,
        -0.0000067133,
        -0.0000037436,
        -0.0000024719,
        -0.0000190384,
        -0.0000063671,
        -0.0000084900,
        -0.0000067133,
        -0.0000037436,
        -0.0000024719,
        -0.0002096102,
        -0.0000773039,
        -0.5600697020,
        -0.0000199676,
        -0.0000610911,
        -1.2613900419,
        -0.0023652036,
        -0.0000272110,
        -0.0002096102,
        -0.0000346158,
        -0.0000479764,
        -0.0000500261,
        -0.0005650258,
        -0.0034082519,
        -0.0000606923,
        -0.0025689291,
        -0.0002421192,
        -0.0002490056,
        -0.0003174064,
        -0.0000378659,
        -0.0000447885,
        -0.0002749821,
        -0.0000791729,
        -0.0000165798,
        -0.0001833732,
        -0.0000269193,
        -0.0000401152,
        -0.0005028912,
        -0.0000623181,
        -0.0000260829,
        -0.0000260829,
        -0.0000694916,
        -0.0002043327,
        -0.0000834514,
        -0.0000222680,
        -0.0000120055,
        -0.0000610911,
       -22.9853780414,
       -22.9853780414,
       -22.1567476341,
       -22.1567476341,
       -19.1065918217,
       -19.1065918217,
       -14.4474167447,
       -14.4474167447,
        -9.8824303418,
        -9.8824303418,
        -5.8462697009,
        -5.8462697009,
        -2.9658756029,
        -2.9658756029,
        -1.2404101115,
        -1.2404101115,
        -0.4202151511,
        -0.4202151511,
        -0.4198278901,
        -0.4202156232,
        -0.4197645374,
        -0.1089010630,
        -0.1089010630,
        -0.1092266384,
        -0.1089021272,
        -0.1090208263,
        -0.0208418730,
        -0.0208418730,
        -0.0215279996,
        -0.0208441904,
        -0.0212511776,
        -0.0025689291,
        -0.0032881155,
        -0.0025735857,
        -0.0030086230,
        -0.0002421192,
        -0.0008189558,
        -0.0005468296,
        -0.0021571094,
        -0.0003174064,
        -0.0000378659,
        -0.0006596125,
        -0.0041488022,
        -0.0000791729,
        -0.0005028912,
        -0.0029318678,
        -0.0000623181,
        -0.0063573213,
        -0.0000854280,
        -0.0494977855,
        -0.0000165798,
        -0.0000269193,
        -0.0002043327,
        -0.0000199676,
        -0.0004952328,
        -0.0000133146,
        -0.0000346158,
        -0.0000500261,
        -0.0000675967,
        -0.0077243539,
        -0.0000501945,
        -0.0000479764,
        -0.0000289821,
        -3.6003766093,
        -0.0005650258,
        -0.0000773039,
        -0.0000269544,
        -0.0004867574,
        -0.0000674222,
        -0.0415289527,
        -0.0000238220,
        -0.0004708701,
        -0.0013931545,
        -0.0002449291,
        -0.0011207016,
        -0.0023451236,
        -0.0003642815,
        -0.0000178538,
        -0.0188594708,
        -0.0668996625,
        -0.2025001447,
        -0.0069502465,
        -0.0257083764,
        -0.0012732225,
        -0.0415289527,
        -0.0045662702,
        -0.0273042020,
        -0.0001521951,
        -0.0000354026,
        -0.5623624376,
        -0.0000872246,
        -0.0000320405,
        -0.0003955841,
        -1.2646862314,
        -0.0001521951,
        -0.0000280605,
        -0.0000936675,
        -0.0000215742,
        -0.0000261126,
        -0.0000367835,
        -0.0003805789,
        -0.0001084618,
        -0.0000736550,
        -0.0000526867,
        -0.0000845860,
        -0.0006043240,
        -0.0000319017,
        -0.0000339708,
        -0.0000561909,
        -0.0002529402,
        -0.0000414908,
        -0.0000179837,
        -0.0000301290,
        -0.0000312218,
        -0.0003999500,
        -0.0000366961,
        -0.0000115020,
        -0.0000203492,
        -0.0000320405,
        -0.0005311651,
        -0.0003420026,
        -0.0025384627,
        -0.0000845860,
        -0.0003126054,
        -0.0000319017,
        -0.4222128804,
        -0.4217327482,
        -0.1101823670,
        -0.1098380443,
        -0.0215729696,
        -0.0213491600,
        -0.0029222335,
        -0.0027971328,
        -0.0003237755,
        -0.0002529402,
        -0.0000414908,
        -0.0000312218,
        -0.0000561909,
        -0.0000179837,
        -0.0000132105,
        -0.0000280605,
        -0.0000215742,
        -0.0000295220,
        -0.0006725059,
        -0.0000334864,
        -0.0000936675,
        -0.0000261126,
        -0.0000354026,
        -0.0000136044,
        -0.0000136044,
        -0.0000078949,
        -0.0000252407,
        -0.0000150705,
        -0.0000058395,
        -0.0000055968,
        -0.0000252407,
        -0.0000150705,
        -0.0000058395,
        -0.0000055968,
        -0.0000078949,
        -0.0000197836,
        -0.0000086793,
        -0.0000045423,
        -0.0000030744,
        -0.0000086793,
        -0.0000258495,
        -0.0000017149,
        -0.0000045423,
        -0.0000025932,
        -0.0000197836,
        -0.0000005760,
        -0.0000002209,
        -0.0000045867,
        -0.0000045867,
        -0.0000006147,
        -0.0000007795,
        -0.0000006147,
        -0.0000007795,
        -0.0000019007,
        -0.0000019619,
        -0.0000019007,
        -0.0000019619,
        -0.0000125843,
        -0.0000171451,
        -0.0000303908,
        -0.0000125843,
        -0.0000266387,
        -0.0000171451,
        -0.0000303908,
        -0.0000610837,
        -0.0000000567,
        -0.0000105724,
        -0.0000120418,
        -0.0000105724,
        -0.0000120418,
        -0.0000084900,
        -0.0000190384,
        -0.0000063671,
        -0.0000190384,
        -0.0000063671,
        -0.0000084900,
        -0.0000067133,
        -0.0000037436,
        -0.0000024719,
        -0.0000154343,
        -0.0000028552,
        -0.0000067133,
        -0.0000037436,
        -0.0000024719,
        -0.0002096102,
        -0.5600697020,
        -0.0000773039,
        -0.0002096102,
        -0.0000272110,
        -1.2613900419,
        -0.0023652036,
        -0.0000610911,
        -0.0004952328,
        -0.0000610911,
        -0.0034082519,
        -0.0000606923,
        -0.0002490056,
        -0.0000447885,
        -0.0002749821,
        -0.0001833732,
        -0.0000401152,
        -0.0000260829,
        -0.0000260829,
        -0.0000694916,
        -0.0000222680,
        -0.0000120055,
        -0.0000834514,
        -0.0008189558,
        -0.0003174064,
        -0.0006596125,
        -0.0000791729,
        -0.0000269193,
        -0.0005028912,
        -0.0000623181,
        -0.0063573213,
        -0.0000854280,
        -0.0002043327,
        -0.0000346158,
        -0.0000500261,
        -0.0000479764,
        -0.0005650258,
        -0.0003174064,
        -0.0000791729,
        -0.0000269193,
        -0.0005028912,
        -0.0000623181,
        -0.0063573213,
        -0.0002043327,
       -22.9851914609,
       -22.1560954130,
       -19.1051798924,
       -14.4449548107,
        -9.8792152204,
        -5.8428867415,
        -2.9631339897,
        -1.2389703616,
        -0.4202344691,
        -0.4198278901,
        -0.4202156232,
        -0.4197645374,
        -0.1100283854,
        -0.1092266384,
        -0.1089021272,
        -0.1090208263,
        -0.0226139957,
        -0.0215279996,
        -0.0208441904,
        -0.0212511776,
        -0.0045794120,
        -0.0032881155,
        -0.0025735857,
        -0.0030086230,
        -0.0040728041,
        -0.0005468296,
        -0.0021571094,
        -0.0144692948,
        -0.0041488022,
        -0.0029318678,
        -0.1787966836,
        -0.0494977855,
        -0.0000289821,
        -3.6003766093,
        -0.0000674222,
        -0.0000346158,
        -0.0000500261,
        -0.0000133146,
        -0.0000773039,
        -0.0000479764,
        -0.0005650258,
        -0.0000675967,
        -0.0077243539,
        -0.0000501945,
        -0.0000269544,
        -0.0045662702,
        -0.0000238220,
        -0.0004708701,
        -0.0002449291,
        -0.0023451236,
        -0.0003642815,
        -0.0000178538,
        -0.0069502465,
        -0.0012732225,
        -0.0069502465,
        -0.0012732225,
        -0.0188594708,
        -0.0668996625,
       -53.9245695704,
        -0.2025001447,
        -0.0257083764,
        -0.0273042020,
        -0.0658109204,
        -0.0658109204,
        -0.0001521951,
        -0.5623624376,
        -0.0000354026,
        -0.0001521951,
        -0.0003955841,
        -1.2646862314,
        -0.0000320405,
        -0.0000320405,
        -0.0000367835,
        -0.0000736550,
        -0.0000526867,
        -0.0000339708,
        -0.0000115020,
        -0.0000301290,
        -0.0000203492,
        -0.0003999500,
        -0.0000366961,
        -0.0000845860,
        -0.0000561909,
        -0.0002529402,
        -0.0000414908,
        -0.0000179837,
        -0.0000312218,
        -0.0000280605,
        -0.0000215742,
        -0.0000936675,
        -0.0000261126,
        -0.0000845860,
        -0.0000561909,
        -0.0002529402,
        -0.0000414908,
        -0.0000179837,
        -0.0000312218,
        -0.0005311651,
        -0.0003420026,
        -0.0025384627,
        -0.0003126054,
        -0.4217327482,
        -0.1098380443,
        -0.0213491600,
        -0.0027971328,
        -0.0003237755,
        -0.0000866269,
        -0.0000280605,
        -0.0000215742,
        -0.0000132105,
        -0.0000354026,
        -0.0000936675,
        -0.0000261126,
        -0.0000295220,
        -0.0006725059,
        -0.0000334864,
        -0.0000866269,
        -0.0000601629,
        -0.0000601629,
        -0.0190284666,
        -0.0043295738,
        -0.0035076419,
        -0.0035076419,
        -0.0037756308,
        -0.0003915322,
        -0.0026428628,
        -0.0004383959,
        -0.0004383959,
        -0.0005592056,
        -0.0000425682,
        -0.0002467329,
        -0.0823390813,
        -0.0832426187,
        -0.0184910504,
        -0.0184910504,
        -0.0187131056,
        -0.0033696544,
        -0.0002993184,
        -0.0002993184,
        -0.0000387861,
        -0.0000414576,
        -0.2153122170,
        -0.0000284522,
      -101.8682311901,
      -101.7188034615,
      -101.8439207859,
       -34.5586197475,
       -34.4440155048,
       -34.4947658320,
       -10.1040907159,
       -10.0270285032,
       -10.0288248159,
        -3.2759413483,
        -3.2261556084,
        -3.2055099275,
        -1.7856903833,
        -1.7515857843,
        -1.7244405374,
        -1.6902812007,
        -1.6480397829,
        -1.7797501728,
        -1.7655268490,
        -1.7797501728,
        -1.7655268490,
        -1.2098224766,
        -1.1885389194,
        -1.1670665213,
        -1.1482952919,
        -1.1147957920,
        -1.2066523217,
        -1.1954160393,
        -1.2066523217,
        -1.1954160393,
        -0.6592691915,
        -0.6476710096,
        -0.6344739477,
        -0.6150430094,
        -0.6696780640,
        -0.6659644991,
        -0.6696780640,
        -0.6659644991,
        -0.6694745757,
        -0.6651425415,
        -0.6534889703,
        -0.6396094370,
        -0.6534889703,
        -0.6396094370,
        -0.6448101394,
        -0.6281352680,
        -0.6659378623,
        -0.6579957529,
        -0.6659378623,
        -0.6579957529,
        -0.6695329435,
        -0.6695329435,
        -0.2780261012,
        -0.2741485010,
        -0.2687787509,
        -0.2621472127,
        -0.2819136749,
        -0.2806872176,
        -0.2819136749,
        -0.2806872176,
        -0.2819035482,
        -0.2804728728,
        -0.2758469078,
        -0.2711781016,
        -0.2758469078,
        -0.2711781016,
        -0.2726281621,
        -0.2669750219,
        -0.2805122776,
        -0.2778765993,
        -0.2805122776,
        -0.2778765993,
        -0.2819786538,
        -0.2819786538,
        -0.0848261103,
        -0.0845970373,
        -0.0828188351,
        -0.0856868353,
        -0.0856289075,
        -0.0856868353,
        -0.0856289075,
        -0.0857596519,
        -0.0856928347,
        -0.0843249907,
        -0.0840203656,
        -0.0843249907,
        -0.0840203656,
        -0.0836476639,
        -0.0853678503,
        -0.0852315737,
        -0.0853678503,
        -0.0852315737,
        -0.0858490651,
        -0.0858490651,
        -0.0178304282,
        -0.0183773257,
        -0.0180407191,
        -0.0177235096,
        -0.0179054404,
        -0.0177235096,
        -0.0179054404,
        -0.0178086416,
        -0.0180226150,
        -0.0178531116,
        -0.0178531116,
        -0.0179501470,
        -0.0177637463,
        -0.0181485179,
        -0.0177637463,
        -0.0181485179,
        -0.0178503344,
        -0.0178503344,
        -0.0030112250,
        -0.0028058300,
        -0.0029209236,
        -0.0028058300,
        -0.0029209236,
        -0.0028741680,
        -0.0030104899,
        -0.0030847381,
        -0.0030847381,
        -0.0032396609,
        -0.0028892943,
        -0.0031369852,
        -0.0028892943,
        -0.0031369852,
        -0.0028772961,
        -0.0028772961,
        -0.0002919284,
        -0.0002225047,
        -0.0002588204,
        -0.0002318745,
        -0.0002318745,
        -0.0002120687,
        -0.0002120687,
        -0.0001853417,
        -0.0002151286,
        -0.0001853417,
        -0.0002151286,
        -0.0000111145,
        -0.0000432002,
        -0.0000432002,
        -0.0000392306,
        -0.0034061460,
        -0.0003231832,
        -0.0003231832,
        -0.0000307668,
        -0.0005005967,
        -0.0003946116,
        -0.0000897626,
        -0.0000072667,
        -0.0000020590,
        -0.0000121986,
        -0.0000072667,
        -0.0000121986,
        -0.0000245273,
        -0.0000055783,
        -0.0000013126,
        -0.0000080353,
        -0.0000080353,
        -0.0000055783,
        -0.0000145732,
        -0.0000601629,
        -0.0035076419,
        -0.0003915322,
        -0.0004383959,
        -0.0006590832,
        -0.0003147918,
        -0.0000440715,
        -0.0000425682,
        -0.0002467329,
        -0.0033696544,
        -0.0002993184,
        -0.0000387861,
        -0.0000414576,
        -0.0190284666,
        -0.0187970872,
        -0.0043295738,
        -0.0035076419,
        -0.0037756308,
        -0.0039015157,
        -0.0026428628,
        -0.0004383959,
        -0.0005592056,
        -0.0006590832,
        -0.0003147918,
        -0.0823390813,
        -0.0832426187,
        -0.0829168959,
        -0.0184910504,
        -0.0184910504,
        -0.0187131056,
        -0.0031952407,
        -0.0002993184,
        -0.0000601629,
        -0.0000440715,
        -0.2153122170,
        -0.0000284522,
        -0.0856868353,
        -0.0856289075,
        -0.0177235096,
        -0.0179054404,
        -0.0177637463,
        -0.0181485179,
        -0.0178503344,
        -0.0030112250,
        -0.0028058300,
        -0.0029209236,
        -0.0028741680,
        -0.0030104899,
        -0.0030847381,
        -0.0028892943,
        -0.0031369852,
        -0.0028772961,
        -0.0002919284,
        -0.0001853417,
        -0.0002151286,
        -0.0002225047,
        -0.0002588204,
        -0.0002318745,
        -0.0002120687,
        -0.0002120687,
     -2370.2437224265,
     -2318.9568761059,
     -2109.4414583515,
     -1699.2816917969,
     -1271.6056544513,
      -843.2437165956,
      -493.9923376821,
      -245.3261773360,
      -101.8682290295,
      -101.6613374015,
      -101.8145455041,
      -101.8682311901,
      -101.7188034615,
      -101.8439207859,
       -34.5586150271,
       -34.3937436765,
       -34.4558724349,
       -34.5586197475,
       -34.4440155048,
       -34.4947658320,
       -10.1040816514,
        -9.9888073681,
        -9.9910097071,
       -10.1040907159,
       -10.0270285032,
       -10.0288248159,
        -3.2759231370,
        -3.1984940294,
        -3.1733023811,
        -3.2759413483,
        -3.2261556084,
        -3.2055099275,
        -1.7856496568,
        -1.7308532503,
        -1.6976837373,
        -1.7856903833,
        -1.7515857843,
        -1.7244405374,
        -1.6902812007,
        -1.6480397829,
        -1.7797501728,
        -1.7655268490,
        -1.7797501728,
        -1.7655268490,
        -1.2097501727,
        -1.1749717759,
        -1.1487103402,
        -1.2098224766,
        -1.1885389194,
        -1.1670665213,
        -1.1482952919,
        -1.1147957920,
        -1.2066523217,
        -1.1954160393,
        -1.2066523217,
        -1.1954160393,
        -0.6694580534,
        -0.6592691915,
        -0.6476710096,
        -0.6344739477,
        -0.6150430094,
        -0.6696780640,
        -0.6659644991,
        -0.6696780640,
        -0.6659644991,
        -0.6694745757,
        -0.6651425415,
        -0.6534889703,
        -0.6396094370,
        -0.6534889703,
        -0.6396094370,
        -0.6448101394,
        -0.6281352680,
        -0.6659378623,
        -0.6579957529,
        -0.6659378623,
        -0.6579957529,
        -0.6695329435,
        -0.6695329435,
        -0.6412621623,
        -0.6235885279,
        -0.6632481672,
        -0.6536426044,
        -0.2818841275,
        -0.2780261012,
        -0.2741485010,
        -0.2687787509,
        -0.2621472127,
        -0.2819136749,
        -0.2806872176,
        -0.2819136749,
        -0.2806872176,
        -0.2819035482,
        -0.2804728728,
        -0.2758469078,
        -0.2711781016,
        -0.2758469078,
        -0.2711781016,
        -0.2726281621,
        -0.2669750219,
        -0.2805122776,
        -0.2778765993,
        -0.2805122776,
        -0.2778765993,
        -0.2819786538,
        -0.2819786538,
        -0.2713017028,
        -0.2652935754,
        -0.2794655804,
        -0.2762672869,
        -0.0857511926,
        -0.0848261103,
        -0.0845970373,
        -0.0828188351,
        -0.0856868353,
        -0.0856289075,
        -0.0857596519,
        -0.0856928347,
        -0.0843249907,
        -0.0840203656,
        -0.0843249907,
        -0.0840203656,
        -0.0836476639,
        -0.0853678503,
        -0.0852315737,
        -0.0853678503,
        -0.0852315737,
        -0.0858490651,
        -0.0858490651,
        -0.0833560446,
        -0.0850881056,
        -0.0849113945,
        -0.0177692853,
        -0.0178304282,
        -0.0183773257,
        -0.0180407191,
        -0.0177235096,
        -0.0179054404,
        -0.0178086416,
        -0.0180226150,
        -0.0178531116,
        -0.0178531116,
        -0.0179501470,
        -0.0177637463,
        -0.0181485179,
        -0.0178503344,
        -0.0179750405,
        -0.0177372157,
        -0.0181951618,
        -0.0028178447,
        -0.0028058300,
        -0.0029209236,
        -0.0030847381,
        -0.0032396609,
        -0.0028892943,
        -0.0031369852,
        -0.0028772961,
        -0.0032906162,
        -0.0028992459,
        -0.0002318745,
        -0.0002343425,
        -0.0000111145,
        -0.0001805725,
        -0.0001853417,
        -0.0002151286,
        -0.0000432002,
        -0.0000392306,
        -0.0003231832,
        -0.0000307668,
        -0.0000897626,
        -0.0034061460,
        -0.0003231832,
        -0.0005005967,
        -0.0003946116,
        -0.0004173984,
        -0.0000432002,
        -0.0000353048,
        -0.0000072667,
        -0.0000121986,
        -0.0000072667,
        -0.0000121986,
        -0.0000055783,
        -0.0000080353,
        -0.0000080353,
        -0.0000055783,
        -0.0000145732,
        -0.0000145732,
        -0.0000098027,
        -0.0000098027,
        -0.0000098027,
        -0.0000098027,
        -0.0000035140,
        -0.0000035140,
        -0.0000007445,
        -0.0000007445,
        -0.0000035140,
        -0.0000007445,
        -0.0000007445,
        -0.0000035140,
        -0.0000000724,
        -0.0000088389,
        -0.0000000725,
        -0.0000000724,
        -0.0000088389,
        -0.0000000725,
        -0.0000000724,
        -0.0000088389,
        -0.0000000725,
        -0.0000000724,
        -0.0000088389,
        -0.0000000725,
        -0.0000072667,
        -0.0000136509,
        -0.0000121986,
        -0.0000020590,
        -0.0000072667,
        -0.0000121986,
        -0.0000020590,
        -0.0000136509,
        -0.0000055783,
        -0.0000245273,
        -0.0000080353,
        -0.0000013126,
        -0.0000080353,
        -0.0000055783,
        -0.0000013126,
        -0.0000245273,
        -0.0000034405,
        -0.0000034405,
        -0.0000601629,
        -0.0035076419,
        -0.0003915322,
        -0.0004383959,
        -0.0000425682,
        -0.0002467329,
        -0.0033696544,
        -0.0002993184,
        -0.0000387861,
        -0.0000414576,
        -0.0190284666,
        -0.0187970872,
        -0.0043295738,
        -0.0035076419,
        -0.0037756308,
        -0.0039015157,
        -0.0026428628,
        -0.0004383959,
        -0.0005592056,
        -0.0006590832,
        -0.0003147918,
        -0.0823390813,
        -0.0832426187,
        -0.0829168959,
        -0.0184910504,
        -0.0184910504,
        -0.0187131056,
        -0.0031952407,
        -0.0002993184,
        -0.0000440715,
        -0.0000601629,
        -0.0000284522,
        -0.0856868353,
        -0.0856289075,
        -0.0177235096,
        -0.0179054404,
        -0.0177637463,
        -0.0181485179,
        -0.0178503344,
        -0.0030112250,
        -0.0028058300,
        -0.0029209236,
        -0.0028741680,
        -0.0030104899,
        -0.0030847381,
        -0.0028892943,
        -0.0031369852,
        -0.0028772961,
        -0.0002919284,
        -0.0001853417,
        -0.0002151286,
        -0.0002225047,
        -0.0002588204,
        -0.0002318745,
        -0.0002120687,
        -0.0002120687,
        -0.0000111145,
     -2370.2437224265,
     -2318.9568761059,
     -2109.4414583515,
     -1699.2816917969,
     -1271.6056544513,
      -843.2437165956,
      -493.9923376821,
      -245.3261773360,
      -101.8682290296,
      -101.6613374015,
      -101.8145455041,
      -101.8682311901,
      -101.7188034615,
      -101.8439207859,
       -34.5586150271,
       -34.3937436765,
       -34.4558724349,
       -34.5586197475,
       -34.4440155048,
       -34.4947658320,
       -10.1040816514,
        -9.9888073681,
        -9.9910097071,
       -10.1040907159,
       -10.0270285032,
       -10.0288248159,
        -3.2759231370,
        -3.1984940294,
        -3.1733023811,
        -3.2759413483,
        -3.2261556084,
        -3.2055099275,
        -1.7856496568,
        -1.7308532503,
        -1.6976837373,
        -1.7856903833,
        -1.7515857843,
        -1.7244405374,
        -1.6902812007,
        -1.6480397829,
        -1.7797501728,
        -1.7655268490,
        -1.7797501728,
        -1.7655268490,
        -1.2097501727,
        -1.1749717759,
        -1.1487103402,
        -1.2098224766,
        -1.1885389194,
        -1.1670665213,
        -1.1482952919,
        -1.1147957920,
        -1.2066523217,
        -1.1954160393,
        -1.2066523217,
        -1.1954160393,
        -0.6694580534,
        -0.6592691915,
        -0.6476710096,
        -0.6344739477,
        -0.6150430094,
        -0.6696780640,
        -0.6659644991,
        -0.6696780640,
        -0.6659644991,
        -0.6694745757,
        -0.6651425415,
        -0.6534889703,
        -0.6396094370,
        -0.6534889703,
        -0.6396094370,
        -0.6448101394,
        -0.6281352680,
        -0.6659378623,
        -0.6579957529,
        -0.6659378623,
        -0.6579957529,
        -0.6695329435,
        -0.6695329435,
        -0.6412621623,
        -0.6235885279,
        -0.6632481672,
        -0.6536426044,
        -0.2818841275,
        -0.2780261012,
        -0.2741485010,
        -0.2687787509,
        -0.2621472127,
        -0.2819136749,
        -0.2806872176,
        -0.2819136749,
        -0.2806872176,
        -0.2819035482,
        -0.2804728728,
        -0.2758469078,
        -0.2711781016,
        -0.2758469078,
        -0.2711781016,
        -0.2726281621,
        -0.2669750219,
        -0.2805122776,
        -0.2778765993,
        -0.2805122776,
        -0.2778765993,
        -0.2819786538,
        -0.2819786538,
        -0.2713017028,
        -0.2652935754,
        -0.2794655804,
        -0.2762672869,
        -0.0857511926,
        -0.0848261103,
        -0.0845970373,
        -0.0828188351,
        -0.0856868353,
        -0.0856289075,
        -0.0857596519,
        -0.0856928347,
        -0.0843249907,
        -0.0840203656,
        -0.0843249907,
        -0.0840203656,
        -0.0836476639,
        -0.0853678503,
        -0.0852315737,
        -0.0853678503,
        -0.0852315737,
        -0.0858490651,
        -0.0858490651,
        -0.0833560446,
        -0.0850881056,
        -0.0849113945,
        -0.0177692853,
        -0.0178304282,
        -0.0183773257,
        -0.0180407191,
        -0.0177235096,
        -0.0179054404,
        -0.0178086416,
        -0.0180226150,
        -0.0178531116,
        -0.0178531116,
        -0.0179501470,
        -0.0177637463,
        -0.0181485179,
        -0.0178503344,
        -0.0179750405,
        -0.0177372157,
        -0.0181951618,
        -0.0028178447,
        -0.0028058300,
        -0.0029209236,
        -0.0030847381,
        -0.0032396609,
        -0.0028892943,
        -0.0031369852,
        -0.0028772961,
        -0.0032906162,
        -0.0028992459,
        -0.0002318745,
        -0.0002343425,
        -0.0001805725,
        -0.0001853417,
        -0.0002151286,
        -0.0000432002,
        -0.0000392306,
        -0.0003231832,
        -0.0000307668,
        -0.0000897626,
        -0.0034061460,
        -0.0003231832,
        -0.0005005967,
        -0.0003946116,
        -0.0004173984,
        -0.0000353048,
        -0.0000432002,
        -0.0000072667,
        -0.0000121986,
        -0.0000072667,
        -0.0000121986,
        -0.0000136509,
        -0.0000020590,
        -0.0000136509,
        -0.0000055783,
        -0.0000080353,
        -0.0000080353,
        -0.0000055783,
        -0.0000013126,
        -0.0000245273,
        -0.0000034405,
        -0.0000145732,
        -0.0000034405,
        -0.0003915322,
        -0.0000425682,
        -0.0002467329,
        -0.0033696544,
        -0.0000387861,
        -0.0000414576,
        -0.0035076419,
        -0.0004383959,
        -0.0006590832,
        -0.0003147918,
        -0.0000440715,
        -0.0031952407,
        -0.0002993184,
        -0.0000601629,
        -0.0035076419,
        -0.0004383959,
        -0.0031952407,
        -0.0002993184,
        -0.0191600949,
        -0.0190284666,
        -0.0187970872,
        -0.0187970872,
        -0.0046823621,
        -0.0043295738,
        -0.0037756308,
        -0.0039015157,
        -0.0039015157,
        -0.0158501052,
        -0.0026428628,
        -0.0005592056,
        -0.0821008297,
        -0.0823390813,
        -0.0832426187,
        -0.0829168959,
        -0.0829168959,
        -0.0184910504,
        -0.0184910504,
        -0.0187131056,
        -0.3210413618,
        -0.2002854694,
        -0.0000601629,
        -0.0000284522,
        -0.0030112250,
        -0.0028741680,
        -0.0030104899,
        -0.0002919284,
        -0.0002225047,
        -0.0002588204,
        -0.0002120687,
        -0.0002120687,
        -0.0000111145,
        -0.0857511926,
        -0.0856868353,
        -0.0856289075,
        -0.0177692853,
        -0.0177235096,
        -0.0179054404,
        -0.0177637463,
        -0.0181485179,
        -0.0178503344,
        -0.0177372157,
        -0.0181951618,
        -0.0028178447,
        -0.0028058300,
        -0.0029209236,
        -0.0030847381,
        -0.0028892943,
        -0.0031369852,
        -0.0028772961,
        -0.0028992459,
        -0.0001805725,
        -0.0001853417,
        -0.0002151286,
        -0.0002318745,
        -0.0002343425,
        -0.0857511926,
        -0.0856868353,
        -0.0856289075,
        -0.0177692853,
        -0.0177235096,
        -0.0179054404,
        -0.0177637463,
        -0.0181485179,
        -0.0178503344,
        -0.0177372157,
        -0.0181951618,
        -0.0028178447,
        -0.0028058300,
        -0.0029209236,
        -0.0030847381,
        -0.0028892943,
        -0.0031369852,
        -0.0028772961,
        -0.0028992459,
        -0.0001805725,
        -0.0001853417,
        -0.0002151286,
        -0.0002318745,
        -0.0002343425,
     -2370.2437224265,
     -2370.2437224265,
     -2370.2348246590,
     -2370.2525395668,
     -2318.9568761059,
     -2318.9568761059,
     -2318.9274602068,
     -2318.9855047820,
     -2109.4414583515,
     -2109.4414583515,
     -2109.3734646567,
     -2109.5049490917,
     -1699.2816917969,
     -1699.2816917969,
     -1699.1460322368,
     -1699.3948806960,
     -1271.6056544513,
     -1271.6056544513,
     -1271.3925853878,
     -1271.7529427212,
      -843.2437165956,
      -843.2437165956,
      -842.9559408097,
      -843.3841905543,
      -493.9923376821,
      -493.9923376821,
      -493.6430657634,
      -494.0706581054,
      -245.3261773360,
      -245.3261773360,
      -244.9466848341,
      -245.2930026478,
      -101.8682290296,
      -101.8682290296,
      -101.4994655400,
      -101.7160127925,
      -101.6613374015,
      -101.8145455041,
      -101.6613374015,
      -101.8145455041,
      -101.8682311901,
      -101.7188034615,
      -101.8439207859,
       -34.5586150271,
       -34.5586150271,
       -34.2473231912,
       -34.3350700485,
       -34.3937436765,
       -34.4558724349,
       -34.3937436765,
       -34.4558724349,
       -34.5586197475,
       -34.4440155048,
       -34.4947658320,
       -10.1040816514,
       -10.1040816514,
        -9.8745630005,
        -9.8776845265,
        -9.9888073681,
        -9.9910097071,
        -9.9888073681,
        -9.9910097071,
       -10.1040907159,
       -10.0270285032,
       -10.0288248159,
        -3.2759231370,
        -3.2759231370,
        -3.1142219520,
        -3.0789981988,
        -3.1984940294,
        -3.1733023811,
        -3.1984940294,
        -3.1733023811,
        -3.2759413483,
        -3.2261556084,
        -3.2055099275,
        -1.7856496568,
        -1.7856496568,
        -1.6667822728,
        -1.6202035607,
        -1.7308532503,
        -1.6976837373,
        -1.7308532503,
        -1.6976837373,
        -1.7856903833,
        -1.7515857843,
        -1.7244405374,
        -1.6902812007,
        -1.6480397829,
        -1.7797501728,
        -1.7655268490,
        -1.7797501728,
        -1.7655268490,
        -1.2097501727,
        -1.2097501727,
        -1.1327840353,
        -1.0958114949,
        -1.1749717759,
        -1.1487103402,
        -1.1749717759,
        -1.1487103402,
        -1.2098224766,
        -1.1885389194,
        -1.1670665213,
        -1.1482952919,
        -1.1147957920,
        -1.2066523217,
        -1.1954160393,
        -1.2066523217,
        -1.1954160393,
        -0.6694580534,
        -0.6694580534,
        -0.6315509302,
        -0.6114162104,
        -0.6592691915,
        -0.6476710096,
        -0.6344739477,
        -0.6150430094,
        -0.6696780640,
        -0.6659644991,
        -0.6696780640,
        -0.6659644991,
        -0.6694745757,
        -0.6651425415,
        -0.6534889703,
        -0.6396094370,
        -0.6534889703,
        -0.6396094370,
        -0.6448101394,
        -0.6281352680,
        -0.6659378623,
        -0.6579957529,
        -0.6659378623,
        -0.6579957529,
        -0.6695329435,
        -0.6695329435,
        -0.6412621623,
        -0.6235885279,
        -0.6632481672,
        -0.6536426044,
        -0.6412621623,
        -0.6235885279,
        -0.6632481672,
        -0.6536426044,
        -0.2818841275,
        -0.2818841275,
        -0.2676900870,
        -0.2608107971,
        -0.2780261012,
        -0.2741485010,
        -0.2687787509,
        -0.2621472127,
        -0.2819136749,
        -0.2806872176,
        -0.2819136749,
        -0.2806872176,
        -0.2819035482,
        -0.2804728728,
        -0.2758469078,
        -0.2711781016,
        -0.2758469078,
        -0.2711781016,
        -0.2726281621,
        -0.2669750219,
        -0.2805122776,
        -0.2778765993,
        -0.2805122776,
        -0.2778765993,
        -0.2819786538,
        -0.2819786538,
        -0.2713017028,
        -0.2652935754,
        -0.2794655804,
        -0.2762672869,
        -0.2713017028,
        -0.2652935754,
        -0.2794655804,
        -0.2762672869,
        -0.0825842268,
        -0.0848261103,
        -0.0845970373,
        -0.0828188351,
        -0.0857596519,
        -0.0856928347,
        -0.0843249907,
        -0.0840203656,
        -0.0843249907,
        -0.0840203656,
        -0.0836476639,
        -0.0853678503,
        -0.0852315737,
        -0.0853678503,
        -0.0852315737,
        -0.0858490651,
        -0.0858490651,
        -0.0833560446,
        -0.0850881056,
        -0.0849113945,
        -0.0833560446,
        -0.0850881056,
        -0.0849113945,
        -0.0180658908,
        -0.0178304282,
        -0.0183773257,
        -0.0180407191,
        -0.0178086416,
        -0.0180226150,
        -0.0178531116,
        -0.0178531116,
        -0.0179501470,
        -0.0179750405,
        -0.0179750405,
        -0.0032396609,
        -0.0032906162,
        -0.0032906162,
        -0.2153122170,
        -0.2153122170,
        -0.0000392306,
        -0.0000307668,
        -0.0000897626,
        -0.0003231832,
        -0.0004173984,
        -0.0000353048,
        -0.0000432002,
        -0.0003231832,
        -0.0004173984,
        -0.0000353048,
        -0.0034578015,
        -0.0034061460,
        -0.0005884554,
        -0.0005005967,
        -0.0003946116,
        -0.0218876453,
        -0.0000432002,
    };

    // Reference

    const double reference_gradient[n_atom * 3]
    {
        0.0, 0.0,  0.2783879717,
        0.0, 0.0, -0.2784345708,
    };

    for (int i_atom = 0; i_atom < n_atom; i_atom++)
        printf("Reference gradient  %2d  %15.10f  %15.10f  %15.10f\n",
            i_atom, reference_gradient[i_atom * 3 + 0], reference_gradient[i_atom * 3 + 1], reference_gradient[i_atom * 3 + 2]);

    // Analytical

    double* analytical_gradient = new double[n_atom * 3];
    memset(analytical_gradient, 0, n_atom * 3 * sizeof(double));

    PeriodicBox::BeckeWeightGradient(n_grid_point, grid_points, n_atom, atom_xyz, atom_radius, periodic_parameter.lattice.unit_cell, epsilon_xc, analytical_gradient);

    for (int i_atom = 0; i_atom < n_atom; i_atom++)
        printf("Analytical gradient %2d  %15.10f  %15.10f  %15.10f\n",
            i_atom, analytical_gradient[i_atom * 3 + 0], analytical_gradient[i_atom * 3 + 1], analytical_gradient[i_atom * 3 + 2]);

    delete[] analytical_gradient;

    // Numerical

    double* numerical_gradient  = new double[n_atom * 3];
    memset(numerical_gradient,  0, n_atom * 3 * sizeof(double));

    double* atom_xyz_copy = new double[n_atom * 3];
    memcpy(atom_xyz_copy, atom_xyz, n_atom * 3 * sizeof(double));
    PeriodicBox::GridPoint* grid_points_copy = new PeriodicBox::GridPoint[n_grid_point];
    memcpy(grid_points_copy, grid_points, n_grid_point * sizeof(PeriodicBox::GridPoint));

    const double dx = 1e-4;

    for (int i_atom = 0; i_atom < n_atom; i_atom++)
        for (int i_xyz = 0; i_xyz < 3; i_xyz++) {
            atom_xyz_copy[i_atom * 3 + i_xyz] = atom_xyz[i_atom * 3 + i_xyz] + dx;
            for (int i_point = 0; i_point < n_grid_point; i_point++) {
                grid_points_copy[i_point].w_total = NAN;
                if (grid_points[i_point].i_atom == i_atom) {
                    grid_points_copy[i_point].x = grid_points[i_point].x + ((i_xyz == 0) ? dx : 0);
                    grid_points_copy[i_point].y = grid_points[i_point].y + ((i_xyz == 1) ? dx : 0);
                    grid_points_copy[i_point].z = grid_points[i_point].z + ((i_xyz == 2) ? dx : 0);
                }
            }
            BeckeWeights(n_grid_point, grid_points_copy, n_atom, atom_xyz_copy, atom_radius, periodic_parameter.lattice.unit_cell);
            const double E_plus = becke_energy(n_grid_point, grid_points_copy, epsilon_xc);

            atom_xyz_copy[i_atom * 3 + i_xyz] = atom_xyz[i_atom * 3 + i_xyz] - dx;
            for (int i_point = 0; i_point < n_grid_point; i_point++) {
                grid_points_copy[i_point].w_total = NAN;
                if (grid_points[i_point].i_atom == i_atom) {
                    grid_points_copy[i_point].x = grid_points[i_point].x - ((i_xyz == 0) ? dx : 0);
                    grid_points_copy[i_point].y = grid_points[i_point].y - ((i_xyz == 1) ? dx : 0);
                    grid_points_copy[i_point].z = grid_points[i_point].z - ((i_xyz == 2) ? dx : 0);
                }
            }
            BeckeWeights(n_grid_point, grid_points_copy, n_atom, atom_xyz_copy, atom_radius, periodic_parameter.lattice.unit_cell);
            const double E_minus = becke_energy(n_grid_point, grid_points_copy, epsilon_xc);

            numerical_gradient[i_atom * 3 + i_xyz] = (E_plus - E_minus) / dx / 2.0;

            atom_xyz_copy[i_atom * 3 + i_xyz] = atom_xyz[i_atom * 3 + i_xyz];
            for (int i_point = 0; i_point < n_grid_point; i_point++) {
                if (grid_points[i_point].i_atom == i_atom) {
                    grid_points_copy[i_point].x = grid_points[i_point].x;
                    grid_points_copy[i_point].y = grid_points[i_point].y;
                    grid_points_copy[i_point].z = grid_points[i_point].z;
                }
            }
        }


    for (int i_atom = 0; i_atom < n_atom; i_atom++)
        printf("Numerical gradient  %2d  %15.10f  %15.10f  %15.10f\n",
            i_atom, numerical_gradient[i_atom * 3 + 0], numerical_gradient[i_atom * 3 + 1], numerical_gradient[i_atom * 3 + 2]);

    double max_abs_diff = 0.0;
    for (int i_atom = 0; i_atom < n_atom; i_atom++)
        for (int i_xyz = 0; i_xyz < 3; i_xyz++) {
            const double abs_diff = fabs(analytical_gradient[i_atom * 3 + i_xyz] - numerical_gradient[i_atom * 3 + i_xyz]);
            if (abs_diff > max_abs_diff) max_abs_diff = abs_diff;
        }
    printf("\nMax abs diff between analytical and numerical: %.5e\n", max_abs_diff);

    delete[] atom_xyz_copy;
    delete[] grid_points_copy;
    delete[] numerical_gradient;

    return 0;
}
