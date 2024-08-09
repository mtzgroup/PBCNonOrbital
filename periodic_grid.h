#ifndef PERIODIC_GRID_H_
#define PERIODIC_GRID_H_

#include "periodic_types.h"

namespace PeriodicBox
{
    struct GridPoint
    {
        double x;
        double y;
        double z;
        double w_total; // w_Becke * w_fixed
        double w_fixed;
        int i_atom;
    };

    class PeriodicGrid
    {
    private:
        enum PeriodicGridMode { Becke, Uniform } grid_mode;
        int n_point;
        GridPoint* points;

        int n_box;
        int* box_start;
        int* box_stop;

        // record info for becke weight gradient
        int n_atom;
        double* atom_xyz;
        double* atom_radius;

        static void place_point_into_box(const int n_point, GridPoint* points, int* n_box, int* box_start, int* box_stop, const LatticeVector unit_cell);

    public:
        PeriodicGrid();
        ~PeriodicGrid();

        void init_uniform_grid(
            const int n_grid_each_dimension[3],
            const LatticeInfo& lattice
        );

        void init_atom_centered_grid(
            const int n_point,
            const double* point_xyz,
            const double* point_fixed_weight,
            const int* point_i_atom,
            const int n_atom,
            const double* atom_xyz,
            const double* atom_radius,
            const LatticeInfo& lattice
        );

        void reorderPointsForGradientCalculation();

        int get_n_point() const { return n_point; }
        const GridPoint& get_point(const int i_point) const { return points[i_point]; }

        int get_n_box() const { return n_box; }
        int get_box_start(const int i_box) const { return box_start[i_box]; }
        int  get_box_stop(const int i_box) const { return  box_stop[i_box]; }

        // // Call this only after a reorderPointsForGradientCalculation() function call!
        // int getBoxCntr(const int box) const { return points[box_start[box]].cntr; }

        // int numAtoms() const { return num_atoms; }
        // const double* AtmXYZ() const { return atmXYZ; }

        // void ComputeBeckeGrad(const double* E, double* Egrad) const;
    };

}

#endif
