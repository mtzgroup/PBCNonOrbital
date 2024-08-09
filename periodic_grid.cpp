#include "periodic_grid.h"

#include "periodic_lattice_export.h"
#include "periodic_lattice_internal.h"
#include "periodic_becke/periodic_becke.h"

#include <distance_cutoff.h> // From gridbox, for BOX_LENGTH
#include <utils.h> // From intbox, for DIE()

#include <string.h> // For memcpy()
#include <math.h>

namespace PeriodicBox
{
    PeriodicGrid::PeriodicGrid()
    {
        grid_mode = PeriodicGridMode::Becke;
        n_point = 0;
        points = nullptr;
        n_box = 0;
        box_start = nullptr;
        box_stop  = nullptr;
        n_atom = 0;
        atom_xyz    = nullptr;
        atom_radius = nullptr;
    }

    PeriodicGrid::~PeriodicGrid()
    {
        if (points != nullptr)
            free(points);
        if (box_start != nullptr)
            free(box_start);
        if (box_stop != nullptr)
            free(box_stop);
        if (atom_xyz != nullptr)
            free(atom_xyz);
        if (atom_radius != nullptr)
            free(atom_radius);
    }

    void PeriodicGrid::init_uniform_grid(const int n_grid_each_dimension[3], const LatticeInfo& lattice)
    {
        if (this->points != nullptr || this->box_start != nullptr || this->box_stop != nullptr)
            DIE("Don't call PeriodicBox::Grid::init_grid() more than once!");

        grid_mode = PeriodicGridMode::Uniform;

        const int n_point = n_grid_each_dimension[0] * n_grid_each_dimension[1] * n_grid_each_dimension[2];
        this->n_point = n_point;
        this->points = (GridPoint*)malloc(sizeof(GridPoint) * n_point);
        for (int i_x = 0; i_x < n_grid_each_dimension[0]; i_x++)
            for (int i_y = 0; i_y < n_grid_each_dimension[1]; i_y++)
                for (int i_z = 0; i_z < n_grid_each_dimension[2]; i_z++) {
                    const int p = i_x * n_grid_each_dimension[1] * n_grid_each_dimension[2] + i_y * n_grid_each_dimension[2] + i_z;
                    const double relative_coordinate[3] { (double)i_x / (double)n_grid_each_dimension[0], (double)i_y / (double)n_grid_each_dimension[1], (double)i_z / (double)n_grid_each_dimension[2] };
                    double absolute_coordinate[3];
                    get_absolute_coord_from_fractional_coord(absolute_coordinate, relative_coordinate, lattice.unit_cell.real_space);
                    this->points[p].x = absolute_coordinate[0];
                    this->points[p].y = absolute_coordinate[1];
                    this->points[p].z = absolute_coordinate[2];
                    this->points[p].w_fixed = lattice.V_real / n_point;
                    this->points[p].i_atom = -1;
                    this->points[p].w_total = lattice.V_real / n_point;
                }

        int* box_start_overallocated = (int*)malloc(sizeof(int) * n_point);
        int* box_stop_overallocated  = (int*)malloc(sizeof(int) * n_point);

        n_box = -1;
        place_point_into_box(n_point, points, &n_box, box_start_overallocated, box_stop_overallocated, lattice.unit_cell);

        box_start = (int*)malloc(sizeof(int) * n_box);
        box_stop  = (int*)malloc(sizeof(int) * n_box);
        memcpy(box_start, box_start_overallocated, sizeof(int) * n_box);
        memcpy(box_stop,  box_stop_overallocated,  sizeof(int) * n_box);
        free(box_start_overallocated);
        free(box_stop_overallocated);
    }

    void PeriodicGrid::init_atom_centered_grid(
        const int n_point,
        const double* point_xyz,
        const double* point_fixed_weight,
        const int* point_i_atom,
        const int n_atom,
        const double* atom_xyz,
        const double* atom_radius,
        const LatticeInfo& lattice
    )
    {
        if (this->points != nullptr || this->box_start != nullptr || this->box_stop != nullptr)
            DIE("Don't call PeriodicBox::Grid::init_grid() more than once!");

        grid_mode = PeriodicGridMode::Becke;

        this->n_point = n_point;
        this->points = (GridPoint*)malloc(sizeof(GridPoint) * n_point);
        for (int p = 0; p < n_point; p++) {
            this->points[p].x = point_xyz[3 * p + 0];
            this->points[p].y = point_xyz[3 * p + 1];
            this->points[p].z = point_xyz[3 * p + 2];
            this->points[p].w_fixed = point_fixed_weight[p];
            this->points[p].i_atom = point_i_atom[p];
            this->points[p].w_total = NAN;
        }

        int* box_start_overallocated = (int*)malloc(sizeof(int) * n_point);
        int* box_stop_overallocated  = (int*)malloc(sizeof(int) * n_point);

        this->n_box = -1;
        place_point_into_box(n_point, points, &(this->n_box), box_start_overallocated, box_stop_overallocated, lattice.unit_cell);

        box_start = (int*)malloc(sizeof(int) * this->n_box);
        box_stop  = (int*)malloc(sizeof(int) * this->n_box);
        memcpy(box_start, box_start_overallocated, sizeof(int) * this->n_box);
        memcpy(box_stop,  box_stop_overallocated,  sizeof(int) * this->n_box);
        free(box_start_overallocated);
        free(box_stop_overallocated);

        this->n_atom = n_atom;
        this->atom_xyz = (double*)malloc(sizeof(double) * n_atom * 3);
        memcpy(this->atom_xyz, atom_xyz, sizeof(double) * n_atom * 3);

        this->atom_radius = (double*)malloc(sizeof(double) * n_atom);
        memcpy(this->atom_radius, atom_radius, sizeof(double) * n_atom);

        BeckeWeights(n_point, this->points, n_atom, atom_xyz, atom_radius, lattice.unit_cell);

        double sum_w = 0.0;
        for (int i = 0; i < n_point; i++)
            sum_w += this->points[i].w_total;
        printf("Total weight summation is %.5f, which suppose to be the same as real space unit cell volume %.5f, the difference indicates how bad your result is.\n",
            sum_w, lattice.V_real); fflush(stdout);
    }

    void PeriodicGrid::reorderPointsForGradientCalculation()
    {
        if (points == nullptr || box_start == nullptr || box_stop == nullptr)
            DIE("Please call PeriodicBox::Grid::reorderPointsForGradientCalculation() after calling initMolecularGrid()!");

        if (grid_mode == PeriodicGridMode::Uniform)
            return;

        DIE("PeriodicGrid::reorderPointsForGradientCalculation() not implemented!\n");
    //     std::stable_sort(points, points + num_points,
    //         [](const GridPoint& p1, const GridPoint& p2) { return p1.cntr < p2.cntr; });

    //     free(box_start);
    //     free(box_stop);
    //     num_boxes = 0;
    //     int* temp_box_start = (int*)malloc(sizeof(int) * num_points);
    //     int* temp_box_stop  = (int*)malloc(sizeof(int) * num_points);
    //     int i_point = 0;
    //     for (int i_atom = 0; i_atom < num_atoms; i_atom++) {
    //         int atmPtCnt = 0;
    //         const int atmPtStart = i_point;
    //         for (int p = 0; p < num_points; p++) {
    //             if (points[p].cntr == i_atom) {
    //                 i_point++;
    //                 atmPtCnt++;
    //             }
    //         }

    //         if (atmPtCnt==0)
    //             continue;

    //         // create box for atomic grid
    //         int atmBoxes = -1;
    //         place_point_into_box(atmPtCnt, points + atmPtStart, &atmBoxes, temp_box_start + num_boxes, temp_box_stop + num_boxes);
    //         for (int i = 0; i < atmBoxes; i++) {
    //             temp_box_start[num_boxes + i] += atmPtStart;
    //             temp_box_stop [num_boxes + i] += atmPtStart;
    //         }
    //         num_boxes += atmBoxes;
    //     }

    //     box_start = (int*)malloc(sizeof(int) * num_boxes);
    //     box_stop  = (int*)malloc(sizeof(int) * num_boxes);
    //     memcpy(box_start, temp_box_start, sizeof(int) * num_boxes);
    //     memcpy(box_stop,  temp_box_stop,  sizeof(int) * num_boxes);

    //     free(temp_box_start);
    //     free(temp_box_stop);
    }

    // void Grid::ComputeBeckeGrad(const double* E, double* grad) const
    // {
    //     BeckeGrad(num_points, points, num_atoms, atmXYZ, atmRad, E, grad);
    // }

    // What this function does is:
    // 1. The bounding box for a periodic system is the bounding box of the unit cell. We use square boxes regardless of unit cell shape.
    // 2. Split the bounding box into BOX_LENGTH*BOX_LENGTH*BOX_LENGTH small boxes
    // 3. Place all points inside a box next to each other, and order the points in terms of boxes
    //
    // See more documentation about the purpose of this function in gridbox/src/grid.cpp
    void PeriodicGrid::place_point_into_box(const int n_point, GridPoint* points, int* n_box, int* box_start, int* box_stop, const LatticeVector unit_cell)
    {
        const double upper_bound_x = fmax(unit_cell.real_space[0], fmax(unit_cell.real_space[3], unit_cell.real_space[6]));
        const double upper_bound_y = fmax(unit_cell.real_space[1], fmax(unit_cell.real_space[4], unit_cell.real_space[7]));
        const double upper_bound_z = fmax(unit_cell.real_space[2], fmax(unit_cell.real_space[5], unit_cell.real_space[8]));
        const double lower_bound_x = fmin(unit_cell.real_space[0], fmin(unit_cell.real_space[3], unit_cell.real_space[6]));
        const double lower_bound_y = fmin(unit_cell.real_space[1], fmin(unit_cell.real_space[4], unit_cell.real_space[7]));
        const double lower_bound_z = fmin(unit_cell.real_space[2], fmin(unit_cell.real_space[5], unit_cell.real_space[8]));
        const double edge_protection = 0.1;
        const double bound_size_x = upper_bound_x - lower_bound_x + edge_protection;
        const double bound_size_y = upper_bound_y - lower_bound_y + edge_protection;
        const double bound_size_z = upper_bound_z - lower_bound_z + edge_protection;
        const int Nx = (int)ceil(bound_size_x / BOX_LENGTH);
        const int Ny = (int)ceil(bound_size_y / BOX_LENGTH);
        const int Nz = (int)ceil(bound_size_z / BOX_LENGTH);

        int* point_to_box_map = (int*)malloc(sizeof(int) * n_point);
        int* point_in_box_count = (int*)malloc(sizeof(int) * Nx*Ny*Nz);
        memset(point_in_box_count, 0, sizeof(int) * Nx*Ny*Nz);
        for (int i = 0; i < n_point; i++) {
            const GridPoint& point = points[i];
            const double3 point_xyz { point.x, point.y, point.z };
            const double3 point_image = get_image_in_origin_cell(unit_cell, point_xyz);

            const int nx = (int)floor((point_image.x - lower_bound_x + edge_protection / 2) / (bound_size_x / Nx));
            const int ny = (int)floor((point_image.y - lower_bound_y + edge_protection / 2) / (bound_size_y / Ny));
            const int nz = (int)floor((point_image.z - lower_bound_z + edge_protection / 2) / (bound_size_z / Nz));

            if (nx < 0 || ny < 0 || nz < 0 || nx >= Nx || ny >= Ny || nz >= Nz) {
                printf("nxyz = %d, %d, %d, Nxyz = %d, %d, %d\n", nx, ny, nz, Nx, Ny, Nz);
                DIE("Logic error when constructing bounding box in PeriodicPairForGrid::place_pair_into_box()!\n");
            }

            const int i_box = nx*Ny*Nz + ny*Nz + nz;
            point_to_box_map[i] = i_box;
            point_in_box_count[i_box] += 1;
        }

        GridPoint* unsorted_points = (GridPoint*)malloc(sizeof(GridPoint) * n_point);
        memcpy(unsorted_points, points, sizeof(GridPoint) * n_point);

        int i_box = 0;
        int i_point  = 0;
        for (int b = 0; b < Nx*Ny*Nz; b++) {
            if (point_in_box_count[b] == 0)
                continue;

            box_start[i_box] = i_point;
            for (int p = 0; p < n_point; p++) {
                if (point_to_box_map[p] == b) {
                    points[i_point] = unsorted_points[p];
                    i_point += 1;
                }
            }
            box_stop[i_box] = i_point;
            i_box += 1;
        }
        *n_box = i_box;

        free(unsorted_points);
        free(point_in_box_count);
        free(point_to_box_map);
    }

}
