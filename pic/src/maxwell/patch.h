#ifndef PATCH_H
#define PATCH_H

/// FIXME /// #include "types.h"

namespace maxwell {

//============================================================================//
//  Patch
//============================================================================//
template<Dim>
class Patch
{
    private:

        // MPI rank
        uint rank;

        // host location identificator
        bool location;

        // Hilbert index
        uint index;

        // Cartesian grid coordinates
        uint * grid_coords;

        // Cartesian grid sizes
        uint * grid_sizes;

        // ghost width
        uint ghost_width;

        // adjacent patches in lexicographic order
        Array<uint> neighbour_indices;

        // device data
        Array<ArithmeticType, Device> data;

        // identify hilbert index
        uint Get_index();

        // identify neighbours
        void Get_neighbours();

    public:

        // initialize
        Patch();
        Patch(const uint *, const uint *);

        // deallocate
        ~Patch();

        // set MPI rank
        void Set_rank(const uint);

        // obtain ghost
        void Get_ghost(const uint);
};

} // namespace maxwell

#include "patch.hpp"

#endif // PATCH_H
