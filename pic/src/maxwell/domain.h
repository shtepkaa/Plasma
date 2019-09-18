#ifndef PATCH_H
#define PATCH_H

/// FIXME /// #include "types.h"

#include <vector>

namespace maxwell {

//============================================================================//
//  Patch
//============================================================================//
template<Dim>
class Domain
{
    private:

        // MPI rank
        uint rank;

        // location
        bool location;

        // Hilbert index
        uint index;

        // Cartesian grid coordinates
        uint * grid_coords;

        // Cartesian grid sizes
        uint * grid_sizes;

        // ghost size
        uint ghost_size;

        // adjacent patches in lexicographic order
        Patch * neighbours;

        // host data
        double * host_data;

        // device data
        void * dev_data;

        // identify hilbert index
        void Get_index();

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
