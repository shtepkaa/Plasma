#ifndef DOMAIN_H
#define DOMAIN_H

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
        Domain();
        Domain(const uint *, const uint *);

        // deallocate
        ~Domain();

        // set MPI rank
        void Set_rank(const uint);

        // obtain ghost
        void Get_ghost(const uint);
};

} // namespace maxwell

#include "domain.hpp"

#endif // DOMAIN_H
