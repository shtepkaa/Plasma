#ifndef PATCH_HPP
#define PATCH_HPP

/// FIXME /// #include "types.h"

namespace maxwell {

template<Dim>
class Patch
{
    private:

        // MPI rank
        uint rank;

        // Hilbert index
        uint index;

        // Cartesian grid coordinates
        uint * coordinates;

        // Cartesian grid sizes
        uint * sizes;

        // ghost size
        uint ghost_size;

        // adjacent patches in lexicographic order
        Patch * neighbours;

        // host data
        double * host_data;

        // device data
        void * dev_data;

        // identify neighbours
        void Get_neighbours();

    public:

        // initialize
        Patch();

        // deallocate
        ~Patch();

        // obtain ghost
        void Get_ghost(const uint);
};

} // namespace maxwell

#endif // PATCH_HPP
