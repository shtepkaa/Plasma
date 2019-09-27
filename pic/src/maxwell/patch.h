#ifndef PATCH_H
#define PATCH_H

/// FIXME /// #include "types.h"

namespace maxwell {

//============================================================================//
//  Patch
//============================================================================//
template<Dim dim, ArithmeticType Type>
class Patch
{
    private:

        // MPI rank
        uint rank;

        // Hilbert index
        uint index;

        // adjacent patches in lexicographic order
        Array<uint> neighbour_indices;

        // Cartesian grid coordinates
        uint grid_coords[dim];

        // Cartesian grid sizes
        uint grid_sizes[dim];

        // ghost width
        uint ghost_width;

        // device data size
        uint data_size;

        // device data
        Type * data;

        // identify own and neighbour hilbert indices
        void Identify_indices();

        /// FIXME /// Either remove or switch to c++11: = delete
        // default initializer
        Patch() {}

    public:

        // initialize
        Patch(const uint, const uint *, const uint *);

        // deallocate
        ~Patch();

        // set MPI rank
        inline void Set_rank(const uint r) { rank = r; }

        // obtain ghost
        void Get_ghost(const uint);
};

} // namespace maxwell

#include "patch.hpp"

#endif // PATCH_H
