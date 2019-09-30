#ifndef PATCH_H
#define PATCH_H

/// FIXME /// #include "types.h"

namespace maxwell {

struct PatchId
{
    uint rank;
    uint index;
    uint size;
};

//============================================================================//
//  Patch
//============================================================================//
template<Dim dim, ArithmeticType Type>
class Patch
{
    private:

            /// // MPI rank
            /// uint rank;

            /// // Hilbert index
            /// uint index;

            /// // adjacent patch type counts
            /// uint neighbour_type_counts[dim];

        // vicinity Hilbert indices in lexicographic order
        Array<PatchId> indices;

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

        // get own Hilbert index
        uint Get_index() const { return indices[indices.Get_size() >> 1]; }
};

} // namespace maxwell

#include "patch.hpp"

#endif // PATCH_H
