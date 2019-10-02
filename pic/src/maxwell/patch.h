#ifndef PATCH_H
#define PATCH_H

/// FIXME /// #include "types.h"

namespace maxwell {

struct PatchId
{
    // MPI rank
    uint rank;

    // Hilbert index
    uint index;

    // data markup
    uint size;
    uint start;
};

//============================================================================//
//  Patch
//============================================================================//
template<Dim dim, ArithmeticType Type>
class Patch
{
    private:

            /// // adjacent patch type counts
            /// uint neighbour_type_counts[dim];

        // vicinity patch identifiers in lexicographic order
        Array<PatchId> ids;

        // Cartesian grid coordinates
        uint coords[dim];

        // Cartesian grid sizes
        uint sizes[dim];

        // ghost width
        uint ghost_width;

        // device data size
        uint data_size;

        // device data
        Type * data;

        // identify own and neighbour Hilbert indices
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

        /// // get own id
        /// const PatchId & Get_id() const { return ids[ids.Get_size() >> 1]; }
};

} // namespace maxwell

#include "patch.hpp"

#endif // PATCH_H
