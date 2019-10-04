#ifndef PATCH_H
#define PATCH_H

/// FIXME /// #include "types.h"

namespace maxwell {

//============================================================================//
//  PatchMarking
//============================================================================//
// implements patch marking
template<Dim dim>
struct PatchMarking
{
    uint index;

    uint sizes[dim];

    uint send_offset;
    uint recv_offset;
};

//============================================================================//
//  Patch
//============================================================================//
template<Dim dim, ArithmeticType Type>
class Patch
{
    private:

        // Cartesian patch coordinates
        uint grid_exps[dim];

        // Cartesian patch coordinates
        uint grid_coords[dim];

        // nominal sizes: excluding ghosts
        uint grid_sizes[dim];

        // ghost width
        uint ghost_width;

        // actual sizes: including ghosts
        uint actual_sizes[dim];

            /// // actual data size: including ghosts
            /// uint data_size;

        // vicinity patch markings in lexicographic order
        Array<GhostMarking> markings;

        // __device__
        // patch data including ghosts
        Type * data;

        // identify own and neighbour Hilbert indices
        void Construct_markings();

        /// FIXME /// Either remove or switch to c++11: = delete
        // default initializer
        Patch() {}

    public:

        // initialize
        Patch(const * uint, const * uint, const * uint, const uint);

        // deallocate
        ~Patch() { cudaFree(data); }

            /// // get own id
            /// const PatchId & Get_id() const { return ids[ids.Get_size() >> 1]; }

        const Array<GhostMarking> & Get_markings() const { return markings; }
};

} // namespace maxwell

#include "patch.hpp"

#endif // PATCH_H
