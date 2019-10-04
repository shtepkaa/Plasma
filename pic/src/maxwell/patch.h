#ifndef PATCH_H
#define PATCH_H

/// FIXME /// #include "types.h"

namespace maxwell {

//============================================================================//
//  GhostMarking
//============================================================================//
// marking structure
template<Dim dim>
struct GhostMarking
{
    uint send_ind;
    uint recv_ind;

    uint sizes[dim];
    uint lead_sizes[dim];

    Type * send_ghost;
    Type * recv_ghost;
};

//============================================================================//
//  Patch
//============================================================================//
template<Dim dim, ArithmeticType Type>
class Patch
{
    private:

        // vicinity patch markings in lexicographic order
        Array<GhostMarking> markings;

        // Cartesian patch coordinates
        uint coords[dim];

        // nominal sizes: excluding ghosts
        uint sizes[dim];

        // ghost width
        uint ghost_width;

        // actual data size: including ghosts
        uint data_size;

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
        Patch(const uint, const uint *, const uint *);

        // deallocate
        ~Patch();

            /// // get own id
            /// const PatchId & Get_id() const { return ids[ids.Get_size() >> 1]; }

        const Array<uint> & Get_markings() const { return markings; }
};

} // namespace maxwell

#include "patch.hpp"

#endif // PATCH_H
