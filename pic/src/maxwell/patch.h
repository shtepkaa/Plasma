#ifndef PATCH_H
#define PATCH_H

/// FIXME /// #include "types.h"
/// FIXME /// #include "utility.h"

namespace maxwell {

/*******************************************************************************
*
*   PatchMarking
*
*******************************************************************************/
// Implements patch marking
template<Dim dim>
struct PatchMarking
{
    uint index;

    Tuple<dim> sizes;

    uint send_offset;
    uint recv_offset;
};

/*******************************************************************************
*
*   Patch
*
*******************************************************************************/
// Implements patch -- hyperrectangular part of a grid in simulation box
// dim -- dimensionaliy of the patch
// Type -- supported arithmetic type
template<Dim dim, typename Type = double>
class Patch
{
    private:

        // Cartesian patch size exponents
        Tuple<dim> grid_exps;

        // Cartesian patch coordinates
        Tuple<dim> grid_coords;

        // nominal sizes: excluding ghosts
        Tuple<dim> grid_sizes;

        // ghost width
        uint ghost_width;

        // actual sizes: including ghosts
        Tuple<dim> actual_sizes;

        // ghost markings in lexicographic order
        Array<GhostMarking> ghost_markings;

        /// consideration ///
        // actual data size: including ghosts
        uint data_size;

        // __device__
        // patch data including ghosts
        Type * data;

        // initializes ghost markings
        void Initialize_markings();

        /// FIXME /// Either remove or switch to c++11: = delete
        // default initialization
        Patch() {}

    public:

        // initialization
        Patch(
            const Tuple<dim> &,
            const Tuple<dim> &,
            const Tuple<dim> &,
            const uint
        );

        // deallocation
        ~Patch() { cudaFree(data); }

        /// useful ///
        /// // get own id
        /// const PatchMarking & Get_id() const { return ids[ids.Get_size() >> 1]; }

        const Array<PatchMarking> & Get_markings() const { return markings; }
};

} // namespace maxwell

#include "patch.hpp"

#endif // PATCH_H
