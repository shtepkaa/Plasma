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
// Implements patch (aka hyperrectangular part) of a grid in the simulation box
// template parameter: dim -- dimensionaliy of the patch
// template parameter: Type -- supported arithmetic type
template<Dim dim, typename Type = double>
class Patch
{
    private:

        // Cartesian grid sizes 
        Tuple<dim> grid_sizes;

        // Cartesian grid patch coordinates
        Tuple<dim> coords;

        // ghost width
        uint ghost_width;

        // nominal sizes: excluding ghosts
        Tuple<dim> nominal_sizes;

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
