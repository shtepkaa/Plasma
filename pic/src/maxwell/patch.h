#ifndef PATCH_H
#define PATCH_H

/// FIXME /// #include "types.h"
/// FIXME /// #include "utility.h"

namespace maxwell {

////////////////////////////////////////////////////////////////////////////////
//
//  GhostMarking
//
////////////////////////////////////////////////////////////////////////////////
// Implements patch marking
template<Dim dim>
struct GhostMarking
{
    uint index;

    Tuple<dim> sizes;

    uint send_offset;
    uint recv_offset;
};

////////////////////////////////////////////////////////////////////////////////
//
//  Patch
//
////////////////////////////////////////////////////////////////////////////////
// Implements patch (aka hyperrectangular part) of a grid in the simulation box
// Template parameter: dim -- dimensionality of the grid
// Template parameter: ord -- order of patch numeration
// Template parameter: Type -- supported arithmetic type
template<Dim dim, Order ord, typename Type = double>
class Patch
{
    private:

        //====================================================================//
        //  Geometry descriptors
        //====================================================================//
        // Layer sizes of the grid in patches
        Tuple<dim> layer_sizes;

        // Layer coordinates of the patch in a layer 
        Tuple<dim> layer_coordinates;

        // Sizes
        Tuple<dim> sizes;

        // Ghost width
        uint ghost_width;

        // Extended sizes including surrounding ghosts
        Tuple<dim> extended_sizes;

        // Ghost markings in lexicographic order
        Array<GhostMarking> ghost_markings;

        //====================================================================//
        //  Data
        //====================================================================//
        // Data size including surrounding ghosts,
        // which is the product of extended sizes
        uint data_size;

        // __device__
        // Patch data including surrounding ghosts
        Type * data;

        //====================================================================//
        //  Initialization methods
        //====================================================================//
        // Initializes ghost markings
        void Initialize_markings();

        /// FIXME /// Either remove or switch to c++11: = delete
        // default initialization
        Patch() {}

    public:

        // construction
        Patch(
            const Tuple<dim> &,
            const Tuple<dim> &,
            const Tuple<dim> &,
            const uint
        );

        // deallocation
        ~Patch() { CUDA_CALL(cudaFree(data)); }

        /// useful ///
        /// // get own id
        /// const PatchMarking & Get_id() const { return ids[ids.Get_size() >> 1]; }

        const Array<PatchMarking> & Get_markings() const { return markings; }
};

        inline Type & operator[](const uint ind) { return data[ind]; }
        inline const Type & operator[](const uint in) const { return data[in]; }

} // namespace maxwell

#include "patch.hpp"

#endif // PATCH_H
