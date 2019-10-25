#ifndef PATCH_H
#define PATCH_H

#include "types.h"
#include "utility.h"

namespace maxwell {

////////////////////////////////////////////////////////////////////////////////
//
//  GhostMarking
//
////////////////////////////////////////////////////////////////////////////////
// Implements patch marking
////////////////////////////////////////////////////////////////////////////////
template<Dim dim, typename Type>
struct GhostMarking
{
    // Receiving patch index
    uint patch_index;

    // Receiving ghost offset
    uint send_offset;

    // Sending ghost offset
    uint recv_offset;

    // Ghost sizes
    Tuple<dim> sizes;

    // Receiving ghost index
    uint8_t recv_ghost_index;

    GhostMarking();
};

////////////////////////////////////////////////////////////////////////////////
//
//  Patch
//
////////////////////////////////////////////////////////////////////////////////
// Implements a patch (aka hyperrectangular part) of the grid corresponding
// to the simulation box
// -----------------------------------------------------------------------------  
// Template parameters:
//     dim -- dimensionality of the grid
//     ord -- order of patch numeration
//     Type -- supported arithmetic type
////////////////////////////////////////////////////////////////////////////////
template<Dim dim, Order ord, typename Type>
class Patch
{
    private:

        //======================================================================
        //  Geometry descriptors
        //======================================================================
        // The layer sizes of the grid (in patches)
        Tuple<dim> layer_sizes;

        // Patch coordinates in the layer 
        Tuple<dim> coordinates;

        // Patch sizes (excluding surrounding ghosts sizes)
        Tuple<dim> sizes;

        // Ghost width
        uint ghost_width;

        // Extended sizes (including surrounding ghosts sizes)
        Tuple<dim> extended_sizes;

        // Ghost markings in lexicographic order (including the patch itself)
        Array<GhostMarking> ghost_markings;

        // Patch marking index of the patch in the ghost markings array
        uint8_t marking_index;

        //======================================================================
        //  Data
        //======================================================================
        // Data size (including surrounding ghosts sizes),
        // which is the product of extended sizes
        uint data_size;

        // __device__
        // Patch data (including surrounding ghosts)
        Type * data;

        //======================================================================
        //  Initialization methods
        //======================================================================
        void Initialize_ghost_markings(const uint);

        // Default
        /// FIXME /// Either remove or switch to c++11: = delete
        Patch() {}

    public:

        // Construction
        Patch(
            const Tuple<dim> &,
            const Tuple<dim> &,
            const Tuple<dim> &,
            const uint,
            const uint
        );

        // Deallocation
        ~Patch() { CUDA_CALL(cudaFree(data)); }

        //======================================================================
        //  Access / mutate methods
        //======================================================================
        /// ??? /// // Gets patch marking index in ghost markings array
        /// ??? /// inline uint Get_marking_index() const { return marking_index; }

        // Gets index of the patch
        uint Get_index() const;

        // Gets ghost markings
        const Array<GhostMarking> & Get_ghost_markings() const;

        // Sets / gets element
        inline Type & operator[](const uint ind) { return data[ind]; }
        inline const Type & operator[](const uint i) const { return data[i]; }

        // Sets / copies ghost data
        void Set_ghost(const uint, const Type *);
        void Copy_ghost(const uint, Type *) const;
};

} // namespace maxwell

#include "patch.hpp"

#endif // PATCH_H
