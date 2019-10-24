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

    // Ghost sizes along all axes
    Tuple<dim> sizes;

    // Sending ghost directions along all axes in { Left, Center, Right }
    Tuple<dim, Dir> directions; 

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

        // Ghost markings in lexicographic order
        Array<GhostMarking> ghost_markings;

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
        // Gets index of the patch
        inline uint Get_index() const { return ghost_markings.Get_size() >> 1; }

        // Gets ghost markings
        const Array<GhostMarking> & Get_ghost_markings() const;

        // Set / get element
        inline Type & operator[](const uint ind) { return data[ind]; }
        inline const Type & operator[](const uint i) const { return data[i]; }

        // Set / get ghost array
        void Set_ghost(const uint, const Type *);
        void Get_ghost(const uint, Type *) const;
};

} // namespace maxwell

#include "patch.hpp"

#endif // PATCH_H
