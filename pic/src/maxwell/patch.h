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
// Contains ghost marking data
////////////////////////////////////////////////////////////////////////////////
template<Dimension dim, typename Type>
struct GhostMarking
{
    // Ghost sizes
    Tuple<dim> sizes;

    // Receiving ghost offset
    uint send_offset;

    // Sending ghost offset
    uint recv_offset;

    // Index of a patch to communicate with
    uint target_patch_index;

    // Index of a ghost in the patch to communicate with
    uint8_t target_ghost_index;

    GhostMarking();
};

/// ??? /// ////////////////////////////////////////////////////////////////////////////////
/// ??? /// //
/// ??? /// //  Ghost index and directions computation functions
/// ??? /// //
/// ??? /// ////////////////////////////////////////////////////////////////////////////////
/// ??? /// //==============================================================================
/// ??? /// //  Compute ghost index
/// ??? /// //==============================================================================
/// ??? /// template<Dimension dim>
/// ??? /// uint8_t Compute_ghost_index(const Tuple<dim, Direction> &);
/// ??? /// 
/// ??? /// //==============================================================================
/// ??? /// //  Compute ghost directions
/// ??? /// //==============================================================================
/// ??? /// template<Dimension dim>
/// ??? /// Tuple<dim, Direction> & Compute_ghost_directions(const uint8_t);
/// ??? /// 
/// ??? /// //==============================================================================
/// ??? /// //  Reflect ghost directions
/// ??? /// //==============================================================================
/// ??? /// template<Dimension dim>
/// ??? /// Tuple<dim, Direction> & Reflect_ghost_directions(Tuple<dim, Direction> &);

////////////////////////////////////////////////////////////////////////////////
//
//  PatchData
//
////////////////////////////////////////////////////////////////////////////////
// Contains a patch data
// -----------------------------------------------------------------------------  
// Template parameters:
//     dim -- dimensionality of the grid
//     ord -- order of patch numeration
//     Type -- supported data type
////////////////////////////////////////////////////////////////////////////////
/// FIXME /// Add a template parameter corresponding to internal patch data
/// FIXME /// storage order
template<Dimension dim, typename Type>
struct PatchData
{
    //==========================================================================
    //  Geometry descriptors
    //==========================================================================
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

    // Ghost marking index of the patch itself in the ghost markings array
    uint8_t ghost_index;

    //==========================================================================
    //  Data
    //==========================================================================
    // Data size (including surrounding ghosts sizes),
    // which is the product of extended sizes
    uint data_size;

    // __device__
    // Patch data (including surrounding ghosts)
    Type * data;

    //==========================================================================
    //  Data management
    //==========================================================================
    /// FIXME /// Either remove or switch to c++11: = delete
    PatchData() {}

    PatchData(
        const Tuple<dim> &, const Tuple<dim> &, const Tuple<dim> &, const uint,
    );

    ~PatchData() { if (data) { CUDA_CALL(cudaFree(data)); } }
};

//==============================================================================
//  Compute patch index
//==============================================================================
// Computes the index of the patch corresponding to a chosen order
template<Dimension dim, Order ord>
uint Compute_patch_index(const Tuple<dim> &, const Tuple<dim> &);

//==============================================================================
//  Compute patch coordinates
//==============================================================================
// Computes the coordinates of the patch corresponding to a chosen order
template<Dimension dim, Order ord>
Tuple<dim> Compute_patch_coordinates(const Tuple<dim> &, const uint);

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
template<Dimension dim, Order ord, typename Type>
class Patch
{
    private:

        PatchData<dim, Type> * data;

        void Set_ghost_markings(const uint);

    public:

        //======================================================================
        //  Data management
        //======================================================================
        Set_data(
            const Tuple<dim> &,
            const Tuple<dim> &,
            const Tuple<dim> &,
            const uint,
            const uint
        );

        Patch(): data(NULL) {}
        ~Patch() { if (data) { delete data; } }

        //======================================================================
        //  Access / mutate methods
        //======================================================================
        uint Get_index() const;

        const Array<GhostMarking> & Get_ghost_markings() const;

            /// ??? /// // Sets / gets an element
            /// ??? /// inline Type & operator[](const uint ind) { return data->data[ind]; }
            /// ??? /// inline const Type & operator[](const uint i) const { return data->data[i]; }

        // Sets / copies ghost data from / to a raw pointer to address
        void Set_ghost(const uint, const Type *);
        void Copy_ghost(const uint, Type *) const;
}

} // namespace maxwell

#include "patch.hpp"

#endif // PATCH_H
