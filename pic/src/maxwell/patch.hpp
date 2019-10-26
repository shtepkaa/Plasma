#ifndef PATCH_HPP
#define PATCH_HPP

#include "types.h"
#include "utility.h"
#include "hilbert.h"

#include <omp.h>

namespace maxwell {

////////////////////////////////////////////////////////////////////////////////
//
//  GhostMarking
//
////////////////////////////////////////////////////////////////////////////////
GhostMarking::GhostMarking():
    sizes(),
    send_offset(0),
    recv_offset(0),
    target_patch_index(~0),
    target_ghost_index(~0)
{}

////////////////////////////////////////////////////////////////////////////////
//
//  Patch
//
////////////////////////////////////////////////////////////////////////////////
//============================================================================//
//  Compute_patch_index
//============================================================================//
// Computes index of the patch corresponding to a chosen order
template<Dim dim, Order ord>
static uint Compute_patch_index(
    const Tuple<dim> & sizes, const Tuple<dim> & coords
)
{
    /// FIXME ///
    if (ord == Cartesian) { return Cartesian_index<dim>(sizes, coords); }
    else if (ord == Hilbertian)
    {
        // Initialize by exponents of the given sizes
        Tuple<dim> exps(sizes, Binary_logarithm);

        return General_Hilbert_index<dim>(exps, coords);
    }
}

//============================================================================//
//  Initialize_ghost_markings
//============================================================================//
template<Dim dim, Order ord, typename Type>
void Patch::Initialize_ghost_markings(const uint index)
{
    // Extended size products
    uint prods[dim] = { 1 };

    for (int d = 1; d < dim; ++d)
    {
        prods[d] = prods[d - 1] * extended_sizes[d];
    }

    // Offsets along axes
    uint send_offsets[dim][3];
    uint recv_offsets[dim][3];
    
    for (int d = 0; d < dim; ++d)
    {
        send_offsets[d][0] = ghost_width * prods[d];
        send_offsets[d][1] = send_offsets[d][0];
        send_offsets[d][2] = sizes[d] * prods[d];

        recv_offsets[d][1] = send_offsets[d][0];
        recv_offsets[d][2] = send_offsets[d][0] + send_offsets[d][2];
    }

    // For all ghosts
#pragma omp parallel for
    for (uint8_t g = 0; g < ghost_markings.Get_size(); ++g)
    {
        GhostMarking & marking = ghost_markings[g];

        // Compute ghost index from multiindex of directions
        Tuple<dim, Dir> & directions = Identify_ghost_directions<dim>(g);

        // Identify patch index
        marking.patch_index
            = (g == Get_index())?
                index:
                Compute_patch_index<dim, ord>(
                    layer_sizes, (coordinates + directions - 1)
                );

        // Identify sizes, sending and receiving offsets
        for (int d = 0; d < dim; ++d)
        {
            marking.send_offset += send_offsets[d][directions[d]];
            marking.recv_offset += recv_offsets[d][directions[d]];

            marking.sizes[d]
                = (directions[d] - 1)? ghost_width: sizes[directions[d]];
        }

        // Identify receiving ghost index
        marking.recv_ghost_index
            = Identify_ghost_index(Reflect_ghost_directions(directions));
    }
}

//============================================================================//
//  Construction
//============================================================================//
template<Dim dim, Order ord, typename Type>
Patch::Patch(
    const Tuple<dim> & relative_layer_sizes,
    const Tuple<dim> & coords,
    const Tuple<dim> & patch_sizes,
    const uint width,
    const uint index
):
    layer_sizes(relative_layer_sizes),
    coordinates(coords),
    sizes(patch_sizes),
    ghost_width(width),
    extended_sizes(sizes + 2 * width),
    ghost_markings(Power(3, dim)),
    marking_index(ghost_markings().Get_size() >> 1),
    data_size(Product<dim>(extended_sizes)),
    data(NULL)
{
    Initialize_markings(index);
    CUDA_CALL(cudaMalloc(&data, data_size * sizeof(Type)));
}

//============================================================================//
// Get_index
//============================================================================//
uint Get_index() const { return ghost_markings[marking_index].patch_index; }

//============================================================================//
//  Get_ghost_markings
//============================================================================//
template<Dim dim, Order ord, typename Type>
const Array<GhostMarking> & Patch::Get_ghost_markings() const
{
    return ghost_markings;
}

//============================================================================//
//  Set_ghost
//============================================================================//
void Patch::Set_ghost(const uint ind, const Type * buf)
{
    /// TODO /// Cuda copy from device to device
}

//============================================================================//
//  Get_ghost
//============================================================================//
void Patch::Copy_ghost(const uint ind, Type * buf) const
{
    /// TODO /// Cuda copy from device to device
}

} // namespace maxwell

#endif // PATCH_HPP
