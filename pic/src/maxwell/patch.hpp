#ifndef PATCH_HPP
#define PATCH_HPP

#include "types.h"
#include "utility.h"
#include "hilbert.h"

#include <omp.h>

namespace maxwell {

//============================================================================//
//  Index
//============================================================================//
// Computes index of the patch corresponding to a chosen order
template<Dim dim, Order ord>
static uint Index(const Tuple<dim> & sizes, const Tuple<dim> & coords)
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
//  Multiindex
//============================================================================//
// 
template<Dim dim, uint base>
static void Multiindex(const uint ind, Tuple<dim> & multiind)
{
    uint tmp = ind;

    for (int d = 0; d < dim - 1; ++d)
    {
        multiind[d] = tmp % base;
        tmp /= base;
    }

    multiind[dim - 1] = tmp;
}

//============================================================================//
//  Initialize_markings
//============================================================================//
template<Dim dim, Order ord, typename Type>
void Patch::Initialize_markings(const uint index)
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

    Tuple<dim> indices;

    // For all ghosts
#pragma omp parallel for
    for (int s = 0; s < ghost_markings.Get_size(); ++s)
    {
        // Calculate indices which are in { 0, 1, 2 }
        Multiindex<dim, 3>(s, indices);

        // Identify index
        ghost_markings[s].index
            = (s == Get_index())?
                index:
                Index<dim, ord>(layer_sizes, (layer_coordinates + indices - 1));

        ghost_markings[s].send_offset = 0;
        ghost_markings[s].recv_offset = 0;

        // Identify sizes, send and receive offsets
        for (int d = 0; d < dim; ++d)
        {
            ghost_markings[s].sizes[d]
                = (indices[d] - 1)? ghost_width: sizes[indices[d]];

            ghost_markings[s].send_offset += send_offsets[d][indices[d]];
            ghost_markings[s].recv_offset += recv_offsets[d][indices[d]];
        }
    }
}

//============================================================================//
//  Construction
//============================================================================//
template<Dim dim, Order ord, typename Type>
Patch::Patch(
    const Tuple<dim> & layer_sizs,
    const Tuple<dim> & layer_coords,
    const Tuple<dim> & sizs,
    const uint width,
    const uint index
):
    layer_sizes(layer_sizs),
    layer_coordinates(layer_coords),
    sizes(sizs),
    ghost_width(width),
    extended_sizes(sizes + 2 * width),
    ghost_markings(Power(3, dim)),
    data_size(Product<dim>(extended_sizes)),
    data(NULL)
{
    Initialize_markings(index);
    CUDA_CALL(cudaMalloc(&data, data_size * sizeof(Type)));
}

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
void Patch::Get_ghost(const uint ind, Type * buf) const
{
    /// TODO /// Cuda copy from device to device
}

} // namespace maxwell

#endif // PATCH_HPP
