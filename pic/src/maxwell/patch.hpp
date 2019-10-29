#ifndef PATCH_HPP
#define PATCH_HPP

#include "types.h"
#include "utility.h"
#include "hilbert_index.h"

#include <omp.h>

namespace maxwell {

////////////////////////////////////////////////////////////////////////////////
//
//  GhostMarking struct
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
//  Ghost index and directions computation functions
//
////////////////////////////////////////////////////////////////////////////////
//==============================================================================
//  Compute ghost index
//==============================================================================
template<Dimension dim>
static uint8_t Compute_ghost_index(const Tuple<dim, Direction> & directions)
{
    uint8_t res = directions[0];

    for (int d = 1; d < dim; ++d)
    {
        res *= 3;
        res += directions[d];
    }

    return res;
}

//==============================================================================
//  Compute ghost directions
//==============================================================================
template<Dimension dim>
static Tuple<dim, Direction> & Compute_ghost_directions(const uint8_t ind)
{
    Tuple<dim, Direction> res;

    uint8_t tmp = ind;

    for (int d = 0; d < dim - 1; ++d)
    {
        res[d] = tmp % 3;
        tmp /= 3;
    }

    res[dim - 1] = tmp;

    return res;
}

//==============================================================================
//  Reflect ghost directions
//==============================================================================
template<Dimension dim>
static Tuple<dim, Direction> & Reflect_ghost_directions(
    Tuple<dim, Direction> & directions
)
{
    for (int d = 0; d < dim; ++d)
    {
        if (directions[d] == LEFT) { directions[d] = RIGHT; }
        else if (directions[d] == RIGHT) { directions[d] = LEFT; }
    }

    return directions;
}

////////////////////////////////////////////////////////////////////////////////
//
//  PatchData
//
////////////////////////////////////////////////////////////////////////////////
template<Dimension dim, typename Type>
PatchData::PatchData(
    const Tuple<dim> & relative_layer_sizes,
    const Tuple<dim> & patch_coordinates,
    const Tuple<dim> & patch_sizes,
    const uint width
):
    layer_sizes(relative_layer_sizes),
    coordinates(patch_coordinates),
    sizes(patch_sizes),
    ghost_width(width),
    extended_sizes(sizes + 2 * width),
    ghost_markings(Power(3, dim)),
    ghost_index(ghost_markings().Get_size() >> 1),
    data_size(Product<dim>(extended_sizes)),
    data(NULL)
{
    CUDA_CALL(cudaMalloc(&data, data_size * sizeof(Type)));
}

////////////////////////////////////////////////////////////////////////////////
//
//  Patch
//
////////////////////////////////////////////////////////////////////////////////
//==============================================================================
//  Set data
//==============================================================================
template<Dimension dim, Order ord, typename Type>
void Set_data(
    const Tuple<dim> & relative_layer_sizes,
    const Tuple<dim> & patch_coordinates,
    const Tuple<dim> & patch_sizes,
    const uint width,
    const uint index
)
{
    /// FIXME /// Probably better to change existing data, than delete/allocate
    if (data) { CUDA_CALL(cudaFree(data)); }

    data
        = new PatchData(
            relative_layer_sizes, patch_coordinates, patch_sizes, width
        );

    Set_ghost_markings(index);
}

//==============================================================================
//  Compute patch index
//==============================================================================
// Computes the index of the patch corresponding to a chosen order
template<Dimension dim, Order ord>
static uint Compute_patch_index(
    const Tuple<dim> & sizes, const Tuple<dim> & coords
)
{
    if (ord == CARTESIAN)
    {
        /// FIXME ///
        return Cartesian_index<dim>(sizes, coords);
    }
    else if (ord == HILBERTIAN)
    {
        // Initialize by exponents of the given sizes
        Tuple<dim> exps(sizes, Binary_logarithm);

        return General_Hilbert_index<dim>(exps, coords);
    }
}

//==============================================================================
//  Compute patch coordinates
//==============================================================================
// Computes the coordinates of the patch corresponding to a chosen order
template<Dimension dim, Order ord>
static uint Compute_patch_coordinates(
    const Tuple<dim> & sizes, const uint index, Tuple<dim> & coords
)
{
    if (ord == CARTESIAN)
    {
        /// FIXME ///
        Inverse_Cartesian_index<dim>(sizes, coords, index);
    }
    else if (ord == HILBERTIAN)
    {
        // Initialize by exponents of the given sizes
        Tuple<dim> exps(sizes, Binary_logarithm);

        // Convert the index to coordinates
        Inverse_general_Hilbert_index<dim>(exps, coords, index);
    }
}

//==============================================================================
//  Set ghost markings
//==============================================================================
template<Dimension dim, Order ord, typename Type>
void Patch::Set_ghost_markings(const uint index)
{
    // Compute extended size products
    uint prods[dim] = { 1 };

    for (int d = 1; d < dim; ++d)
    {
        prods[d] = prods[d - 1] * data->extended_sizes[d];
    }

    // Compute offsets along all axes
    uint send_offsets[dim][3];
    uint recv_offsets[dim][3];
    
    for (int d = 0; d < dim; ++d)
    {
        send_offsets[d][0] = data->ghost_width * prods[d];
        send_offsets[d][1] = send_offsets[d][0];
        send_offsets[d][2] = data->sizes[d] * prods[d];

        recv_offsets[d][1] = send_offsets[d][0];
        recv_offsets[d][2] = send_offsets[d][0] + send_offsets[d][2];
    }

    // For all ghosts
#pragma omp parallel for
    for (uint8_t g = 0; g < data->ghost_markings.Get_size(); ++g)
    {
        GhostMarking & marking = data->ghost_markings[g];

        // Convert the ghost index to a multiindex of directions
        Tuple<dim, Direction> & directions = Compute_ghost_directions<dim>(g);

        // Set target patch index
        marking.target_patch_index
            = (g == data->ghost_index)?
                index:
                Compute_patch_index<dim, ord>(
                    data->layer_sizes, (data->coordinates + directions - 1)
                );

        // Set sizes, sending and receiving offsets
        for (int d = 0; d < dim; ++d)
        {
            marking.send_offset += send_offsets[d][directions[d]];
            marking.recv_offset += recv_offsets[d][directions[d]];

            marking.sizes[d]
                = (directions[d] - 1)?
                    data->ghost_width:
                    data->sizes[directions[d]];
        }

        // Set target ghost index
        marking.target_ghost_index
            = Compute_ghost_index<dim>(
                Reflect_ghost_directions<dim>(directions)
            );
    }
}

//==============================================================================
//  Get index
//==============================================================================
template<Dimension dim, Order ord, typename Type>
uint Patch::Get_index() const
{
    return data->ghost_markings[data->ghost_index].target_patch_index;
}

//==============================================================================
//  Get ghost markings
//==============================================================================
template<Dimension dim, Order ord, typename Type>
const Array<GhostMarking> & Patch::Get_ghost_markings() const
{
    return data->ghost_markings;
}

//==============================================================================
//  Set ghost
//==============================================================================
template<Dimension dim, Order ord, typename Type>
void Patch::Set_ghost(const uint ind, const Type * buf)
{
    /// TODO /// Cuda copy from device to device
}

//==============================================================================
//  Get ghost
//==============================================================================
template<Dimension dim, Order ord, typename Type>
void Patch::Copy_ghost(const uint ind, Type * buf) const
{
    /// TODO /// Cuda copy from device to device
}

} // namespace maxwell

#endif // PATCH_HPP
