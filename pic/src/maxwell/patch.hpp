#ifndef PATCH_HPP
#define PATCH_HPP

/// FIXME /// #include "types.h"
/// FIXME /// #include "utility.h"
/// FIXME /// #include "hilbert.h"

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
        Tuple<dim> exps(sizes, Logarithm);

        return General_Hilbert_index<dim>(exps, coords);
    }
}

//============================================================================//
//  Initialize markings
//============================================================================//
template<Dim dim, Order ord, typename Type>
void Patch::Initialize_markings()
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
    for (int s = 0; s < ghost_markings.Get_size(); ++s)
    {
        uint index = s;

        // Calculate indices which are in { 0, 1, 2 }
        Tuple<dim> indices;

        for (int d = 0; d < dim - 1; ++d)
        {
            indices[d] = index % 3;
            index /= 3;
        }

        indices[dim - 1] = index;

        // Identify index
        ghost_markings[s].index
            = Index<dim, ord>(layer_sizes, (layer_coordinates + indices - 1));

        // Identify sizes
        for (int d = 0; d < dim; ++d)
        {
            ghost_markings[s].sizes[d]
                = (indices[d] - 1)? ghost_width: sizes[indices[d]];
        }

        // Identify send and receive offsets
        ghost_markings[s].send_offset = 0;
        ghost_markings[s].recv_offset = 0;

        for (int d = 0; d < dim; ++d)
        {
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
    const uint width
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
    Initialize_markings();
    CUDA_CALL(cudaMalloc(&data, data_size * sizeof(Type)));
}

} // namespace maxwell

#endif // PATCH_HPP
