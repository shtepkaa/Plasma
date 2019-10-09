#ifndef PATCH_HPP
#define PATCH_HPP

/// FIXME /// #include "types.h"
/// FIXME /// #include "utility.h"
/// FIXME /// #include "hilbert.h"

#include <omp.h>

namespace maxwell {

//============================================================================//
//  Initialize markings
//============================================================================//
template<Dim dim>
void Patch::Initialize_markings()
{
    // Actual size products
    uint prods[dim] = { 1 };

    for (uint d = 1; d < dim; ++d)
    {
        prods[d] = prods[d - 1] * actual_sizes[d];
    }

    // Offsets along axes
    uint send_offsets[dim][3];
    uint recv_offsets[dim][3];
    
    for (uint d = 0; d < dim; ++d)
    {
        send_offsets[d][0] = ghost_width * prods[d];
        send_offsets[d][1] = send_offsets[d][0];
        send_offsets[d][2] = nominal_sizes[d] * prods[d];

        recv_offsets[d][1] = send_offsets[d][0];
        recv_offsets[d][2] = send_offsets[d][0] + send_offsets[d][2];
    }

    // For all ghosts
#pragma omp parallel for
    for (uint s = 0; s < ghost_markings.Get_size(); ++s)
    {
        uint index = s;

        // Calculate indices which are in { 0, 1, 2 }
        Tuple<dim> indices;

        for (uint d = 0; d < dim - 1; ++d)
        {
            indices[d] = index % 3;
            index /= 3;
        }

        indices[dim - 1] = index;

        // Identify Hilbert index
        ghost_markings[s].index
            = General_Hilbert_index<dim>(
                grid_exps, (grid_coords + indices - 1)
            );

        // Identify sizes
        for (uint d = 0; d < dim; ++d)
        {
            ghost_markings[s].sizes[d]
                = (indices[d] - 1)? ghost_width: nominal_sizes[indices[d]];
        }

        // Identify send and receive offsets
        ghost_markings[s].send_offset = 0;
        ghost_markings[s].recv_offset = 0;

        for (uint d = 0; d < dim; ++d)
        {
            ghost_markings[s].send_offset += send_offsets[d][indices[d]];
            ghost_markings[s].recv_offset += recv_offsets[d][indices[d]];
        }
    }
}

//============================================================================//
//  Construction
//============================================================================//
template<Dim dim, typename Type>
Patch::Patch(
    const Tuple<dim> & g_sizes,
    const Tuple<dim> & coords,
    const Tuple<dim> & sizes,
    const uint width
):
    grid_sizes(exps),
    grid_coords(coords),
    ghost_width(width),
    nominal_sizes(sizes),
    actual_sizes(sizes + 2 * width),
    markings(Total_hypercube_count(dim)),
    data_size(Product<dim>(actual_sizes)),
    data(NULL)
{
    Initialize_markings();
    cudaMalloc(&data, data_size * sizeof(Type));
}

} // namespace maxwell

#endif // PATCH_HPP
