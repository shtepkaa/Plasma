#ifndef PATCH_HPP
#define PATCH_HPP

/// FIXME /// #include "types.h"
/// FIXME /// #include "utility.h"

/// FIXME /// #include <omp.h>

namespace maxwell {

//============================================================================//
//  Construct markings
//============================================================================//
template<Dim dim>
void Patch::Initialize_markings()
{
    // actual size products
    uint prods[dim] = { 1 };

    for (uint d = 1; d < dim; ++d)
    {
        prods[d] = prods[d - 1] * actual_sizes[d];
    }

    // offsets along axes
    uint send_offsets[dim][3];
    uint recv_offsets[dim][3];
    
    for (uint d = 0; d < dim; ++d)
    {
        send_offsets[d][0] = ghost_width * prods[d];
        send_offsets[d][1] = send_offsets[d][0];
        send_offsets[d][2] = grid_sizes[d] * prods[d];

        recv_offsets[d][1] = send_offsets[d][0];
        recv_offsets[d][2] = send_offsets[d][0] + send_offsets[d][2];
    }

    // for all ghosts
#pragma omp parallel for
    for (uint s = 0; s < ghost_markings.Get_size(); ++s)
    {
        uint index = s;

        // calculate indices which are in { 0, 1, 2 }
        Tuple<dim> indices;

        for (uint d = 0; d < dim - 1; ++d)
        {
            indices[d] = index % 3;
            index /= 3;
        }

        indices[dim - 1] = index;

        // identify Hilbert index
        ghost_markings[s].index
            = General_Hilbert_index<dim>(
                grid_exps, (grid_coords + indices - 1)
            );

        // identify sizes
        for (uint d = 0; d < dim; ++d)
        {
            ghost_markings[s].sizes[d]
                = (indices[d] - 1)? ghost_width: grid_sizes[indices[d]];
        }

        // identify send and receive offsets
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
//  Constructor
//============================================================================//
template<Dim dim, typename Type>
Patch::Patch(
    const Tuple<dim> & exps,
    const Tuple<dim> & coords,
    const Tuple<dim> & sizes,
    const uint width
):
    grid_exps(exps),
    grid_coords(coords),
    grid_sizes(sizes),
    ghost_width(width),
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
