#ifndef DOMAIN_HPP
#define DOMAIN_HPP

/// FIXME /// #include "types.h"

#include <mpi.h>

namespace maxwell {

//============================================================================//
//  Set neighbour ranks
//============================================================================//
// Identifies neighbour ranks
template<Dim dim>
void Domain::Set_neighbour_ranks()
{
}

//============================================================================//
//  Set border patch indices
//============================================================================//
// Identifies border patch indices
template<Dim dim>
void Domain::Set_border_patch_indices()
{
}

//============================================================================//
//  Constructor
//============================================================================//
template<Dim dim>
Domain::Domain(
    const Tuple<dim> & sizes
):
    rank(0),
    range(0),
    domain_bounds(),
    neighbour_ranks(),
    patch_sizes(sizes),
    patches(Get_patch_max_index() - Get_patch_min_index()),
    local_buffer_markings(),
    recv_buffer_markings(),
    recv_buffers(),
    send_buffer_markings(),
    send_buffers()
{
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &range);

    domain_bounds.Reallocate(range + 1);
    neighbour_ranks.Reallocate(range);

    for (int p = 0; p < patches.Get_size(); ++p)
    {
        const Array<GhostMarkings> & markings = patches[p].Get_ghost_markings();

        for (int m = 0; m < markings.Get_size(); ++m)
        {
            if (markings.index != Get_patch_min_index + p)
            {
            }
            else
            {
            }
        }
    }
}

//============================================================================//
//  Deallocation
//============================================================================//
template<Dim dim, Order ord, typename Type>
Domain::~Domain()
{
    for (uint p = 0; p < patches.Get_size(); ++p)
    {
        delete patches[p];

        CUDA_CALL(cudaFree(recv_buffers[p]));
        CUDA_CALL(cudaFree(send_buffers[p]));
    }
}

} // namespace maxwell

#endif // DOMAIN_HPP
