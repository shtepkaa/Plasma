#ifndef DOMAIN_HPP
#define DOMAIN_HPP

/// FIXME /// #include "types.h"

#include <mpi.h>

namespace maxwell {

// identify neighbour ranks
template<Dim dim>
void Domain::Set_neighbour_ranks()
{
}

// identify border patch indices
template<Dim dim>
void Domain::Set_border_patch_indices()
{
}

// initialize
template<Dim dim>
Domain::Domain(const uint min_index, const uint max_index, const uint * sizes):
    rank(0),
    range(0),
    neighbour_ranks(),
    patch_min_index(min_index),
    patch_max_index(max_index),
    patches(max_index - min_index, NULL),
    patch_sizes(dim, sizes),
    border_patch_indices(),
    recv_buffer_markups(),
    recv_buffers(),
    send_buffer_markups(),
    send_buffers()
{
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &range);

    for (int p = 0; p < patches.Get_size(); ++p)
    {
    }
}

// deallocate
template<Dim dim>
Domain::~Domain()
{
}

} // namespace maxwell

#endif // DOMAIN_HPP
