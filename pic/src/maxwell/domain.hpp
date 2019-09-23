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
    neighbour_count(0),
    neighbour_ranks(NULL),
    patch_min_index(min_index),
    patch_max_index(max_index),
    patch_count(max_index - min_index),
    patches(NULL),
    patch_sizes(NULL),
    border_patch_count(0),
    border_patch_indices(NULL),
    recv_buffer_markups(NULL),
    recv_buffers(NULL),
    send_buffer_markups(NULL),
    send_buffers(NULL)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &range);

    Init_array(dim, sizes, patch_sizes);
    Init_array(patch_count, patches, NULL);

    for (int p = 0; p < patch_count; ++p)
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
