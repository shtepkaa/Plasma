#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include <mpi.h>

namespace maxwell {

////////////////////////////////////////////////////////////////////////////////
//
//  TransferDescriptor
//
////////////////////////////////////////////////////////////////////////////////
TransferDescriptor::TransferDescriptor(
    const uint transfer_rank, const uint patch_count
):
    rank(transfer_rank),
    buffer_markings(patch_count),
    send_buffer(NULL),
    recv_buffer(NULL)
{} 

////////////////////////////////////////////////////////////////////////////////
//
//  Domain
//
////////////////////////////////////////////////////////////////////////////////
//============================================================================//
//  Initialize domain bounds
//============================================================================//
// Sets initial domain bounds
template<Dim dim, Order ord, typename Type>
void Domain::Initialize_domain_bounds()
{
    domain_bounds.Reallocate(range + 1);
    domain_bounds[0] = 0;

    if (ord == Cartesian)
    {
        /// TODO ///
    }
    else if (ord == Hilbertian)
    {
        uint total_patch_count
            = Product<dim>(grid_sizes) / Product<dim>(patch_sizes);

        uint patch_per_domain = total_patch_count / range;
        uint rem = total_patch_count % range;

        for (uint r = 1; r < range + 1; ++r)
        {
            domain_bounds[r]
                = domain_bounds[r - 1] + patch_per_domain + (r < rem);
        }
    }
}

//============================================================================//
//  Identify domain bounds
//============================================================================//
// Identifies domain bounds
template<Dim dim, Order ord, typename Type>
void Domain::Identify_domain_bounds()
{
    /// TODO ///
}

//============================================================================//
//  Identify transfer descriptors
//============================================================================//
// Identifies transfer descriptors
template<Dim dim, Order ord, typename Type>
void Domain::Identify_transfer_descriptors()
{
    // Allocate
    neighbour_ranks.Reallocate(range);

    uint neighbour_count = 0;
    uint trial_rank;
    uint trial_position;

    uint index;

    // For each patch in the domain
    for (int p = 0; p < patches.Get_size(); ++p)
    {
        const Array<GhostMarkings> & markings = patches[p].Get_ghost_markings();

        // Examine each ghost marking
        for (int m = 0; m < markings.Get_size(); ++m)
        {
            // If marking does not correspond to patch itself
            if (m != markings.Get_index())
            {
                // Identify MPI rank associated with the marking
                trial_rank = Binary_search(domain_bounds, markings.index);

                // Identify the position of the rank in neighbour ranks array
                trial_position = Binary_search(neighbour_ranks, trial_rank);

                // Try to add the trial rank to array of neighbour ranks
                if (
                    !neighbour_count
                    || trial_rank != neighbour_ranks[trial_position] 
                )
                {
                    neighbour_ranks.Insert(
                        trial_rank, trial_position, neighbour_count
                    );

                    ++neighbour_count;
                }
            }
        }
    }

    // Shrink to fit
    neighbour_ranks.Reallocate(neighbour_count);
}

//============================================================================//
//  Constructor
//============================================================================//
template<Dim dim, Order ord, typename Type>
Domain::Domain(const Tuple<dim> & grid_sizs, const Tuple<dim> & patch_sizs):
    rank(0),
    range(0),
    domain_bounds(),
    grid_sizes(grid_sizs),
    patch_sizes(patch_sizs),
    patches(),
    local_markings(),
    transfer_descriptor()
{
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &range);

    Initialize_domain_bounds();
    patches.Reallocate(Get_patch_max_index() - Get_patch_min_index());

    if (range)
    {
        Identify_neighbour_ranks();
        Identify_comm_descriptors();
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
