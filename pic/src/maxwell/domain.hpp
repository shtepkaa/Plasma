#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include <mpi.h>

namespace maxwell {

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
//  Identify communication descriptors
//============================================================================//
// Identifies communication descriptors
template<Dim dim, Order ord, typename Type>
void Domain::Identify_comm_descriptors()
{
    neighbour_ranks.Reallocate(range); /// !!! /// Overkill
    patches.Reallocate(Get_patch_max_index() - Get_patch_min_index());

    uint neighbour_count = 0;

    uint index;
    uint neighbour_rank;

    for (int p = 0; p < patches.Get_size(); ++p)
    {
        const Array<GhostMarkings> & markings = patches[p].Get_ghost_markings();

        for (int m = 0; m < markings.Get_size(); ++m)
        {
            index = markings.index;

            if (index != Get_patch_min_index() + p)
            {
                for (int n = 0; n < neighbour_count; ++n)
                {
                    neighbour_rank
                        = Binary_search(
                            domain_bounds.Get_size(),
                            domain_bounds.Export_data(),
                            index
                        );
                }
            }
            else
            {
            }
        }
    }
}

//============================================================================//
//  Constructor
//============================================================================//
template<Dim dim, Order ord, typename Type>
Domain::Domain(const Tuple<dim> & grid_sizs, const Tuple<dim> & patch_sizs):
    rank(0),
    range(0),
    domain_bounds(),
    neighbour_ranks(),
    grid_sizes(grid_sizs),
    patch_sizes(patch_sizs),
    patches(),
    local_markings(),
    recv_buffer_markings(),
    recv_buffers(),
    send_buffer_markings(),
    send_buffers()
{
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &range);

    Initialize_domain_bounds();

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
