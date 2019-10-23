#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include <mpi.h>

namespace maxwell {

////////////////////////////////////////////////////////////////////////////////
//
//  TransferDescriptorData
//
////////////////////////////////////////////////////////////////////////////////
TransferDescriptorData::TransferDescriptorData(const uint transfer_rank):
    rank(transfer_rank),
    buffer_markings(),
    size(0),
    send_buffer(NULL),
    recv_buffer(NULL)
{}

TransferDescriptorData::~TransferDescriptorData()
{
    CUDA_CALL(cudaFree(recv_buffer));
    CUDA_CALL(cudaFree(send_buffer));
}

////////////////////////////////////////////////////////////////////////////////
//
//  TransferDescriptor
//
////////////////////////////////////////////////////////////////////////////////
/// ??? /// TransferDescriptor::TransferDescriptor(const uint transfer_rank):
/// ??? ///     data(new TransferDescriptorData(transfer_rank))
/// ??? /// {}

//==============================================================================
//  Data management
//==============================================================================
void TransferDescriptor::Set_data(const uint transfer_rank)
{
    data = new TransferDescriptorData(transfer_rank);
}

//==============================================================================
//  Access / mutate methods
//==============================================================================
Array<BufferMarking> & TransferDescriptor::Get_buffer_markings()
{
    return data->buffer_markings;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Domain
//
////////////////////////////////////////////////////////////////////////////////
//==============================================================================
//  Initialize domain bounds
//==============================================================================
/// FIXME ///
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

        for (uint r = 1; r <= range; ++r)
        {
            domain_bounds[r]
                = domain_bounds[r - 1] + patch_per_domain + (r < rem);
        }
    }
}

//==============================================================================
//  Identify domain bounds
//==============================================================================
// Collects the new domain bounds from all processes
template<Dim dim, Order ord, typename Type>
void Domain::Identify_domain_bounds()
{
    MPI_Allgather(
        (void *)(domain_bounds + rank), 1, MPI_UNSIGNED, (void *)domain_bounds,
        1, MPI_UNSIGNED, MPI_COMM_WORLD
    );
}

//==============================================================================
//  Comparison
//==============================================================================
// This is required for the binary search routine
template<typename Type>
static bool operator<=(const TransferDescriptor<Type> & desc, const uint val)
{
    return desc.data->rank <= val;
}

// This is required for the binary search routine
template<typename Type>
static bool operator>(const TransferDescriptor<Type> & desc, const uint val)
{
    return desc.data->rank > val;
}

//==============================================================================
//  Add_descriptor
//==============================================================================
template<Dim dim, Order ord, typename Type>
void Domain::Add_descriptor(const uint pos, const uint rank)
{
    // Create a new descriptor at the end of the array
    transfer_descriptors.Append();

    if (pos)
    {
        for (int d = transfer_descriptors.Get_size() - 1; d > pos; --d)
        {
            transfer_descriptors[d] = transfer_descriptors[d - 1];
        }
    }

    transfer_descriptors[pos].Set_data(rank);
}

//==============================================================================
//  Append_buffer_markings
//==============================================================================
template<Dim dim, Order ord, typename Type>
void Domain::Append_buffer_markings(
    const uint pos,
    const uint send_index,
    const uint recv_index,
    const uint size
)
{
    Array<BufferMarkings> & markings
        = transfer_descriptors[pos].Get_buffer_markings();

    uint offset = pos? markings.End().offset: 0;

    markings.Append();
    markings.End().Set(send_index, recv_index, size, offset + size);
}

//==============================================================================
//  Identify transfer descriptors
//==============================================================================
template<Dim dim, Order ord, typename Type>
void Domain::Identify_transfer_descriptors()
{
    uint trial_rank;
    uint trial_pos;

    // For each patch in the domain
    for (int p = 0; p < patches.Get_size(); ++p)
    {
        const Array<GhostMarkings> & markings = patches[p].Get_ghost_markings();

        // Examine each ghost marking
        for (int m = 0; m < markings.Get_size(); ++m)
        {
            // If marking does not correspond to patch itself
            if (m != markings[m].index)
            {
                // Identify trial MPI rank associated with the marking
                trial_rank = Binary_search(domain_bounds, markings.index);

                // Identify the correct position in transfer descriptor array
                trial_pos = Binary_search(transfer_descriptors, trial_rank);

                // If the trial rank descriptor is not found the array
                if (
                    !transfer_descriptors.Get_size()
                    || trial_rank != transfer_descriptors[trial_pos].Get_rank()
                )
                {
                    // Insert a new descriptor at the given position with rank
                    Add_descriptor(trial_pos, trial_rank);
                }

                Append_buffer_markings(
                    trial_pos, patches[p].Get_index(), markings[m].index,
                    Product<dim>(markings[m].sizes)
                );
            }
        }
    }

    for (int d = 0; d < transfer_descriptors.Get_size(); ++d)
    {
        for (int b = 0; b < transfer_descriptors[s].Get_buffer_markings().Get_size(); ++b)
        {
        }
    }

    // Shrink to fit
    transfer_descriptors.Truncate();
}

//==============================================================================
//  Constructor
//==============================================================================
template<Dim dim, Order ord, typename Type>
Domain::Domain(const Tuple<dim> & grid_sizs, const Tuple<dim> & patch_sizs):
    rank(0),
    range(0),
    domain_bounds(),
    grid_sizes(grid_sizs),
    patch_sizes(patch_sizs),
    patches(),
    local_markings(),
    transfer_descriptors()
{
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &range);

    Initialize_domain_bounds();

    patches.Reallocate(Get_patch_count());
    /// TODO /// Patch allocation

    if (range)
    {
        Identify_neighbour_ranks();
        Identify_comm_descriptors();
    }
}

//==============================================================================
//  Deallocation
//==============================================================================
template<Dim dim, Order ord, typename Type>
Domain::~Domain()
{
    for (uint p = 0; p < patches.Get_size(); ++p) { delete patches[p]; }

    for (uint t = 0; t < transfer_descriptors.Get_size(); ++t)
    {
        delete transfer_descriptors[t];
    }
}

//==============================================================================
//  Get patches index range
//==============================================================================
uint Get_patch_min_index(const uint ind) const { return domain_bounds[ind]; }

uint Get_patch_max_index(const uint ind) const
{
    return domain_bounds[ind + 1];
}

//==============================================================================
//  Get patch count
//==============================================================================
uint Get_patch_count(const uint ind) const
{
    return Get_patch_max_index(ind) - Get_patch_min_index(ind);
}

} // namespace maxwell

#endif // DOMAIN_HPP
