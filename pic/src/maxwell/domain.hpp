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
//  Set_data
//==============================================================================
void TransferDescriptor::Set_data(const uint transfer_rank)
{
    data = new TransferDescriptorData(transfer_rank);
}

//==============================================================================
//  Get_buffer_markings
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
//  Initialize_domain_bounds
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
//  Identify_domain_bounds
//==============================================================================
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
//  Add_transfer_descriptor
//==============================================================================
template<Dim dim, Order ord, typename Type>
void Domain::Add_transfer_descriptor(const uint pos, const uint rank)
{
    // Create a new descriptor at the end of the array
    transfer_descriptors.Append();

    if (pos)
    {
        for (int t = transfer_descriptors.Get_size() - 1; t > pos; --t)
        {
            transfer_descriptors[t] = transfer_descriptors[t - 1];
        }
    }

    transfer_descriptors[pos].Set_data(rank);
}

//==============================================================================
//  Find_transfer_descriptor_index
//==============================================================================
template<Dim dim, Order ord, typename Type>
uint Domain::Find_transfer_descriptor_index(const GhostMarking & ghost_marking)
{
    // Identify transfer MPI rank associated with the ghost marking
    const uint transfer_rank
        = domain_bounds[Find_index(domain_bounds, ghost_marking.patch_index)];

    if (transfer_rank == rank) { return ~0; }
    else
    {
        // Identify the correct position in transfer descriptor array
        const uint ind = Find_index(transfer_descriptors, transfer_rank);

        // If the transfer rank descriptor is not found in the array
        if (
            !transfer_descriptors.Get_size()
            || transfer_rank != transfer_descriptors[ind].Get_rank()
        )
        {
            // Insert a new descriptor at the given position with the new rank
            Add_transfer_descriptor(ind, transfer_rank);
        }

        return ind;
    }
}

//==============================================================================
//  Identify_transfer_descriptors 
//==============================================================================
template<Dim dim, Order ord, typename Type>
void Domain::Identify_transfer_descriptors()
{
    const uint ghost_markings_size = patch[0].Get_ghost_markings().Get_size();

    // For each patch in the domain
    for (int p = 0; p < patches.Get_size(); ++p)
    {
        const uint send_patch_index = patches[p].Get_index();

        // Examine each ghost marking
        for (int g = 0; g < ghost_markings_size; ++g)
        {
            const GhostMarking & ghost_marking
                = patches[p].Get_ghost_markings()[g];

            const uint recv_patch_index = ghost_marking.patch_index;

            // If marking does not correspond to patch itself
            if (send_patch_index != recv_patch_index)
            {
                // Find a buffer markings array of a transfer descriptor
                // correspoding to a receiving patch index
                Array<BufferMarking> & buffer_markings
                    = Find_buffer_markings(ghost_marking);

                const uint size = Product<dim>(ghost_marking.sizes);

                const uint offset
                    = buffer_markings.Get_size()?
                        buffer_markings.End().offset:
                        0;

                // Add new marking
                buffer_markings.Append(
                    BufferMarking(
                        size, offset + size, send_patch_index, recv_patch_index,
                        g, ghost_marking.recv_ghost_index
                    )
                );
            }
        }
    }

    /// FIXME /// Allocate buffers here

    /// FIXME /// Also consider local transfers here

    for (int t = 0; t < transfer_descriptors.Get_size(); ++t)
    {
        for (
            int b = 0;
            b < transfer_descriptors[t].Get_buffer_markings().Get_size();
            ++b
        )
        {
            const BufferMarking & buffer_marking
                = transfer_descriptors[t].Get_buffer_markings()[b];

            patches[
                buffer_marking.send_patch_index - Get_patch_min_index()
            ].Copy_ghost(buffer_marking.send_ghost_index, );
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
// [unsafe]
uint Get_patch_min_index(const uint ind) const { return domain_bounds[ind]; }

// [unsafe]
uint Get_patch_max_index(const uint ind) const
{
    return domain_bounds[ind + 1];
}

//==============================================================================
//  Get_patch_count
//==============================================================================
// [unsafe]
uint Get_patch_count(const uint ind) const
{
    return domain_bounds[ind + 1] - domain_bounds[ind];
}

} // namespace maxwell

#endif // DOMAIN_HPP
