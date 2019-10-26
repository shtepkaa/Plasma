#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include <mpi.h>

namespace maxwell {

////////////////////////////////////////////////////////////////////////////////
//
//  TransferDescriptorData
//
////////////////////////////////////////////////////////////////////////////////
TransferDescriptorData::TransferDescriptorData(const uint target_rank):
    rank(target_rank),
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
/// ??? /// TransferDescriptor::TransferDescriptor(const uint target_rank):
/// ??? ///     data(new TransferDescriptorData(target_rank))
/// ??? /// {}

//==============================================================================
//  Set_data
//==============================================================================
void TransferDescriptor::Set_data(const uint target_rank)
{
    data = new TransferDescriptorData(target_rank);
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

    if (ord == CARTESIAN)
    {
        /// TODO ///
    }
    else if (ord == HILBERTIAN)
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
//  Create_transfer_descriptor
//==============================================================================
template<Dim dim, Order ord, typename Type>
void Domain::Create_transfer_descriptor(const uint ind, const uint rank)
{
    // Create a new descriptor at the end of the array
    transfer_descriptors.Append();

    // Shift the right part of the array
    for (int t = transfer_descriptors.Get_size() - 1; t > ind; --t)
    {
        transfer_descriptors[t] = transfer_descriptors[t - 1];
    }

    transfer_descriptors[ind].Set_data(rank);
}

//==============================================================================
//  Get_transfer_descriptor
//==============================================================================
// Finds the corresponding descriptor in array by rank comparisons, creates
// a new descriptor if there is no an appropriate one
template<Dim dim, Order ord, typename Type>
TransferDescriptor<Type> & Domain::Get_transfer_descriptor(
    const uint target_rank
)
{
    // Identify the correct position in transfer descriptor array
    const uint ind
        = transfer_descriptors.Get_size()?
            Find_index(transfer_descriptors, target_rank):
            0;

    // If the receiver rank descriptor is not found in the array
    if (
        !transfer_descriptors.Get_size()
        || target_rank != transfer_descriptors[ind].Get_rank()
    )
    {
        // Insert a new descriptor at the given index
        Create_transfer_descriptor(ind, target_rank);
    }

    return transfer_descriptors[ind];
}

//==============================================================================
//  Create_buffer_marking
//==============================================================================
template<Dim dim, Order ord, typename Type>
void Domain::Create_buffer_marking(
    const uint local_patch_index,
    const uint8_t local_ghost_index,
    const GhostMarking & ghost_marking
)
{
    const uint extern_patch_index = ghost_marking.target_patch_index;

    // If ghost marking does not correspond to the patch itself
    if (local_patch_index != extern_patch_index)
    {
        const uint size = Product<dim>(ghost_marking.sizes);
        const uint8_t extern_ghost_index = ghost_marking.target_ghost_index;

        // Identify receiving MPI rank associated with the ghost marking
        const uint target_rank
            = domain_bounds[Find_index(domain_bounds, extern_patch_index)];

        // If the transfer is local
        if (target_rank == rank)
        {
            // Construct a new marking at the end of the local transfers array
            local_markings.Append(
                BufferMarking(
                    size, UNDEFINED, local_patch_index, extern_patch_index,
                    local_ghost_index, extern_ghost_index
                )
            );
        }
        // Otherwise, the transfer is an interprocess one
        else
        {
            TransferDescriptor<Type> & desc
                = Get_transfer_descriptor(target_rank);

            // Construct a new marking at the end of corresponding array
            desc.Get_buffer_markings().Append(
                BufferMarking(
                    size, desc.Get_size(), local_patch_index,
                    extern_patch_index, local_ghost_index, extern_ghost_index
                )
            );

            // Enlarge the size of the buffer by the one of the ghost
            desc.Update_size(size);
        }
    }
}

//==============================================================================
//  Allocate_transfer_buffers 
//==============================================================================
template<Dim dim, Order ord, typename Type>
void Allocate_transfer_buffers()
{
    for (int t = 0; t < transfer_descriptors.Get_size(); ++t)
    {
        /// TODO /// CudaMalloc here
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
        const uint local_patch_index = patches[p].Get_index();

        const Array<GhostMarking> & ghost_markings
            = patches[p].Get_ghost_markings();

        // For each ghost marking
        for (int g = 0; g < ghost_markings_size; ++g)
        {
            Create_buffer_marking(local_patch_index, g, ghost_markings[g]);
        }
    }

    // Allocate buffers here
    Allocate_transfer_buffers();

    /// FIXME /// Also consider local transfers here, like extern ones below

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
                buffer_marking.local_patch_index - Get_patch_min_index()
            ].Copy_ghost(buffer_marking.local_ghost_index, );
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

    /// TODO /// Patch allocation
    /// ??? /// patches.Reallocate(Get_patch_count());

    if (range)
    {
        /// ??? /// Identify_neighbour_ranks();
        Identify_transfar_descriptors();
    }
}

//==============================================================================
//  Deallocation
//==============================================================================
template<Dim dim, Order ord, typename Type>
Domain::~Domain()
{
    /// ??? /// for (uint p = 0; p < patches.Get_size(); ++p) { delete patches[p]; }
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
