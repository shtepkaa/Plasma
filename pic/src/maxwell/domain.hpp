#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include <mpi.h>

namespace maxwell {

////////////////////////////////////////////////////////////////////////////////
//
//  TransferDescriptorData
//
////////////////////////////////////////////////////////////////////////////////
template<typename Type>
TransferDescriptorData::TransferDescriptorData(const uint target_rank):
    rank(target_rank),
    buffer_markings(),
    size(0),
    sending_buffer(NULL),
    receiving_buffer(NULL)
{}

template<typename Type>
TransferDescriptorData::~TransferDescriptorData()
{
    if (sending_buffer) { CUDA_CALL(cudaFree(sending_buffer)); }
}

////////////////////////////////////////////////////////////////////////////////
//
//  TransferDescriptor
//
////////////////////////////////////////////////////////////////////////////////
//==============================================================================
//  Set data
//==============================================================================
template<typename Type>
void TransferDescriptor::Set_data(const uint target_rank)
{
    data = new TransferDescriptorData(target_rank);
}

//==============================================================================
//  Allocate buffers
//==============================================================================
template<typename Type>
void TransferDescriptor::Allocate_buffers()
{
    CUDA_CALL(cudaMalloc(sending_buffer, 2 * buffer_size * sizeof(Type)));

    receiving_buffer = sending_buffer + buffer_size;
}

//==============================================================================
//  Add buffer marking
//==============================================================================
void Add_buffer_marking(const BufferMarking & marking)
{
    buffer_markings.Append(marking);
}

//==============================================================================
//  Get buffer markings
//==============================================================================
template<typename Type>
const Array<BufferMarking> & TransferDescriptor::Get_buffer_markings() const
{
    return data->buffer_markings;
}

//==============================================================================
//  Update buffer size
//==============================================================================
template<typename Type>
void TransferDescriptor::Update_buffer_size(const uint marking_size)
{
    data->buffer_size += marking_size;
}

//==============================================================================
//  Access buffer records
//==============================================================================
template<typename Type>
Type * TransferDescriptor::Sending_buffer_record(const uint ind)
{
    return data->sending_buffer + data->buffer_markings[ind].offset;
}

template<typename Type>
Type * TransferDescriptor::Receiving_buffer_record(const uint ind)
{
    return data->receiving_buffer + data->buffer_markings[ind].offset;
}

//==============================================================================
//  Pack transfer data
//==============================================================================
template<typename Type>
void TransferDescriptor::Pack_transfer_data(const Array<Patch> & patches)
{
    uint min_ind = patches[0].Get_index();

    // For each buffer marking
    for (int b = 0; b < buffer_markings.Get_size(); ++b)
    {
        const uint ind = buffer_markings[b].local_patch_index - min_ind;

        // Copy the associated ghost to the correct sending buffer record
        patches[ind].Copy_ghost(
            buffer_markings[b].local_ghost_index, Sending_buffer_record(b)
        );
    }
}

//==============================================================================
//  Unpack transfer data
//==============================================================================
template<typename Type>
void TransferDescriptor::Unpack_transfer_data(
    const Array<Patch> & patches
) const
{
    uint min_ind = patches[0].Get_index();

    // For each buffer marking
    for (int b = 0; b < buffer_markings.Get_size(); ++b)
    {
        const uint ind = buffer_markings[b].local_patch_index - min_ind;

        // Copy the associated receiving buffer record to the correct ghost
        patches[ind].Set_ghost(
            buffer_markings[b].local_ghost_index, Receiving_buffer_record(b)
        );
    }
}

//==============================================================================
//  Send
//==============================================================================
template<typename Type>
void TransferDescriptor::Send() const
{
    /// MPI_Isend(
    ///     (void *)(domain_bounds + rank), 1, MPI_UNSIGNED, (void *)domain_bounds,
    ///     1, MPI_UNSIGNED, MPI_COMM_WORLD
    /// );

    MPI_Request request;

    int MPI_Isend(
        (void *)data->rankbuffer,
        data->rankize,
            MPI_Datatype, /// Create MPI datatype
                          /// for sending structure somewhere else
        data->rank,
            int tag, /// Create tag /// ??? send + recv + iteration_number
        MPI_COMM_WORLD,
        request
    );
}

//==============================================================================
//  Receive
//==============================================================================
template<typename Type>
void TransferDescriptor::Receive()
{
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
template<Dimension dim, Order ord, typename Type>
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

        uint patches_per_domain = total_patch_count / range;
        uint rem = total_patch_count % range;

        for (uint r = 1; r <= range; ++r)
        {
            domain_bounds[r]
                = domain_bounds[r - 1] + patches_per_domain + (r < rem);
        }
    }
}

//==============================================================================
//  Set domain bounds
//==============================================================================
template<Dimension dim, Order ord, typename Type>
void Domain::Set_domain_bounds()
{
    MPI_Allgather(
        (void *)(domain_bounds + rank), 1, MPI_UNSIGNED, (void *)domain_bounds,
        1, MPI_UNSIGNED, MPI_COMM_WORLD
    );
}

//==============================================================================
//  Allocate patches
//==============================================================================
template<Dimension dim, Order ord, typename Type>
void Allocate_patches(const uint ghost_width);
{
    patches.Reallocate(Get_patch_count());

    for (uint p = 0; p < patches.Get_size(); ++p)
    {
        /// Coordinates or index should be provided? (Coordinates probably)
        /// TO BE CONTINUED
        patches[p].Set_data(
            grid_sizes / patch_sizes,
                const Tuple<dim> & patch_coordinates,
            patch_sizes,
            ghost_width,
                const uint index
        );
    }
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
//  Create transfer descriptor
//==============================================================================
template<Dimension dim, Order ord, typename Type>
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
//  Get transfer descriptor
//==============================================================================
// Finds the corresponding descriptor in array by rank comparisons, creates
// a new descriptor if there is no an appropriate one
template<Dimension dim, Order ord, typename Type>
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
//  Create buffer marking
//==============================================================================
template<Dimension dim, Order ord, typename Type>
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
        // Otherwise, the transfer is an inter-domain one
        else
        {
            TransferDescriptor<Type> & desc
                = Get_transfer_descriptor(target_rank);

            // Construct a new marking at the end of corresponding array
            desc.Add_buffer_marking(
                BufferMarking(
                    size, desc.Get_buffer_size(), local_patch_index,
                    extern_patch_index, local_ghost_index, extern_ghost_index
                )
            );

            // Enlarge the size of the buffer by the one of the ghost
            desc.Update_buffer_size(size);
        }
    }
}

//==============================================================================
//  Allocate transfer buffers 
//==============================================================================
template<Dimension dim, Order ord, typename Type>
void Domain::Allocate_transfer_buffers()
{
    for (uint t = 0; t < transfer_descriptors.Get_size(); ++t)
    {
        transfer_descriptors[t].Allocate_buffers();
    }
}

//==============================================================================
//  Set transfer descriptors 
//==============================================================================
template<Dimension dim, Order ord, typename Type>
void Domain::Set_transfer_descriptors()
{
    // For each patch in the domain
    for (uint p = 0; p < patches.Get_size(); ++p)
    {
        const uint local_patch_index = patches[p].Get_index();

        const Array<GhostMarking> & ghost_markings
            = patches[p].Get_ghost_markings();

        // For each ghost marking
        for (uint g = 0; g < ghost_markings.Get_size(); ++g)
        {
            Create_buffer_marking(local_patch_index, g, ghost_markings[g]);
        }
    }

    // Allocate buffers here
    Allocate_transfer_buffers();

    // Shrink to fit
    transfer_descriptors.Truncate();
}

//==============================================================================
//  Perform local transfers
//==============================================================================
template<Dimension dim, Order ord, typename Type>
void Domain::Perform_local_transfers()
{
    for (uint ; ; )
    {
        /// FIXME /// Also consider local transfers here, like extern ones below
    }
}

//==============================================================================
//  Perform global transfers
//==============================================================================
template<Dimension dim, Order ord, typename Type>
void Domain::Perform_global_transfers()
{
    for (int t = 0; t < transfer_descriptors.Get_size(); ++t)
    {
        transfer_descriptors[t].Pack_transfer_data(patches);
        transfer_descriptors[t].Send();

        transfer_descriptors[t].Receive();
        transfer_descriptors[t].Unpack_transfer_data(patches);
    }
}

//==============================================================================
//  Data management
//==============================================================================
template<Dimension dim, Order ord, typename Type>
Domain::Domain(
    const Tuple<dim> & data_grid_sizes,
    const Tuple<dim> & data_patch_sizes,
    const uint ghost_width
):
    rank(0),
    range(0),
    domain_bounds(),
    grid_sizes(data_grid_sizes),
    patch_sizes(data_patch_sizes),
    patches(),
    local_markings(),
    transfer_descriptors()
{
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &range);

    Initialize_domain_bounds();

    Allocate_patches(ghost_width);

    // Configure communication descriptors in case of multi-process computation 
    if (range) { Set_transfer_descriptors(); }
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
//  Get patch count
//==============================================================================
// [unsafe]
uint Get_patch_count(const uint ind) const
{
    return domain_bounds[ind + 1] - domain_bounds[ind];
}

} // namespace maxwell

#endif // DOMAIN_HPP
