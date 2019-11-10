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
    transfer_markings(),
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
//  Access buffer records
//==============================================================================
template<typename Type>
Type * TransferDescriptor::Sending_buffer_record(const uint ind)
{
    return data->sending_buffer + data->transfer_markings[ind].offset;
}

template<typename Type>
Type * TransferDescriptor::Receiving_buffer_record(const uint ind)
{
    return data->receiving_buffer + data->transfer_markings[ind].offset;
}

//==============================================================================
//  Set
//==============================================================================
template<typename Type>
void TransferDescriptor::Set(const uint target_rank)
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
//  Add transfer marking
//==============================================================================
void TransferDescriptor::Add_transfer_marking(const TransferMarking & marking)
{
    // Enlarge the size of the buffer by the one of the marking
    Update_buffer_size(marking.size);

    transfer_markings.Append(marking);
}

//==============================================================================
//  Get transfer markings
//==============================================================================
template<typename Type>
const Array<TransferMarking> & TransferDescriptor::Get_transfer_markings() const
{
    return data->transfer_markings;
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
//  Pack transfer data
//==============================================================================
template<typename Type>
template<Dimension dim, Order ord>
void TransferDescriptor::Pack_transfer_data(
    const Array< Patch<dim, ord, Type> > & patches
)
{
    uint min_ind = patches[0].Get_index();

    // For each transfer marking
    for (uint b = 0; b < transfer_markings.Get_size(); ++b)
    {
        const uint ind = transfer_markings[b].sending_patch_index - min_ind;

        // Send the ghost to the appropriate buffer record
        patches[ind].Send_ghost(
            transfer_markings[b].sending_ghost_index, Sending_buffer_record(b)
        );
    }
}

//==============================================================================
//  Unpack transfer data
//==============================================================================
template<typename Type>
template<Dimension dim, Order ord>
void TransferDescriptor::Unpack_transfer_data(
    Array< Patch<dim, ord, Type> > & patches
) const
{
    uint min_ind = patches[0].Get_index();

    // For each transfer marking
    for (uint b = 0; b < transfer_markings.Get_size(); ++b)
    {
        const uint ind = transfer_markings[b].sending_patch_index - min_ind;

        // Receive the ghost from the correct buffer record
        patches[ind].Receive_ghost(
            transfer_markings[b].receiving_ghost_index,
            Receiving_buffer_record(b)
        );
    }
}

//==============================================================================
//  Send
//==============================================================================
/// FIXME /// Add MPI datatype
/// FIXME /// Add request handling
/// FIXME /// Add RDMA checking
template<typename Type>
void TransferDescriptor::Send() const
{
    MPI_Request request;

    /* int */ MPI_Isend(
        (void *)data->sending_buffer,
        data->size,
            MPI_Datatype, /// Create MPI datatype
                          /// for sending structure somewhere else
        data->rank, // receiving rank
        rank, // tag coincides with sending rank
        MPI_COMM_WORLD,
        request
    );
}

//==============================================================================
//  Receive
//==============================================================================
/// FIXME /// Add MPI datatype
/// FIXME /// Add status handling
/// FIXME /// Add RDMA checking
template<typename Type>
void TransferDescriptor::Receive()
{
    MPI_Status status;

    /* int */ MPI_Recv(
        (void *)data->receiving_buffer,
        data->size,
            MPI_Datatype,
        data->rank, // sending rank
        data->rank, // tag coincides with sending rank
        MPI_COMM_WORLD,
        &status
    );
}

////////////////////////////////////////////////////////////////////////////////
//
//  Domain
//
////////////////////////////////////////////////////////////////////////////////
//==============================================================================
//  Initialize domain bounds
//==============================================================================
/// FIXME /// Cartesian index
template<Dimension dim, Order ord, typename Type>
void Domain::Initialize_domain_bounds()
{
    domain_bounds.Allocate(range + 1);
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
        (void *)(domain_bounds + rank), 1, MPI_UNSIGNED,
        (void *)domain_bounds, 1, MPI_UNSIGNED,
        MPI_COMM_WORLD
    );
}

//==============================================================================
//  Allocate patches
//==============================================================================
template<Dimension dim, Order ord, typename Type>
void Allocate_patches(const uint ghost_width);
{
    Tuple<dim> layer_sizes = grid_sizes / patch_sizes;

    patches.Allocate(Get_patch_count());

// #pragma omp parallel for
    for (uint p = 0; p < patches.Get_size(); ++p)
    {
        uint patch_index = Get_min_patch_index() + p;

        const Tuple<dim> patch_coords
            = Compute_patch_coordinates<dim, ord>(layer_sizes, patch_index);

        patches[p].Set(
            layer_sizes, patch_coords, patch_sizes, ghost_width, patch_index
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
    for (uint t = transfer_descriptors.Get_size() - 1; t > ind; --t)
    {
        transfer_descriptors[t] = transfer_descriptors[t - 1];
    }

    // Set a new descriptor at the correct position
    transfer_descriptors[ind].Set(rank);
}

//==============================================================================
//  Transfer descriptor
//==============================================================================
// Finds the corresponding descriptor in array by rank comparisons, creates
// a new descriptor if there is no an appropriate one
template<Dimension dim, Order ord, typename Type>
TransferDescriptor<Type> & Domain::Transfer_descriptor(const uint target_rank)
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
//  Create transfer marking
//==============================================================================
template<Dimension dim, Order ord, typename Type>
void Domain::Create_transfer_marking(
    const uint sending_patch_index,
    const uint8_t sending_ghost_index,
    const GhostMarking<dim, Type> & ghost_marking
)
{
    const uint size = Product<dim>(ghost_marking.sizes);

    const uint receiving_patch_index = ghost_marking.target_patch_index;
    const uint8_t receiving_ghost_index = ghost_marking.target_ghost_index;

    // Identify receiving MPI rank associated with the ghost marking
    const uint target_rank
        = domain_bounds[Find_index(domain_bounds, receiving_patch_index)];

    // If the transfer is local (i. e. intra-domain)
    if (target_rank == rank)
    {
        // Construct a new marking at the end of the local transfers array
        local_markings.Append(
            TransferMarking(
                size, UNDEFINED,
                sending_patch_index, receiving_patch_index,
                sending_ghost_index, receiving_ghost_index
            )
        );
    }
    // Otherwise, the transfer is an inter-domain one
    else
    {
        // Find the correct descriptor, insert a new one if needed
        TransferDescriptor<Type> & desc = Transfer_descriptor(target_rank);

        // Construct a new marking at the end of corresponding array
        desc.Add_transfer_marking(
            TransferMarking(
                size, desc.Get_buffer_size(),
                sending_patch_index, receiving_patch_index,
                sending_ghost_index, receiving_ghost_index
            )
        );
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
//  Set transfer markings 
//==============================================================================
template<Dimension dim, Order ord, typename Type>
void Domain::Set_transfer_markings()
{
    // For each patch in the domain
    for (uint p = 0; p < patches.Get_size(); ++p)
    {
        const uint sending_patch_index = patches[p].Get_index();

        const Array< GhostMarking<dim, Type> > & ghost_markings
            = patches[p].Get_ghost_markings();

        // For each ghost marking
        for (uint g = 0; g < ghost_markings.Get_size(); ++g)
        {
            const uint receiving_patch_index
                = ghost_markings[g].target_patch_index;

            // If ghost marking does not correspond to the patch itself
            if (
                receiving_patch_index != ~0
                && sending_patch_index != receiving_patch_index
            )
            {
                Create_transfer_marking(
                    sending_patch_index, g, ghost_markings[g]
                );
            }
        }
    }

    // Shrink to fit
    transfer_descriptors.Truncate();

    Allocate_transfer_buffers();
}

//==============================================================================
//  Perform local transfers
//==============================================================================
template<Dimension dim, Order ord, typename Type>
void Domain::Perform_local_transfers()
{
// #pragma omp parallel for
    for (uint l = 0; l < local_markings; ++l)
    {
        const TransferMarking & marking = local_markings[l];

        // Transfer ghost
        patches[marking.sending_patch_index].Transfer_ghost(
            marking.sending_ghost_index, patches[marking.receiving_patch_index],
            marking.receiving_ghost_index
        );
    }
}

//==============================================================================
//  Perform global transfers
//==============================================================================
template<Dimension dim, Order ord, typename Type>
void Domain::Perform_global_transfers()
{
    for (uint t = 0; t < transfer_descriptors.Get_size(); ++t)
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
    Set_transfer_markings();
}

//==============================================================================
//  Get patches index range
//==============================================================================
// [unsafe]
uint Get_min_patch_index(const uint ind) const { return domain_bounds[ind]; }

// [unsafe]
uint Get_max_patch_index(const uint ind) const
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
