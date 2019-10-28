#ifndef DOMAIN_H
#define DOMAIN_H

#include "patch.h"
#include "utility.h"
#include "types.h"

namespace maxwell {

////////////////////////////////////////////////////////////////////////////////
//
//  BufferMarking
//
////////////////////////////////////////////////////////////////////////////////
// Contains buffer marking data corresponding to a single ghost transfer
////////////////////////////////////////////////////////////////////////////////
struct BufferMarking
{
    // Transfer data size
    uint size;

    // Transfer data offset in buffer
    uint offset;

    uint local_patch_index;
    uint extern_patch_index;

    // Indices of corresponding ghosts in a patch
    uint8_t local_ghost_index;
    uint8_t extern_ghost_index;

    BufferMarking();

    BufferMarking(
        const uint,
        const uint,
        const uint,
        const uint,
        const uint8_t,
        const uint8_t
    );
};

////////////////////////////////////////////////////////////////////////////////
//
//  TransferDescriptorData
//
////////////////////////////////////////////////////////////////////////////////
// Contains process communication descriptor data
// -----------------------------------------------------------------------------
// Template parameter:
//     Type -- the supported arithmetic type
////////////////////////////////////////////////////////////////////////////////
template<typename Type>
struct TransferDescriptorData
{
    //==========================================================================
    //  Data
    //==========================================================================
    // Extern MPI rank to communicate with
    uint rank;

    // Buffer markings
    Array<BufferMarking> buffer_markings;

    // Size of buffer
    uint buffer_size;

    // __device__
    // Buffers combined location, direct pointer to the sending buffer
    Type * sending_buffer;

    // __device__
    // Direct pointer to the receiving buffer
    Type * receiving_buffer;

    //==========================================================================
    //  Data management
    //==========================================================================
    TransferDescriptorData(const uint);
    ~TransferDescriptorData();
};

////////////////////////////////////////////////////////////////////////////////
//
//  TransferDescriptor
//
////////////////////////////////////////////////////////////////////////////////
// Implements a process communication descriptor
// -----------------------------------------------------------------------------
// Template parameter:
//     Type -- the supported arithmetic type
//////////////////////////////////////////////////////////////////////////////// 
template<typename Type>
class TransferDescriptor
{
    //==========================================================================
    //  Friends functions
    //==========================================================================
    friend bool operator<=(const TransferDescriptor &, const uint);
    friend bool operator>(const TransferDescriptor &, const uint);

    private:

        TransferDescriptorData<Type> * data;

    public:

        //======================================================================
        //  Data management
        //======================================================================
        void Set_data(const uint);
        void Allocate_buffers();

        void Add_buffer_marking(const BufferMarking &);

        TransferDescriptor(): data(NULL) {}
        ~TransferDescriptor() { delete data; }

        //======================================================================
        //  Access / mutate methods
        //======================================================================
        inline uint Get_rank() const { return data? data->rank: ~0; }

        inline uint Get_buffer_size() const { return data->buffer_size; }
        void Update_buffer_size(const uint);

        const Array<BufferMarking> & Get_buffer_markings() const;

        // Returns a raw pointer to the buffer location corresponding
        // to the record with a given index
        Type * Sending_buffer_record(const uint = 0);
        Type * Receiving_buffer_record(const uint = 0);

        //======================================================================
        //  Communication methods
        //======================================================================
        void Send() const;
        void Receive();
};

////////////////////////////////////////////////////////////////////////////////
//
//  Domain
//
////////////////////////////////////////////////////////////////////////////////
// Implements numerical grid data structure per process corresponding to
// a single MPI rank
// -----------------------------------------------------------------------------
// Template parameters:
//     dim -- the dimensionality of the grid
//     ord -- the order of patch numeration
//     Type -- the supported arithmetic type
// -----------------------------------------------------------------------------
/// !!! /// Implementation should probably be changed to a singleton class
// -----------------------------------------------------------------------------
/// !!! /// Type should probably be changed to some data structure, representing
/// !!! /// fields and particles
////////////////////////////////////////////////////////////////////////////////
template<Dimension dim, Order ord, typename Type>
class Domain
{
    private:

        //======================================================================
        //  Geometry descriptors
        //======================================================================
        // MPI rank
        uint rank;

        // Global range of ranks
        uint range;

        // Grid Hilbert index decomposition
        Array<uint> domain_bounds;

        // Grid sizes
        Tuple<dim> grid_sizes;

        // Patch sizes
        Tuple<dim> patch_sizes;

        //======================================================================
        //  Grid data in patches
        //======================================================================
        // Patches arranged in specified Order
        Array< Patch<dim, ord, Type> > patches;

        //======================================================================
        //  Intra-domain communication markings
        //======================================================================
        Array<BufferMarking> local_markings;

        //======================================================================
        //  Inter-domain communication descriptors
        //======================================================================
        Array<TransferDescriptor> transfer_descriptors;

        //======================================================================
        //  Initialization methods
        //======================================================================
        // Sets initial domain bounds
        void Initialize_domain_bounds();

        // Collects the domain bounds from all the processes
        void Set_domain_bounds();

        //======================================================================
        //  Intra-domain communication methods
        //======================================================================
        void Perform_local_transfers()

        //======================================================================
        //  Inter-domain communication methods
        //======================================================================
        // Inserts a new transfer descriptor corresponding to the given MPI rank
        // at a given position to the descriptor array already sorted
        // in ascending order
        void Create_transfer_descriptor(const uint, const uint);

        // Returns a transfer descriptor corresponding to a target MPI rank
        TransferDescriptor<Type> & Get_transfer_descriptor(const uint);

        // Creates the buffer marking in the correct transfer descriptor array
        // corrseponding to a given sender date and a ghost marking
        void Create_buffer_marking(
            const uint, const uint8_t, const GhostMarking &
        );

        // Allocates inter-domain transfer buffers on device 
        void Allocate_transfer_buffers();

        void Set_transfer_descriptors();

        // Prepare a transfer package to send
        void Pack_transfer_data();

        // Unpack a received transfer package
        void Unpack_transfer_data();

        /// FIXME /// Either remove or switch to c++11: = delete
        Domain() {}

    public:

        //======================================================================
        //  Data management
        //======================================================================
        /// FIXME /// Patch sizes, width etc. should probably be included here
        Domain(const Tuple<dim> &, const Tuple<dim> &);

        ~Domain() {}

        //======================================================================
        //  Access / mutate methods
        //======================================================================
        // Returns MPI rank for the domain
        inline uint Get_rank() const { return rank; }

        // Returns patch indices range for a given domain
        uint Get_patch_min_index(const uint = rank) const;
        uint Get_patch_max_index(const uint = rank) const;

        // Computes patch count for a given domain
        uint Get_patch_count(const uint = rank) const;
};

} // namespace maxwell

#include "domain.hpp"

#endif // DOMAIN_H
