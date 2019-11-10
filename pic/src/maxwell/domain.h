#ifndef DOMAIN_H
#define DOMAIN_H

#include "patch.h"
#include "utility.h"
#include "types.h"

namespace maxwell {

////////////////////////////////////////////////////////////////////////////////
//
//  TransferMarking
//
////////////////////////////////////////////////////////////////////////////////
// Contains marking data corresponding to a transfer of a single ghost
////////////////////////////////////////////////////////////////////////////////
struct TransferMarking
{
    // Buffer record size
    uint size;

    // Buffer record offset
    uint offset;

    uint sending_patch_index;
    uint receiving_patch_index;

    // Indices of corresponding ghosts in a patch
    uint8_t sending_ghost_index;
    uint8_t receiving_ghost_index;

    //==========================================================================
    //  Data management
    //==========================================================================
    TransferMarking();

    TransferMarking(
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
    Array<TransferMarking> transfer_markings;

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

        //======================================================================
        //  Access / mutate methods
        //======================================================================
        // Returns a raw pointer to the buffer location corresponding
        // to the record with a given index
        Type * Sending_buffer_record(const uint = 0);
        Type * Receiving_buffer_record(const uint = 0);

    public:

        //======================================================================
        //  Data management
        //======================================================================
        void Set(const uint);
        void Allocate_buffers();
        void Add_transfer_marking(const TransferMarking &);

        TransferDescriptor(): data(NULL) {}
        ~TransferDescriptor() { delete data; }

        //======================================================================
        //  Access / mutate methods
        //======================================================================
        inline uint Get_rank() const { return data? data->rank: ~0; }

        inline uint Get_buffer_size() const { return data->buffer_size; }
        void Update_buffer_size(const uint);

        const Array<TransferMarking> & Get_transfer_markings() const;

        //======================================================================
        //  Communication methods
        //======================================================================
        template<Dimension dim, Order ord>
        void Pack_transfer_data(const Array< Patch<dim, ord, Type> > &);

        template<Dimension dim, Order ord>
        void Unpack_transfer_data(Array< Patch<dim, ord, Type> > &) const;

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

        // Grid domain decomposition bounds
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
        Array<TransferMarking> local_markings;

        //======================================================================
        //  Inter-domain communication descriptors
        //======================================================================
        Array<TransferDescriptor> transfer_descriptors;

        //======================================================================
        //  Data management
        //======================================================================
        /// FIXME /// Either remove or switch to c++11: = delete
        // Domain() {}

        // Sets initial domain bounds
        void Initialize_domain_bounds();

        // Collects the domain bounds from all the processes
        void Set_domain_bounds();

        void Allocate_patches();

        //======================================================================
        //  Intra-domain communication methods
        //======================================================================
        void Perform_local_transfers();

        //======================================================================
        //  Inter-domain communication methods
        //======================================================================
        // Inserts a new transfer descriptor corresponding to a given MPI rank
        // at a given position to the descriptor array already sorted
        // in the ascending order
        void Create_transfer_descriptor(const uint, const uint);

        // Returns a transfer descriptor corresponding to a target MPI rank
        TransferDescriptor<Type> & Transfer_descriptor(const uint);

        // Creates the transfer marking in the correct transfer descriptor array
        // corrseponding to a given sender date and a ghost marking
        void Create_transfer_marking(
            const uint, const uint8_t, const GhostMarking<dim, Type> &
        );

        // Allocates inter-domain transfer buffers on device 
        void Allocate_transfer_buffers();

        // Configures inter-domain communications
        void Set_transfer_markings();

        void Perform_global_transfers();

    public:

        //======================================================================
        //  Data management
        //======================================================================
        Domain(const Tuple<dim> &, const Tuple<dim> &, const uint);
        ~Domain() {}

        //======================================================================
        //  Access / mutate methods
        //======================================================================
        // Returns MPI rank for the domain
        inline uint Get_rank() const { return rank; }

        // Returns patch indices range for a given domain
        uint Get_min_patch_index(const uint = rank) const;
        uint Get_max_patch_index(const uint = rank) const;

        // Computes patch count for a given domain
        uint Get_patch_count(const uint = rank) const;
};

} // namespace maxwell

#include "domain.hpp"

#endif // DOMAIN_H
