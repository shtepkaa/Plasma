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
// Implements buffer marking 
////////////////////////////////////////////////////////////////////////////////
/// ??? /// template<Dim dim>
struct BufferMarking
{
    uint send_index;
    uint recv_index;

    /// Tuple<dim> sizes;

    uint size;
    uint offset;
};

////////////////////////////////////////////////////////////////////////////////
//
//  TransferDescriptorData
//
////////////////////////////////////////////////////////////////////////////////
// Implements process communication descriptor data
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
    // Communication rank
    uint rank;

    // Buffers marking
    Array<BufferMarking> buffer_markings;

    // __device__
    // Incoming buffer
    Type * send_buffer;

    // __device__
    // Outcoming buffer
    Type * recv_buffer;

    //==========================================================================
    //  Data management
    //==========================================================================
    TransferDescriptorData(const uint, const uint);
    ~TransferDescriptorData();
};

////////////////////////////////////////////////////////////////////////////////
//
//  TransferDescriptor
//
////////////////////////////////////////////////////////////////////////////////
// Implements process communication descriptor
// -----------------------------------------------------------------------------  
// Template parameter:
//     Type -- the supported arithmetic type
////////////////////////////////////////////////////////////////////////////////
template<typename Type>
class TransferDescriptor
{
    private:

        TransferDescriptorData<Type> * data;

    public:

        //======================================================================
        //  Data management
        //======================================================================
        TransferDescriptor(const uint, const uint);
        ~TransferDescriptor() { delete data; }
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
template<Dim dim, Order ord, typename Type>
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
        Array< Patch<dim, ord, Type> * > patches;

        //======================================================================
        //  Intra-domain communication descriptors
        //======================================================================
        // Local buffer marking
        Array<GhostMarking> local_markings;

        //======================================================================
        //  Inter-domain communication descriptors
        //======================================================================
        Array<TransferDescriptor> transfer_descriptors;

        //======================================================================
        //  Initialization methods
        //======================================================================
        // Sets initial domain bounds
        void Initialize_domain_bounds();

        // Identifies domain bounds
        void Identify_domain_bounds();

            /// // Allocates patches
            /// void Allocate_patches();

        // Identifies transfer descriptors
        void Identify_transfer_descriptors();

        // Default
        /// FIXME /// Either remove or switch to c++11: = delete
        Domain() {}

    public:

        //======================================================================
        //  Data management
        //======================================================================
        // Construction
        /// FIXME /// Data initialization required
        Domain(const Tuple<dim> &, const Tuple<dim> &);

        // Deallocation
        ~Domain();

        //======================================================================
        //  Get methods
        //======================================================================
        // Get rank
        inline uint Get_rank() const { return rank; }

        // Get patches index range
        inline uint Get_patch_min_index() const { return domain_bounds[rank]; }
        uint Get_patch_max_index() const { return domain_bounds[rank + 1]; }

        //======================================================================
        //  Get methods
        //======================================================================
        // Compute patch count
        inline uint Get_patch_count() const { return patches.Get_size(); }

            /// // copy patch to CPU
            /// void Copy_patch(const uint, const Type *) const;

};

} // namespace maxwell

#include "domain.hpp"

#endif // DOMAIN_H
