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
//  Domain
//
////////////////////////////////////////////////////////////////////////////////
// Implements data structure per process corresponding to a single MPI rank
// Template parameter: dim -- the dimensionality of the grid
// Template parameter: ord -- the order of patch numeration
// Template parameter: Type -- the supported arithmetic type
//
/// !!! /// Implementation should probably be changed to a singleton class
//
/// !!! /// Type should probably be changed to some data structure, representing
/// !!! /// fields and particles
template<Dim dim, Order ord, typename Type = double>
class Domain
{
    private:

        //====================================================================//
        //  Geometry descriptors
        //====================================================================//
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

        //====================================================================//
        //  Data
        //====================================================================//
        // Patches arranged in specified Order
        Array< Patch<dim, ord, Type> * > patches;

        //====================================================================//
        //  Domain communication descriptors
        //====================================================================//
        // Local buffer marking
        Array<GhostMarking> local_markings;

        // Neighbour ranks
        Array<uint> neighbour_ranks;

        // Neighbours incoming buffers marking
        Array< Array<GhostMarking> > recv_buffer_markings;

        // __device__
        // neighbours incoming buffers
        Array<Type *> recv_buffers;

        // Neighbours outcoming buffers marking
        Array< Array<GhostMarking> > send_buffer_markings;

        // __device__
        // Neighbours outcoming buffers
        Array<Type *> send_buffers;

        //====================================================================//
        //  Initialization methods
        //====================================================================//
        // Sets initial domain bounds
        void Initialize_domain_bounds();

        // Identifies domain bounds
        void Identify_domain_bounds()

        // Identifies communication descriptors
        void Identify_comm_descriptors();

        // Default
        /// FIXME /// Either remove or switch to c++11: = delete
        Domain() {}

    public:

        // Construction
        /// FIXME /// Data initialization required
        Domain(const Tuple<dim> &, const Tuple<dim> &);

        // Deallocation
        ~Domain();

        //====================================================================//
        //  Get methods
        //====================================================================//
        // Get rank
        inline uint Get_rank() const { return rank; }

        // Get patches index range
        inline uint Get_patch_min_index() const { return domain_bounds[rank]; }
        uint Get_patch_max_index() const { return domain_bounds[rank + 1]; }

        //====================================================================//
        //  Get methods
        //====================================================================//
        // Compute patch count
        inline uint Get_patch_count() const { return patches.Get_size(); }

            /// // copy patch to CPU
            /// void Copy_patch(const uint, const Type *) const;
};

} // namespace maxwell

#include "domain.hpp"

#endif // DOMAIN_H
