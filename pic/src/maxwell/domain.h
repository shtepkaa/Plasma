#ifndef DOMAIN_H
#define DOMAIN_H

/// FIXME /// #include "patch.h"
/// FIXME /// #include "utility.h"

namespace maxwell {

////////////////////////////////////////////////////////////////////////////////
//
//  BufferMarking
//
////////////////////////////////////////////////////////////////////////////////
// Implements buffer marking 
template<Dim dim>
struct BufferMarking
{
    uint send_index;
    uint recv_index;

    Tuple<dim> sizes;

    uint offset;
};

////////////////////////////////////////////////////////////////////////////////
//
//  Domain
//
////////////////////////////////////////////////////////////////////////////////
// Implements data structure per process corresponding to a single MPI rank
// Template parameter: dim -- dimensionality of the grid
// Template parameter: ord -- order of patch numeration
// Template parameter: Type -- supported arithmetic type
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

        // Neighbour ranks
        Array<uint> neighbour_ranks;

        // Patch Cartesian grid sizes
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
        Array<GhostMarking> local_buffer_markings;

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
        // Identify neighbour ranks
        void Identify_neighbour_ranks();

        /// FIXME /// Either remove or switch to c++11: = delete
        // default initializer
        Domain() {}

    public:

        // Construction
        Domain(const uint, const uint, const Tuple<dim> &);

        // Deallocation
        ~Domain();

        //====================================================================//
        //  Get methods
        //====================================================================//
        // Get rank
        inline uint Get_rank() const { return rank; }

        // Get Hilbert index range
        uint Get_patch_min_index() const { return domain_bounds[rank]; }
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
