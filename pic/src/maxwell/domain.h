#ifndef DOMAIN_H
#define DOMAIN_H

/// #include ""

/// #include <>

namespace maxwell {

//============================================================================//
//  Domain
//============================================================================//
// Implements data structure per process
// corresponds to single MPI rank
/// >>> /// implementation should probably be changed to singleton class
template<Dim dim, /* ??? Order ord,*/ ArithmeticType Type>
class Domain
{
    private:

        // MPI rank
        uint rank;

        // global range of ranks
        uint range;

        // grid Hilbert index decomposition
        Array<uint> domain_bounds;

        // neighbour ranks
        Array<uint> neighbour_ranks;

            // patch Cartesian grid sizes
            uint patch_sizes[dim];

        // patches arranged in specified Order
        Array<Patch> patches;

        // local buffer marking
        Array<GhostMarking> local_buffer_markings;

        // neighbours incoming buffers marking
        Array< Array<GhostMarking> > recv_buffer_markings;

        // should probably (re)allocate independantly
        // for different buffers
            /// // neighbours incoming buffers
            /// Array<uint> recv_buffers;
        // __device__
        // neighbours incoming buffers
        Array<Type *> recv_buffers;

        // neighbours outcoming buffers marking
        Array< Array<GhostMarking> > send_buffer_markings;

        // should probably (re)allocate independantly
        // for different buffers
            /// // neighbours outcoming buffers
            /// Array<uint> send_buffers;
        // __device__
        // neighbours outcoming buffers
        Array<Type *> send_buffers;

        // identify neighbour ranks
        void Identify_neighbour_ranks();

        /// FIXME /// Either remove or switch to c++11: = delete
        // default initializer
        Domain() {}

    public:

        // initialization
        Domain(const uint, const uint, const uint *);

        // deallocation
        ~Domain();

        // get rank
        inline uint Get_rank() const { return rank; }

        // get Hilbert index range
        uint Get_patch_min_index() const { return domain_bounds[rank]; }
        uint Get_patch_max_index() const { return domain_bounds[rank + 1]; }

        // copy patch to CPU
        void Copy_patch(const uint, const Type *) const;
};

} // namespace maxwell

#include "domain.hpp"

#endif // DOMAIN_H
