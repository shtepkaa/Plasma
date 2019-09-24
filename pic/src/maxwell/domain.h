#ifndef DOMAIN_H
#define DOMAIN_H

/// #include <vector>

namespace maxwell {

//============================================================================//
//  GhostIdentifier
//============================================================================//
// markup structure
struct GhostIdentifier
{
    uint send_ind;
    uint recv_ind;

    uint size;
    uint start;
};

//============================================================================//
//  Domain
//============================================================================//
template<Dim, Order, ArithmeticType>
class Domain
{
    private:

        // rank
        uint rank;

        // global range of ranks
        uint range;

        // neighbour ranks
        Array<uint> neighbour_ranks;

            /// // Cartesian grid coordinates
            /// uint * grid_coords;

        // patch Cartesian grid sizes
        uint * patch_sizes;

        // patches Hilbert index range
        uint patch_min_index;
        uint patch_max_index;

        // border patch indices
        Array<uint> border_patch_indices;

        // patches arranged in specified Order
        Array<Patch> patches;

        // neighbours incoming buffers markup
        Array< Array<GhostIdentifier> > recv_buffer_markups;

        // neighbours incoming buffers
        Array<ArithmeticType> recv_buffers;

        // neighbours outcoming buffers markup
        Array< Array<GhostIdentifier> > send_buffer_markups;

        // neighbours outcoming buffers
        Array<ArithmeticType> send_buffers;

        // identify neighbour ranks
        void Set_neighbour_ranks();

        // identify border patch indices
        void Set_border_patch_indices();

        /// FIXME /// Either remove or switch to c++11: = delete
        // default initializer
        Domain() {}

    public:

        // initializion
        Domain(const uint, const uint, const uint *);

        // deallocation
        ~Domain();

        // get rank
        inline uint Get_rank() const { return rank; }

        // get Hilbert index range
        inline uint Get_patch_min_index() const { return patch_min_index; }
        inline uint Get_patch_max_index() const { return patch_max_index; }

        // copy patch to CPU
        void Copy_patch(const uint, ArithmeticType *) const;
};

} // namespace maxwell

#include "domain.hpp"

#endif // DOMAIN_H
