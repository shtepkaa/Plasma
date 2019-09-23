#ifndef DOMAIN_H
#define DOMAIN_H

/// #include <vector>

namespace maxwell {

struct GhostId
{
    uint send_ind;
    uint recv_ind;

    uint size;
    uint start;
};

//============================================================================//
//  Domain
//============================================================================//
template<Dim, Order>
class Domain
{
    private:

        // rank
        uint rank;

        // global range of ranks
        uint range;

        // neighbour count of ranks
        uint neighbour_count;

        // neighbour ranks
        uint * neighbour_ranks;

            /// // Cartesian grid coordinates
            /// uint * grid_coords;

        // patches Hilbert index range
        uint patch_min_index;
        uint patch_max_index;

        // patches count
        uint patch_count;

        // patches arranged in specified Order
        Patch ** patches;

        // patch Cartesian grid sizes
        uint * patch_sizes;

        // border patch count
        uint border_patch_count;

        // border patch indices
        uint * border_patch_indices;

        // neighbours incoming buffers markup
        GhostId ** recv_buffer_markups;

        // neighbours incoming buffers
        double * recv_buffers;

        // neighbours outcoming buffers markup
        GhostId ** send_buffer_markups;

        // neighbours outcoming buffers
        double * send_buffers;

        // identify neighbour ranks
        void Set_neighbour_ranks();

        // identify border patch indices
        void Set_border_patch_indices();

        // default initializer
        Domain() {}

    public:

        // initialize
        Domain(const uint, const uint, const uint *);

        // deallocate
        ~Domain();

        // get rank
        inline uint Get_rank() const { return rank; }

        // get Hilbert index range
        inline uint Get_patch_min_index() const { return patch_min_index; }
        inline uint Get_patch_max_index() const { return patch_max_index; }

        // copy patch to CPU
        void Copy_patch(const uint, double *) const;
};

} // namespace maxwell

#include "domain.hpp"

#endif // DOMAIN_H
