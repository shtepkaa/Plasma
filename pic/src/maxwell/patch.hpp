#ifndef PATCH_HPP
#define PATCH_HPP

/// FIXME /// #include "types.h"
/// FIXME /// #include "utility.h"

namespace maxwell {

//============================================================================//
//  Construct markings
//============================================================================//
template<Dim dim>
void Patch::Initialize_markings()
{
    Tuple<dim> shifts;

    for (uint s = 0; s < ghost_markings.Get_size(); ++s)
    {
        uint index = s;

        for (uint d = 0; d < dim - 1; ++d)
        {
            shifts[d] += (index % 3) - 1;
            index /= 3;
        }

        shifts[dim - 1] += index - 1;

        // now all shifts are in { ~0, 0, 1 }

        ghost_markings[s].index
            = General_Hilbert_index<dim>(grid_exps, grid_coords + shifts);

        for (uint d = 0; d < dim - 1; ++d)
        {
            ghost_markings[s].sizes[d]
                = (shifts[s])? grid_sizes[shifts[s] + 1]: ghost_width;
        }
    }
}

//============================================================================//
//  Constructor
//============================================================================//
template<Dim dim, typename Type>
Patch::Patch(
    const Tuple<dim> & exps,
    const Tuple<dim> & coords,
    const Tuple<dim> & sizes,
    const uint width
):
    grid_exps(exps),
    grid_coords(coords),
    grid_sizes(sizes),
    ghost_width(width),
    actual_sizes(sizes += 2 * width),
    markings(Total_hypercube_count(dim)),
    data_size(Product<dim>(actual_sizes)),
    data(NULL)
{
    Construct_markings();
    cudaMalloc(&data, data_size * sizeof(Type));
}

} // namespace maxwell

#endif // PATCH_HPP
