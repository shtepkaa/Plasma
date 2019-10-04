#ifndef PATCH_HPP
#define PATCH_HPP

/// FIXME /// #include "types.h"

namespace maxwell {

//============================================================================//
//  Identify_indices
//============================================================================//
template<Dim>
void Patch::Identify_ids()
{
    /// TODO ///
}

//============================================================================//
//  Initialization
//============================================================================//
template<Dim dim>
Patch::Patch(
    const uint * exps,
    const uint * coords,
    const uint * sizes,
    const uint width
):
    ghost_width(width),
    data(NULL),
    markings(Total_hypercube_count(dim))
{
    /// FIXME /// -> use Tuple in initializer list /// copy(dim, exps, 1, grid_exps, 1);
    /// FIXME /// -> use Tuple in initializer list /// copy(dim, coords, 1, grid_coords, 1);
    /// FIXME /// -> use Tuple in initializer list /// copy(dim, sizes, 1, grid_sizes, 1);
    /// FIXME /// -> use Tuple in initializer list /// axpy(dim, 2 * width, sizes, 1, actual_sizes, 1);

    for (uint i = 0; i < indices.Get_size(); ++i)
    {
        markings.index = Get_hilbert_index<dim>(grid_exps, grid_coords);
    }
}

} // namespace maxwell

#endif // PATCH_HPP
