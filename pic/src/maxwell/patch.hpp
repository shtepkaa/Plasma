#ifndef PATCH_HPP
#define PATCH_HPP

/// FIXME /// #include "types.h"

namespace maxwell {

//============================================================================//
//  Identify_indices
//============================================================================//
template<Dim>
void Patch::Identify_indices()
{
    /// TODO ///
}

//============================================================================//
//  Initialization
//============================================================================//
template<Dim dim>
Patch::Patch(
    const uint rk,
    const uint * dim_exps,
    const uint * coords,
    const uint * sizes,
    const uint width
):
    rank(rk),
    /// index(Get_hilbert_index<dim>(dim_exps, coords)),
    indices(Total_hypercube_count(dim)),
    ghost_width(width),
    data_size(1),
    data(NULL)
{

    for (uint i = 0; i < indices.Get_size(); ++i)
    {
        Get_hilbert_index<dim>(dim_exps, coords)
    }

    for (uint d = 0; d < dim; ++d) { data_size *= sizes[d] + 2 * width; }

    copy(dim, grid_coords, 1, coords, 1);
    copy(dim, grid_sizes, 1, sizes, 1);
}

//============================================================================//
//  Deallocation
//============================================================================//
template<Dim>
Patch::~Patch()
{
    /// TODO ///
}

//============================================================================//
//  Get_ghost
//============================================================================//
template<Dim dim>
void Patch::Get_ghost(const uint pos)
{
    /// TODO ///
}

} // namespace maxwell

#endif // PATCH_HPP
