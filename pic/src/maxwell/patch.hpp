#ifndef PATCH_HPP
#define PATCH_HPP

/// FIXME /// #include "types.h"

namespace maxwell {

//============================================================================//
template<Dim>
void Patch::Get_index()
{
    /// TODO ///
}

//============================================================================//
template<Dim>
void Patch::Get_neighbours()
{
    /// TODO ///
}

//============================================================================//
template<Dim>
Patch::Patch()
{
    /// IS REQUIRED? ///
}

//============================================================================//
template<Dim dim>
Patch::Patch(const uint * coords, const uint * sizes)
{
    Init_array(dim, grid_coords, coords);
    Init_array(dim, grid_sizes, sizes);
}

//============================================================================//
template<Dim>
Patch::~Patch()
{
    /// TODO ///
}

//============================================================================//
template<Dim>
void Patch::Set_rank(const uint r) { rank = r; }

//============================================================================//
template<Dim dim>
void Patch::Get_ghost(const uint pos)
{
    /// TODO ///
}

} // namespace maxwell

#endif // PATCH_HPP
