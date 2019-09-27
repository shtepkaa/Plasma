#ifndef PATCH_HPP
#define PATCH_HPP

/// FIXME /// #include "types.h"


template<int N>
struct Factorial {
  enum { value = N * Factorial<N - 1>::value };
};

template<>
struct Factorial<0> {
  enum { value = 1 };
};

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
Patch::Patch(const uint rank, const uint * coords, const uint * sizes):
    neighbour_indices(), grid_coords(dim, coords), grid_sizes(dim, sizes)
{}

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
