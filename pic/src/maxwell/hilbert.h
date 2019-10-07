#ifndef HILBERT_H
#define HILBERT_H

namespace maxwell {

//============================================================================//
//  Identify descend order
//============================================================================//
static void Identify_descend_order(const Tuple<3D> &, Tuple<3D> &);

//============================================================================//
//  Bit access
//============================================================================//
static uint Get_bit(const uint n, const uint pos) { return (n >> pos) & 1U; }
static void Set_bit(uint &, const uint, const uint);

//============================================================================//
//  Bitwise cyclic left rotation 
//============================================================================//
static uint Rotate_left(const uint, const uint, const uint);
static uint Rotate_right(const uint, const uint, const uint);

//============================================================================//
//  Binary reflected Gray code
//============================================================================//
static uint Gray_code(const uint num) { return num ^ (num >> 1); }
static uint Gray_code_inverse(const uint);

//============================================================================//
//  Trailing set bit
//============================================================================//
static uint Trailing_set_bit(const uint);

//============================================================================//
//  Direction and entry for intra sub-hypercube
//============================================================================//
static uint Direction(const uint, const uint);
static uint Entry(const uint);

//============================================================================//
//  Standard transformation of Gray code with given entry and direction
//============================================================================//
static void Transform(const uint, const uint, uint &, const uint);
static void Transform_inverse(const uint, const uint, uint &, const uint);

//============================================================================//
//  Hilbert index
//============================================================================//
static uint Hilbert_index_orientation(
    const uint, const Tuple<2D> &, uint &, uint &
);

template<Dim dim>
uint Hilbert_index(const uint, const Tuple<dim> &, const uint, const uint);

template<Dim dim>
void Hilbert_index_inverse(
    const uint, Tuple<dim> &, const uint, const uint, const uint
);

//============================================================================//
//  General Hilbert index
//============================================================================//
static uint General_Hilbert_index_orientation(
    const Tuple<2D> &, const Tuple<dim> &, uint &, uint &
);

template<Dim dim>
uint General_Hilbert_index(const Tuple<dim> &, const Tuple<dim> &);

template<Dim dim>
void General_Hilbert_index_inverse(
    const Tuple<dim> &, Tuple<dim> &, const uint
);

} // namespace maxwell

#include "hilbert.hpp"

#endif // HILBERT_H
