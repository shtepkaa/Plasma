#ifndef HILBERT_H
#define HILBERT_H

namespace maxwell {

//============================================================================//
//  Identify descend order
//============================================================================//
static void Identify_descend_order(const uint *, uint *);

//============================================================================//
//  Bit access
//============================================================================//
static uint Get_bit(const uint i, const uint k) { return (i >> k) & 1U; }
static void Set_bit(uint &, const uint, const uint);

//============================================================================//
//  Bitwise cyclic left rotation 
//============================================================================//
static uint Rotate_left(const uint, const uint, const uint);
static uint Rotate_right(const uint, const uint, const uint);

//============================================================================//
//  Binary reflected Gray code
//============================================================================//
static uint Gray_code(const uint i) { return i ^ (i >> 1); }
static uint Gray_code_inverse(const uint);

//============================================================================//
//  Trailing set bit
//============================================================================//
static uint Trailing_set_bit(const uint);

//============================================================================//
//  Direction and entry for intra su-hypercube
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
static uint Hilbert_index_orientation(const uint, const uint *, uint &, uint &);

template<Dim>
uint Hilbert_index(const uint, const uint *, const uint, const uint);

template<Dim>
void Hilbert_index_inverse(
    const uint, uint *, const uint, const uint, const uint
);

//============================================================================//
//  General Hilbert index
//============================================================================//
static uint General_hilbert_index_orientation(
    const uint *, const int *, uint &, uint &
);

template<Dim>
uint General_hilbert_index(const uint *, const int *);

template<Dim>
void General_hilbert_index_inverse(const uint *, uint *, const uint);

} // namespace maxwell

#include "hilbert.hpp"

#endif // HILBERT_H
