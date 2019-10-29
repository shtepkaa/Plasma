#ifndef HILBERT_INDEX_H
#define HILBERT_INDEX_H

namespace maxwell {

//==============================================================================
//  Identify descend order
//==============================================================================
static void Identify_descend_order(const uint *, uint *);

//==============================================================================
//  Bit access
//==============================================================================
static uint Get_bit(const uint num, const uint pos) { return (num >> pos) & 1; }
static void Set_bit(uint &, const uint, const uint);

//==============================================================================
//  Bitwise cyclic left rotation 
//==============================================================================
static uint Rotate_left(const uint, const uint, const uint);
static uint Rotate_right(const uint, const uint, const uint);

//==============================================================================
//  Binary reflected Gray code
//==============================================================================
static uint Gray_code(const uint num) { return num ^ (num >> 1); }
static uint Inverse_Gray_code(const uint);

//==============================================================================
//  Trailing set bit
//==============================================================================
static uint Trailing_set_bit(const uint);

//==============================================================================
//  Direction and entry for intra sub-hypercube
//==============================================================================
static uint Direction(const uint, const uint);
static uint Entry(const uint);

//==============================================================================
//  Standard transformation of Gray code with given entry and direction
//==============================================================================
static void Transform(const uint, const uint, uint &, const uint);
static void Inverse_transform(const uint, const uint, uint &, const uint);

//==============================================================================
//  Hilbert index
//==============================================================================
static uint Hilbert_index_orientation(const uint, const uint *, uint &, uint &);

template<Dimension>
static uint Hilbert_index(const uint, const uint *, const uint, const uint);

template<Dimension>
static void Inverse_Hilbert_index(
    const uint, uint *, const uint, const uint, const uint
);

//==============================================================================
//  General Hilbert index
//==============================================================================
static uint General_Hilbert_index_orientation(
    const uint *, const uint *, uint &, uint &
);

template<Dimension>
uint General_Hilbert_index(const uint *, const uint *);

template<Dimension>
void Inverse_general_Hilbert_index(const uint *, uint *, const uint);

} // namespace maxwell

#endif // HILBERT_INDEX_H
