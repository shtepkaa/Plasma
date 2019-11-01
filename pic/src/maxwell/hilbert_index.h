#ifndef HILBERT_INDEX_H
#define HILBERT_INDEX_H

namespace maxwell {

//==============================================================================
//  Identify descend order
//==============================================================================
// Finds the descending order of argument entries
// -----------------------------------------------------------------------------
static void Identify_descend_order(const uint *, uint *);

//==============================================================================
//  Bit access / mutate
//==============================================================================
static uint Get_bit(const uint num, const uint pos) { return (num >> pos) & 1; }
static void Set_bit(uint &, const uint, const uint);

//==============================================================================
//  Bitwise cyclic left rotation 
//==============================================================================
// After rotation retains only the given number of the first bits
// -----------------------------------------------------------------------------
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
// Computes number of trailing set bits in the binary representation of number,
// which is also the inter sub-hypercube direction
// -----------------------------------------------------------------------------
static uint Trailing_set_bit(const uint);

//==============================================================================
//  Direction and entry for intra sub-hypercube
//==============================================================================
// Computes intra sub-hypercube entry and direction for a given index
// -----------------------------------------------------------------------------
static uint Entry(const uint);
static uint Direction(const uint, const uint);

//==============================================================================
//  Standard transformation of a Gray code with the given entry and direction
//==============================================================================
// Computes the transformation such that the Gray code ordering
// of sub-hypercubes in the Hilbert curve defined by the entry and direction
// will map to the standard binary reflected Gray code
// -----------------------------------------------------------------------------
static void Transform(const uint, const uint, uint &, const uint);
static void Inverse_transform(const uint, const uint, uint &, const uint);

//==============================================================================
//  Hilbert index
//==============================================================================
// Calculates the Hilbert index and the orientation of a patch with given
// coordinates for a simulation box with 2^{exp} patches per side
// (i. e. 2^{2 * exp} patches in total)
// -----------------------------------------------------------------------------
static uint Hilbert_index_orientation(const uint, const uint *, uint &, uint &);

// -----------------------------------------------------------------------------
// Calculates the Hilbert index of a patch with given coordinates
// for a simulation box with 2^{exp} patches per side (i. e. 2^{dim * exp}
// patches in total)
// -----------------------------------------------------------------------------
template<Dimension>
static uint Hilbert_index(const uint, const uint *, const uint, const uint);

// -----------------------------------------------------------------------------
// Calculates the coordinates of a patch with a given Hilbert index
// in a simulation box with 2^{exp} patches per side (i. e. 2^{dim * exp}
// patches in total)
// -----------------------------------------------------------------------------
template<Dimension>
static void Inverse_Hilbert_index(
    const uint, uint *, const uint, const uint, const uint
);

//==============================================================================
//  General Hilbert index
//==============================================================================
// Calculates the Hilbert index and orientation of a patch with given
// coordinates for a simulation box with 2^{exps[d]} patches per side
// (i. e. 2^{exps[0] + exps[1]} patches in total)
// -----------------------------------------------------------------------------
static uint General_Hilbert_index_orientation(
    const uint *, const uint *, uint &, uint &
);

// -----------------------------------------------------------------------------
// Calculates the compact Hilbert index of a patch with given coordinates
// for a simulation box with 2^{exps[d]} patches per side
// -----------------------------------------------------------------------------
template<Dimension>
uint General_Hilbert_index(const uint *, const uint *);

// -----------------------------------------------------------------------------
// Calculates the coordinates of a patch for a given Hilbert index
// in a simulation box with 2^{exps[d]} patches per side
// -----------------------------------------------------------------------------
template<Dimension>
void Inverse_general_Hilbert_index(const uint *, uint *, const uint);

} // namespace maxwell

#endif // HILBERT_INDEX_H
