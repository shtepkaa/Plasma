// src/.../hilbert_index.cpp

#include "hilbert_index.h"

namespace maxwell {

//==============================================================================
//  Identify descend order in 3D
//==============================================================================
void Identify_descend_order(const uint * arg, uint * order)
{
    const uint comp[THREE_DIM]
        = { arg[1] < arg[2], arg[0] < arg[2], arg[0] < arg[1] };

    order[!comp[0] + !comp[1]] = 2;
    order[comp[0] && comp[1]] = comp[2];
    order[1 + (comp[0] || comp[1])] = !comp[2];
}

//==============================================================================
//  Set bit
//==============================================================================
void Set_bit(uint & num, const uint pos, const uint val)
{
    num = (num & ~(1U << pos)) | (val << pos);
}

//==============================================================================
//  Bitwise cyclic left rotation 
//==============================================================================
// Works only for (shift <= exp)
// -----------------------------------------------------------------------------
uint Rotate_left(const uint val, const uint shift, const uint exp)
{
    return (val << shift | val >> (exp - shift)) & ((1U << exp) - 1);
}

//==============================================================================
//  Bitwise cyclic right rotation 
//==============================================================================
// Works only for (shift <= exp)
// -----------------------------------------------------------------------------
uint Rotate_right(const uint val, const uint shift, const uint exp)
{
    return (val >> shift | val << (exp - shift)) & ((1U << exp) - 1);
}

//==============================================================================
//  Inverse Gray code
//==============================================================================
uint Inverse_Gray_code(const uint code)
{
    uint num = code;

    for (uint c = 1; (1U << c) <= code; ++c) { num ^= code >> c; }

    return num;
}

//==============================================================================
//  Trailing set bit
//==============================================================================
uint Trailing_set_bit(const uint num)
{
    uint count = 0;

    for (uint n = num; n & 1; n >>= 1) { ++count; }

    return count;
}

//==============================================================================
//  Entry
//==============================================================================
uint Entry(const uint ind)
{
    if (ind) { return Gray_code(((ind - 1) >> 1) << 1); }
    else { return 0; }
}

//==============================================================================
//  Direction
//==============================================================================
// (0 <= ind < 2^{exp})
// -----------------------------------------------------------------------------
uint Direction(const uint ind, const uint exp)
{
    if (!ind) { return 0; }
    else if (ind & 1) { return Trailing_set_bit(ind) % exp; }
    else { return Trailing_set_bit(ind - 1) % exp; }
}

//==============================================================================
//  Transform
//==============================================================================
void Transform(const uint entry, const uint dir, uint & bin, const uint exp)
{
    bin = Rotate_right(bin ^ entry, dir + 1, exp);
}

void Inverse_transform(
    const uint entry, const uint dir, uint & bin, const uint exp
)
{
    bin = Rotate_left(bin, dir + 1, exp) ^ entry;
}

//==============================================================================
//  Hilbert index with orientation in 2D
//==============================================================================
uint Hilbert_index_orientation(
    const uint exp, const uint * coords, uint & entry, uint & dir
)
{
    uint index = 0;

    uint loc;
    uint loc_ind;

    for (int i = exp - 1; i >= 0; --i)
    {
        loc = (Get_bit(coords[1], i) << 1) + Get_bit(coords[0], i);
        Transform(entry, dir, loc, 2);

        loc_ind = Inverse_Gray_code(loc);

        entry ^= Rotate_left(Entry(loc_ind), dir + 1, 2);
        dir = (dir + Direction(loc_ind, 2) + 1) & 1;
        index = (index << 2) | loc_ind;
    }

    return index;
}

//==============================================================================
//  Hilbert index in 2D
//==============================================================================
template<>
uint Hilbert_index<TWO_DIM>(const uint exp, const uint * coords)
{
    uint entry = 0;
    uint dir = 0;

    return Hilbert_index_orientation(exp, coords, entry, dir);
}

//==============================================================================
//  Hilbert index in 3D
//==============================================================================
template<>
uint Hilbert_index<THREE_DIM>(
    const uint exp, const uint * coords, const uint entry, const uint dir
)
{
    uint index = 0;

    uint loc;
    uint loc_ind;

    for (int i = exp - 1; i >= 0; --i)
    {
        loc = (Get_bit(coords[2], i) << 2) + (Get_bit(coords[1], i) << 1)
            + Get_bit(coords[0], i);

        Transform(entry, dir, loc, 3);
        loc_ind = Inverse_Gray_code(loc);

        entry ^= Rotate_left(Entry(loc_ind), dir + 1, 3);
        dir = (dir + Direction(loc_ind, 3) + 1) % 3;
        index = (index << 3) | loc_ind;
    }

    return index;
}

//==============================================================================
//  Inverse Hilbert index in 2D
//==============================================================================
template<>
void Inverse_Hilbert_index<TWO_DIM>(
    const uint exp,
    uint * coords,
    const uint index,
    const uint entry,
    const uint dir
)
{
    uint loc;
    uint loc_ind;

    coords[1] = coords[0] = 0;

    for (int i = exp - 1; i >= 0; --i)
    {
        loc_ind = (Get_bit(index, (i << 1) + 1) << 1) + Get_bit(index, i << 1);
        loc = Gray_code(loc_ind);
        Inverse_transform(entry, dir, loc, 2);

        Set_bit(coords[0], i, Get_bit(loc, 0));
        Set_bit(coords[1], i, Get_bit(loc, 1));

        entry ^= Rotate_left(Entry(loc_ind), dir + 1, 2);
        dir = (dir + Direction(loc_ind, 2) + 1) & 1U;
    }
}

//==============================================================================
//  Inverse Hilbert index in 3D
//==============================================================================
template<>
void Inverse_Hilbert_index<THREE_DIM>(
    const uint exp,
    uint * coords,
    const uint index,
    const uint entry,
    const uint dir
)
{
    uint loc;
    uint loc_ind;

    coords[2] = coords[1] = coords[0] = 0;

    for (int i = exp - 1; i >= 0; --i)
    {
        loc_ind
            = (Get_bit(index, 3 * i + 2) << 2)
                + (Get_bit(index, 3 * i + 1) << 1) + Get_bit(index, 3 * i);

        loc = Gray_code(loc_ind);
        Inverse_transform(entry, dir, loc, 3);

        Set_bit(coords[0], i, Get_bit(loc, 0));
        Set_bit(coords[1], i, Get_bit(loc, 1));
        Set_bit(coords[2], i, Get_bit(loc, 2));

        entry ^= Rotate_left(Entry(loc_ind), dir + 1, 3);
        dir = (dir + Direction(loc_ind, 3) + 1) % 3;
    }
}

//==============================================================================
//  General Hilbert index with orientation in 2D
//==============================================================================
uint General_Hilbert_index_orientation(
    const uint * exps, const uint * coords, uint & entry, uint & dir
)
{
    if (coords[0] >= (1U << exps[0]) || coords[1] >= (1U << exps[1]))
    {
        return ~0;
    }

    dir = (exps[0] < exps[1]);

    uint index = 0;
    uint min_exp = exps[1 - dir];

    uint loc;

    for (int i = exps[dir] - 1; i >= int(min_exp); --i)
    {
        loc = Get_bit(coords[dir], i);
        index += loc * (1 << (i + min_exp));
        coords[dir] -= loc * (1 << i);
    }

    // Calculate entry and direction
    if (min_exp)
    {
        index += Hilbert_index_orientation(min_exp, coords, entry, dir);
    }

    return index;
}

//==============================================================================
//  General Hilbert index in 1D 
//==============================================================================
uint General_Hilbert_index<ONE_DIM>(const uint * exps, uint * coords)
{
    return (coords[0] >= (1U << exps[0])? ~0: coords[0];
}

//==============================================================================
//  General Hilbert index in 2D 
//==============================================================================
uint General_Hilbert_index<TWO_DIM>(const uint * exps, uint * coords)
{
    uint entry = 0;
    uint dir = 0;

    return General_Hilbert_index_orientation(exps, coords, entry, dir);
}

//==============================================================================
//  General Hilbert index 3D
//==============================================================================
template<>
uint General_Hilbert_index<THREE_DIM>(const uint * exps, const uint * coords)
{
    if (
        coords[0] >= (1U << exps[0]) || coords[1] >= (1U << exps[1])
        || coords[2] >= (1U << exps[2])
    )
    {
        return ~0;
    }

    uint order[THREE_DIM];
                      
    // compare exponents
    Identify_descend_order(exps, order);

    uint ord_coords[THREE_DIM]
        = { coords[order[0]], coords[order[1]], coords[order[2]] };

    const uint min_exp = exps[order[2]];

    // Approach on flattened 2D grid along max and med exps
    // 3D grid is projected along min_exp axis
    // Erase last min_exp bits, not relevant for this phase
    const uint flat_exps[TWO_DIM]
        = { exps[order[0]] - min_exp, exps[order[1]] - min_exp };

    const uint flat_coords[TWO_DIM]
        = { ord_coords[0] >> min_exp, ord_coords[1] >> min_exp };

    uint entry = 0;
    uint dir = 0;

    uint index
        = General_Hilbert_index_orientation(flat_exps, flat_coords, entry, dir)
            * (1U << (3 * min_exp));

    // In the local cube of side 2^{min_exp} the local entry point and initial
    // direction of the local Hilbert curve has already been determined
    // The relative position in the local cube is given by last min_exp bits
    // of the position
    for (int d = 0; d < 3; ++d) { ord_coords[d] &= (1U << min_exp) - 1; }

    return index + Hilbert_index<THREE_DIM>(min_exp, ord_coords, entry, dir);
}

//==============================================================================
//  Inverse general Hilbert index in 1D
//==============================================================================
template<>
void Inverse_general_Hilbert_index<ONE_DIM>(
    const uint *, uint * coords, const uint index
)
{
    coords[0] = index;
}

//==============================================================================
//  Inverse general Hilbert index in 2D
//==============================================================================
template<>
void Inverse_general_Hilbert_index<TWO_DIM>(
    const uint * exps, uint * coords, const uint index
)
{
    // Compare the exponents and target a coordinate which must be shifted
    uint dir = (exps[0] < exps[1]);
    uint min_exp = exps[1 - dir];
    uint ind = index;

    uint shift = 0;
    uint loc;

    // Define in which sub-hypercube of side 2^{min_exp} the point is
    for (int i = exps[0] + exps[1] - 1; i >= int(2 * min_exp); --i)
    {
        loc = Get_bit(ind, i);
        shift += loc * (1 << (i - min_exp));
        ind -= loc * (1 << i);
    }

    // Run the 2D inversion algorithm in the sub-hypercube
    Inverse_Hilbert_index<TWO_DIM>(min_exp, coords, ind, 0, dir);

    // Shift the appropriate coordinate by the necessary value
    coords[dir] += shift;
}

//==============================================================================
//  Inverse general Hilbert index in 3D
//==============================================================================
template<>
void Inverse_general_Hilbert_index<THREE_DIM>(
    const uint * exps, uint * coords, const uint index
)
{
    uint order[THREE_DIM];
                      
    // Compare exponents
    Identify_descend_order(exps, order);

    const uint min_exp = exps[order[2]];

    const uint flat_exps[TWO_DIM]
        = { exps[order[0]] - min_exp, exps[order[1]] - min_exp };

    uint loc_coords[THREE_DIM];

    // Localize in which sub-hypercube the point is
    // Do not account for the first (3 * min_exp) bits of index
    uint loc_ind = index >> (3 * min_exp);

    // Run the 2D inversion algorithm on reduced domain
    Inverse_general_Hilbert_index<TWO_DIM>(flat_exps, loc_coords, loc_ind);

    // Run the 3D inversion algorithm on this cube with the correct entry point
    // and direction
    uint entry = 0;
    uint dir = 0;

    // Run the 2D indexgenerator in order to evaluate the entry and direction
    loc_ind
        = General_Hilbert_index_orientation(flat_exps, loc_coords, entry, dir);

    // Store transformed coordinates in the resulting global frame
    coords[order[0]] = loc_coords[0] * (1U << min_exp);
    coords[order[1]] = loc_coords[1] * (1U << min_exp);
    coords[order[2]] = 0;

    // Use only first bits of index for the local hypercube
    loc_ind = index & ((1U << (3 * min_exp)) - 1);

    // Run the cubic inversion algorithm in the local sub hypercube
    Inverse_Hilbert_index<THREE_DIM>(min_exp, loc_coords, loc_ind, entry, dir);

    // Store the resulting coordinates
    coords[order[0]] += loc_coords[0];
    coords[order[1]] += loc_coords[1];
    coords[order[2]] += loc_coords[2];
}

} // namespace maxwell

// src/.../hilbert_index.cpp
