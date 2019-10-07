// hilbert.cpp

namespace maxwell {

//============================================================================//
//  Identify descend order
//============================================================================//
void Identify_descend_order(const uint * arg, uint * order)
{
    const uint comp[3] = { arg[1] < arg[2], arg[0] < arg[2], arg[0] < arg[1] };

    order[!comp[0] + !comp[1]] = 2;
    order[comp[0] && comp[1]] = comp[2];
    order[1 + (comp[0] || comp[1])] = !comp[2];
}

//============================================================================//
//  Set bit
//============================================================================//
// Sets the bit of a number with a given position to a value
void Set_bit(uint & num, const uint pos, const uint val)
{
    num = (num & ~(1 << pos)) | (val << pos);
}

//============================================================================//
//  Bitwise cyclic left rotation 
//============================================================================//
// Retains only the first exp bits
// works only for (shift <= exp)
uint Rotate_left(const uint val, const uint shift, const uint exp)
{
    return (val << shift | val >> (exp - shift)) & ((1 << exp) - 1);
}

//============================================================================//
//  Bitwise cyclic right rotation 
//============================================================================//
// Retains only the first exp bits
// works only for (shift <= exp)
uint Rotate_right(const uint val, const uint shift, const uint exp)
{
    return (val >> shift | val << (exp - shift)) & ((1 << exp) - 1);
}

//============================================================================//
//  Gray code inverse
//============================================================================//
// Calculates the non-negative integer i such that Gray_code(i) = code,
// for non-negative integer g
uint Gray_code_inverse(const uint code)
{
    uint i = code;

    for (uint c = 1; (1U << c) <= code; ++c) { i ^= code >> c; }

    return i;
}

//============================================================================//
//  Trailing set bit
//============================================================================//
// Computes number of trailing set bits in the binary representation of number,
// which is also the inter sub-hypercube direction
uint Trailing_set_bit(const uint num)
{
    uint count = 0;

    for (uint n = num; n & 1; n >>= 1) { ++count; }

    return count;
}

//============================================================================//
//  Direction
//============================================================================//
// Computes intra sub-hypercube direction for (0 <= ind < 2^{exp})
uint Direction(const uint ind, const uint exp)
{
    if (!ind) { return 0; }
    else if (ind & 1) { return Trailing_set_bit(ind) % exp; }
    else { return Trailing_set_bit(ind - 1) % exp; }
}

//============================================================================//
//  Entry
//============================================================================//
// Computes entry point for (0 <= ind < 2^{exp})
uint Entry(const uint ind)
{
    if (ind) { return Gray_code(((ind - 1) >> 1) << 1); }
    else { return 0; }
}

//============================================================================//
//  Transform
//============================================================================//
// Computes transformation such that the Gray code ordering of sub-hypercubes
// in the Hilbert curve defined by "entry" and "dir"
// will map to the standard binary reflected Gray code
void Transform(const uint entry, const uint dir, uint & bin, const uint exp)
{
    bin = Rotate_right(bin ^ entry, dir + 1, exp);
}

void Transform_inverse(
    const uint entry, const uint dir, uint & bin, const uint exp
)
{
    bin = Rotate_left(bin, dir + 1, exp) ^ entry;
}

//============================================================================//
//  Hilbert index 2D with orientation
//============================================================================//
// Calculates the Hilbert index and orientation of a patch
// with given coordinates for simulation box with 2^{exp} patches per side
// (2^(2 * exp) patches in total)
uint Hilbert_index_orientation(
    const uint exp, const uint * coordinates, uint & entry, uint & dir
)
{
    uint index = 0;

    uint loc;
    uint loc_ind;

    for (uint i = exp - 1; i >= 0; --i)
    {
        loc = (Get_bit(coordinates[1], i) << 1) + Get_bit(coordinates[0], i);
        Transform(entry, dir, loc, 2);

        loc_ind = Gray_code_inverse(loc);

        entry ^= Rotate_left(Entry(loc_ind), dir + 1, 2);
        dir = (dir + Direction(loc_ind, 2) + 1) & 1;
        index = (index << 2) | loc_ind;
    }

    return index;
}

//============================================================================//
//  Hilbert index 2D
//============================================================================//
// Calculates the Hilbert index of patch of given coordinates
// for a simulation box with 2^{exp} patches
// per side (2^(2 * exp) patches in total)
template<>
uint Hilbert_index<2D>(const uint exp, const uint * coordinates)
{
    uint entry = 0;
    uint dir = 0;

    return Hilbert_index_orientation(exp, coordinates, entry, dir);
}

//============================================================================//
//  Hilbert index 3D
//============================================================================//
// Calculates the Hilbert index of patch of given coordinates
// for a simulation box with 2^{exp} patches
// per side (2^(3 * exp) patches in total)
template<>
uint Hilbert_index<3D>(
    const uint exp, const uint * coordinates, const uint entry, const uint dir
)
{
    uint index = 0;

    uint loc;
    uint loc_ind;

    for (uint i = exp - 1; i >= 0; --i)
    {
        loc = (Get_bit(coordinates[2], i) << 2)
            + (Get_bit(coordinates[1], i) << 1) + Get_bit(coordinates[0], i);

        Transform(entry, dir, loc, 3);
        loc_ind = Gray_code_inverse(loc);

        entry ^= Rotate_left(Entry(loc_ind), dir + 1, 3);
        dir = (dir + Direction(loc_ind, 3) + 1) % 3;
        index = (index << 3) | loc_ind;
    }

    return index;
}

//============================================================================//
//  Hilbert index 2D inverse
//============================================================================//
// Calculates the coordinates of a patch with given Hilbert index
// in a simulation box with 2^{exp} patches per side
// (2^(2 * exp) patches in total)
template<>
void Hilbert_index_inverse<2D>(
    const uint exp,
    uint * coordinates,
    const uint index,
    const uint entry,
    const uint dir
)
{
    uint loc;
    uint loc_ind;

    coordinates[0] = coordinates[1] = 0;

    for (uint i = exp - 1; i >= 0; --i)
    {
        loc_ind = (Get_bit(index, (i << 1) + 1) << 1) + Get_bit(index, i << 1);
        loc = Gray_code(loc_ind);
        Transform_inverse(entry, dir, loc, 2);

        Set_bit(coordinates[0], i, Get_bit(loc, 0));
        Set_bit(coordinates[1], i, Get_bit(loc, 1));

        entry ^= Rotate_left(Entry(loc_ind), dir + 1, 2);
        dir = (dir + Direction(loc_ind, 2) + 1) & 1U;
    }
}

//============================================================================//
//  Hilbert index 3D inverse
//============================================================================//
// Calculates the coordinates of a patch with given Hilbert index
// in a simulation box with 2^{exp} patches per side
// (2^(3 * exp) patches in total)
template<>
void Hilbert_index_inverse<3D>(
    const uint exp,
    uint * coordinates,
    const uint index,
    const uint entry,
    const uint dir
)
{
    uint loc;
    uint loc_ind;

    coordinates[0] = coordinates[1] = coordinates[2] = 0;

    for (uint i = exp - 1; i >= 0; --i)
    {
        loc_ind
            = (Get_bit(index, 3 * i + 2) << 2)
                + (Get_bit(index, 3 * i + 1) << 1) + Get_bit(index, 3 * i);

        loc = Gray_code(loc_ind);
        Transform_inverse(entry, dir, loc, 3);

        Set_bit(coordinates[0], i, Get_bit(loc, 0));
        Set_bit(coordinates[1], i, Get_bit(loc, 1));
        Set_bit(coordinates[2], i, Get_bit(loc, 2));

        entry ^= Rotate_left(Entry(loc_ind), dir + 1, 3);
        dir = (dir + Direction(loc_ind, 3) + 1) % 3;
    }
}

//============================================================================//
//  General Hilbert index with orientation
//============================================================================//
// Calculates Hilbert index and orientation of a patch with given coordinates
// for a simulation box with 2^{exps[d]} patches per side
// (2^{exps[0] + exps[1]} patches in total)
uint General_Hilbert_index_orientation(
    const uint * exps, const uint * coordinates, uint & entry, uint & dir
)
{
    if (coordinates[0] >= (1 << exps[0]) || coordinates[1] >= (1 << exps[1]))
    {
        return -1;
    }

    uint coords[2] = { coordinates[0], coordinates[1] }; 

    dir = (exps[0] < exps[1]);

    uint index = 0;
    uint min_exp = exps[1 - dir];

    uint loc;

    for (uint i = exps[dir] - 1; i >= min_exp; --i)
    {
        loc = Get_bit(coords[dir], i);
        index += loc * (1 << (i + min_exp));
        coords[dir] -= loc * (1 << i);
    }

    // calculate entry and direction
    if (min_exp)
    {
        index += Hilbert_index_orientation(min_exp, coords, entry, dir);
    }

    return index;
}

//============================================================================//
//  General Hilbert 2D index
//============================================================================//
// Calculates the compact Hilbert index of a patch with given coordinates
// for a simulation box with 2^{exps[d]} patches per side
// (2^(exps[0] + exps[1]) patches in total)
uint General_Hilbert_index<2D>(const uint * exps, uint * coords)
{
    uint entry = 0;
    uint dir = 0;

    return General_Hilbert_index_orientation(exps, coords, entry, dir);
}

//============================================================================//
//  General Hilbert 3D index
//============================================================================//
// Calculates the compact Hilbert index of a patch with given coordinates
// for a simulation box with 2^{exps[d]} patches per side
// (2^(exps[0] + exps[1] + exps[2]) patches in total)
template<>
uint General_Hilbert_index<3D>(const uint * exps, const uint * coordinates)
{
    if (
        coordinates[0] >= (1 << exps[0]) || coordinates[1] >= (1 << exps[1])
        || coordinates[2] >= (1 << exps[2])
    )
    {
        return -1;
    }

    uint order[3];
                      
    // compare exponents
    Identify_descend_order(exps, order);

    uint coords[3]
        = {
            coordinates[order[0]], coordinates[order[1]], coordinates[order[2]]
        };

    const uint & min_exp = exps[order[2]];

    // approach on flattened 2D grid along max and med exps
    // 3D grid is projected along min_exp axis
    // erase last min_exp bits, not relevent for this phase
    const uint flat_exps[2]
        = { exps[order[0]] - min_exp, exps[order[1]] - min_exp };

    const uint flat_coords[2] = { coords[0] >> min_exp, coords[1] >> min_exp };

    uint entry = 0;
    uint dir = 0;

    uint index
        = General_Hilbert_index_orientation(flat_exps, flat_coords, entry, dir)
            * (1 << (3 * min_exp));

    // in local cube of side min_exp
    // local entry point "entry" and initial direction "dir" of local hilbert
    // curve has been determined by the previous call of General_Hilbert_index
    // relative position in local cube is given by last min_exp bits of position
    // only keep the last min_exp bits
    const uint mask = (1 << min_exp) - 1;

    for (uint d = 0; d < 3; ++d) { coords[d] &= mask; }

    // return overall index
    return index + Hilbert_index<3D>(min_exp, coords, entry, dir);
}

//============================================================================//
//  General Hilbert index 2D inverse
//============================================================================//
// Calculates coordinates of a patch for a given Hilbert index
// in a simulation box with 2^{exps[d]} patches
// per side (2^(exps[0] + exps[1]) patches in total)
template<>
void General_Hilbert_index_inverse<2D>(
    const uint * exps, uint * coordinates, const uint index
)
{
    // compare exponents and target coordinate which must be shifted
    uint dir = (exps[0] < exps[1]);
    uint min_exp = exps[1 - dir];
    uint ind = index;

    uint shift = 0;
    uint loc;

    // define in which sub-hypercube of side 2^{min_exp} the point is
    for (uint i = exps[0] + exps[1] - 1; i >= 2 * min_exp; --i)
    {
        loc = Get_bit(ind, i);
        shift += loc * (1 << (i - min_exp));
        ind -= loc * (1 << i);
    }

    // run the cubic inversion algorithm in the sub-hypercube
    Hilbert_index_inverse<2D>(min_exp, coordinates, ind, 0, dir);

    // shift the appropriate coordinate by the necessary value
    coordinates[dir] += shift;
}

//============================================================================//
//  General Hilbert index inverse
//============================================================================//
// Calculates coordinates of a patch for a given Hilbert index
// in a simulation box with 2^{exps[d]} patches
// per side (2^(exps[0] + exps[1] + exps[2]) patches in total)
template<>
void General_Hilbert_index_inverse<3D>(
    const uint * exps, uint * coordinates, const uint index
)
{
    uint order[3];
                      
    // compare exponents
    Identify_descend_order(exps, order);

    const uint & min_exp = exps[order[2]];

    const uint flat_exps[2]
        = { exps[order[0]] - min_exp, exps[order[1]] - min_exp };

    uint coords[3];

    // localize in which sub-hypercube the point is
    // do not account for the first 3 * min_exp bits of index
    uint loc_ind = index >> (3 * min_exp);

    // run 2D inversion algorithm on reduced domain
    General_Hilbert_index_inverse<2D>(flat_exps, coords, loc_ind);

    // now coordinates store the position of the cube in the 2D domain
    // we need to run the 3D inversion algorithm on this cube
    // with the correct entry point and direction
    uint entry = 0;
    uint dir = 0;

    // run the 2D indexgenerator in order to evaluate entry and dir
    loc_ind = General_Hilbert_index_orientation(flat_exps, coords, entry, dir);

    // store transformed coordinates in the resulting global frame
    coordinates[order[0]] = coords[0] * (1 << min_exp);
    coordinates[order[1]] = coords[1] * (1 << min_exp);
    coordinates[order[2]] = 0;

    // use only first bits of index for the local hypercube
    loc_ind = index & ((1 << (3 * min_exp)) - 1);

    // run the cubic inversion algorithm in the local sub hypercube
    Hilbert_index_inverse<3D>(min_exp, coords, loc_ind, entry, dir);

    // store the resulting coordinates
    coordinates[order[0]] += coords[0];
    coordinates[order[1]] += coords[1];
    coordinates[order[2]] += coords[2];
}

} // namespace maxwell

// hilbert.cpp
