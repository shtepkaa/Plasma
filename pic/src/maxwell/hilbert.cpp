// hilbert.cpp

namespace maxwell {

//============================================================================//
//  Bit
//============================================================================//
// Get k-th Bit of i
uint Bit(const uint i, const uint k) { return (i >> k) & 1U; }

//============================================================================//
//  Set_bit
//============================================================================//
// Set k-th Bit of i to value
void Set_bit(uint & i, const uint k, const uint value)
{
    i = (i & ~(1 << k)) | (value << k);
}


//============================================================================//
//  Bitwise cyclic left rotation 
//============================================================================//
// Works only for shift <= dim
// Evaluate left side of the rotated value then right side of it
// Finally retain only the first dim bits
uint Rotate_left(const uint value, const uint shift, const uint dim)
{
    return (value << shift | value >> (dim - shift)) & ((1 << dim) - 1);
}

//============================================================================//
//  Bitwise cyclic right rotation 
//============================================================================//
// Works only for shift <= dim
// Evaluate right side of the rotated value then left side of it
// Finally retain only the first dim bits
uint Rotate_right(const uint value, const uint shift, const uint dim)
{
    return (value >> shift | value << (dim - shift)) & ((1 << dim) - 1);
}

//============================================================================//
//  Gray_code
//============================================================================//
// Generates the binary reflected Gray Code
uint Gray_code(const uint i) { return i ^ (i >> 1); }

//============================================================================//
//  Gray_code_inverse
//============================================================================//
// Given a non-negative integer g,
// calculates the non-negative integer i such that Gray_code(i) = g
uint Gray_code_inverse(const uint g)
{
    uint i = g;

    for (uint j = 1; (1U << j) <= g; ++j) { i ^= g >> j; }

    return i;
}

//============================================================================//
//  Trailing_set_bit
//============================================================================//
// It is the number of trailing set bits in the binary representation of i
// Trailing_set_bit is also the inter sub-hypercube direction, g(i)
uint Trailing_set_bit(const uint i)
{
    uint k = 0;

    for (uint j = i; j & 1; j >>= 1) { ++k; }

    return k;
}

//============================================================================//
//  Direction
//============================================================================//
// computes the sequence of intra sub-hypercube Direction,
// d(i) for 0 <= i < 2^dim
uint Direction(const uint i, const uint dim)
{
    if (!i) { return 0; }
    else if (i & 1) { return Trailing_set_bit(i) % dim; }
    else { return Trailing_set_bit(i - 1) % dim; }
}

//============================================================================//
//  Entry
//============================================================================//
// computes the sequence of Entry points, e(i) for 0 <= i < 2^dim
uint Entry(const uint i)
{
    if (i) { return Gray_code(((i - 1) >> 1) << 1); }
    else { return 0; }
}

//============================================================================//
//  Transform
//============================================================================//
// is the transformation such that
// the Gray_code ordering of sub-hypercubes in the Hilbert curve
// defined by "entry" and "dir" will map to the standard binary
// reflected Gray_code
void Transform(const uint entry, const uint dir, uint & b, const uint dim)
{
    b = Rotate_right(b ^ entry, dir + 1, dim);
}

void Transform_inverse(
    const uint entry, const uint dir, uint & b, const uint dim
)
{
    b = Rotate_left(b, dir + 1, dim) ^ entry;
}

//============================================================================//
//  Hilbert index with orientation
//============================================================================//
// 2D
// calculates the Hilbert index index of patch of given coordinates
// for simulation box with 2^{dim} patches per side (2^(2*dim) patches in total)
uint Hilbert_index_orientation(
    const uint dim, const uint * coords, uint & entry, uint & dir
)
{
    uint index = 0;

    uint loc;
    uint loc_ind;

    for (uint i = dim - 1; i >= 0; --i)
    {
        // i-th bit of coords[1] at the leftmost position of loc,
        // and i-th bit of coords[0] at the rightmost position of loc
        loc = (Bit(coords[1], i) << 1) + Bit(coords[0], i);
        Transform(entry, dir, loc, 2);

        loc_ind = Gray_code_inverse(loc);

        entry ^= Rotate_left(Entry(loc_ind), dir + 1, 2);
        dir = (dir + Direction(loc_ind, 2) + 1) & 1;
        index = (index << 2) | loc_ind;
    }

    return index;
}

//============================================================================//
//  Hilbert index
//============================================================================//
// 3D
// calculates the Hilbert index of patch of given coordinates
// for a simulation box with 2^{dim} patches
// per side (2^(3 * dim) patches in total)
uint Hilbert_index(
    const uint dim, const uint * coords, const uint entry, const uint dir
)
{
    uint index = 0;

    uint loc;
    uint loc_ind;

    for (uint i = dim - 1; i >= 0; --i)
    {
        loc = (Bit(coords[2], i) << 2) + (Bit(coords[1], i) << 1)
            + Bit(coords[0], i);

        Transform(entry, dir, loc, 3);
        loc_ind = Gray_code_inverse(loc);

        entry = entry ^ (Rotate_left(Entry(loc_ind), dir + 1, 3));
        dir = (dir + Direction(loc_ind, 3) + 1) % 3;
        index = (index << 3) | loc_ind;
    }

    return index;
}

//============================================================================//
//  Hilbert index
//============================================================================//
// calculates the coordinates of patch of given Hilbert index
// in a simulation box with 2^{dim} patches per side
// (2^(2*dim) patches in total)
void Hilbert_index_inverse(
    const uint dim,
    uint * coords,
    const uint index,
    const uint entry,
    const uint dir
)
{
    uint loc;
    uint loc_ind;

    coords[0] = coords[1] = 0;

    for (uint i = dim - 1; i >= 0; --i)
    {
        loc_ind = (Bit(index, (i << 1) + 1) << 1) + Bit(index, i << 1);
        loc = Gray_code(loc_ind);
        Transform_inverse(entry, dir, loc, 2);

        Set_bit(coords[0], i, Bit(loc, 0));
        Set_bit(coords[1], i, Bit(loc, 1));

        entry ^= Rotate_left(Entry(loc_ind), dir + 1, 2);
        dir = (dir + Direction(loc_ind, 2) + 1) & 1U;
    }
}

//============================================================================//
//  Hilbert index inverse
//============================================================================//
// 3D
void Hilbert_index_inverse(
    const uint dim,
    uint * coords,
    const uint index,
    const uint entry,
    const uint dir
)
{
    uint loc;
    uint loc_ind;

    coords[0] = coords[1] = coords[2] = 0;

    for (uint i = dim - 1; i >= 0; --i)
    {
        loc_ind = (Bit(index, 3 * i + 2) << 2) + (Bit(index, 3 * i + 1) << 1)
            + Bit(index, 3 * i);

        loc = Gray_code(loc_ind);
        Transform_inverse(entry, dir, loc, 3);

        Set_bit(coords[0], i, Bit(loc, 0));
        Set_bit(coords[1], i, Bit(loc, 1));
        Set_bit(coords[2], i, Bit(loc, 2));

        entry ^= Rotate_left(Entry(loc_ind), dir + 1, 3);
        dir = (dir + Direction(loc_ind, 3) + 1) % 3;
    }
}

//============================================================================//
//  General Hilbert index
//============================================================================//
// calculates Hilbert index of patch of given coordinates
// for a simulation box with 2^{dims[i]} patches per side
// (2^{dims[0] + dims[1]} patches in total)
uint General_hilbert_index_orientation(
    const uint * dims, const int * coords, uint & entry, uint & dir
)
{
    if (
        coords[0] < 0 || coords[0] >= (1 << dims[0])
        || coords[1] < 0 || coords[1] >= (1 << dims[1])
    )
    {
        return MPI_PROC_NULL;
    }

    const uint crds[2] = { uint(coords[0]), uint(coords[1]) }; 

    dir = (dims[0] < dims[1]);

    uint index = 0;
    uint min_dim = dims[1 - dir];

    uint loc;

    for (uint i = dims[dir] - 1; i >= min_dim; --i)
    {
        loc = Bit(coords[dir], i);
        index += loc * (1 << (i + min_dim));
        coords[dir] -= loc * (1 << i);
    }

    // calculate entry and direction
    if (min_dim)
    {
        index += Hilbert_index_orientation(min_dim, crds, entry, dir);
    }

    return index;
}

//============================================================================//
//  General Hilbert index
//============================================================================//
uint General_hilbert_index(const uint * dims, int * coords)
{
    uint entry = 0;
    uint dir = 0;

    return General_hilbert_index_orientation(dims, coords, entry, dir);
}


template<>
void Arrange_array<2>(const uint * in, uint * out)
{
    const uint imax = (in[0] < in[1]);

    out[0] = in(imax);
    out[1] = in(1 - imax);
}

template<>
void Arrange_array<3>(const uint * in, uint * out)
{
    const uint comp[3] = { (in[1] < in[2]), (in[0] < in[2]), (in[0] < in[1]) };

    const uint imed = comp[0] ^ comp[1] ^ comp[2];

    comp[imed] <<= imed & 1;

    out[0] = in(1 - (imed > 0) + comp[imed]);
    out[1] = in(imed);
    out[2] = in(1 + (imed < 2) - comp[imed]);
}

//============================================================================//
//  General Hilbert index
//============================================================================//
// Calculates the compact Hilbert index of a patch of given coordinates
// for a simulation box with 2^{dims[i]} patches per side
// (2^(dims[0] + dims[1] + dims[2]) patches in total)
uint General_hilbert_index(const uint * dims, int * coords)
{
    if (
        coords[0] < 0 || coords[0] >= (1 << dims[0])
        || coords[1] < 0 || coords[1] >= (1 << dims[1])
        || coords[2] < 0 || coords[2] >= (1 << dims[2])
    )
    {
        return MPI_PROC_NULL;
    }

    uint crds[3] = { uint(coords[0]), uint(coords[1]), uint(coords[2]) };

    uint index = 0;
    uint entry = 0;
    uint dir = 0;

    uint imin;
    uint imax;
    uint imed;

    uint arraged_dims[3];
                      
    // compare dimension sizes
    if (dims[0] >= dims[1] && dims[0] >= dims[2]) { imax = 0; }
    else if ((dims[1] > dims[0]) && (dims[1] >= dims[2])) { imax = 1; }
    else { imax = 2; }

    if (dims[(imax + 1) % 3] >= dims[(imax + 2) % 3])
    {
        imed = (imax + 1) % 3;
        imin = (imax + 2) % 3;
    }
    else { imed = (imax + 2) % 3; imin = (imax + 1) % 3; }

    // approach on a flattened 2D grid along imax and imed
    // the 3D grid is projected along imin axis
    // erase last dims[imin] bits. Not relevent for this phase
    tempp[imax] = crds[imax] >> dims[imin];

    // erase last dims[imin] bits. Not relevent for this phase
    tempp[imed] = crds[imed] >> dims[imin];

    index += General_hilbert_index(
        dims[imax] - dims[imin], dims[imed] - dims[imin], tempp[imax],
        tempp[imed], entry, dir
    ) * (1 << (3 * dims[imin]));

    // In local cube of side dims[imin]
    // The local Entry point "entry" and initial Direction "dir" of the local
    // hilbert curve has been determined by the previous call
    // to compacthilbertindex2
    // Relative position in the local cube is given by the last dims[imin]
    // bits of the position

    // Only keep the last dims[imin] bits
    tempp[imax] = crds[imax] & ((1 << dims[imin]) - 1);
    tempp[imed] = crds[imed] & ((1 << dims[imin]) - 1);
    tempp[imin] = crds[imin] & ((1 << dims[imin]) - 1);

    // Add local index to the previously calculated one
    index
        += Hilbert_index(
            dims[imin], tempp[imax], tempp[imed], tempp[imin], entry, dir
        );

    return index;
}

//============================================================================//
// 2D version
//============================================================================//
// General Hilbert index inverse calculates given coordinates
// of a patch for a given Hilbert index in a simulation box with 2^{dims[i]} patches
// per side (2^(dims[0] + dims[1]) patches in total)
void General_hilbert_index_inverse(
    const uint * dims, uint * coords, const uint index
)
{
    // compare dimensions and target coordinate which must be shifted
    uint dir = (dims[0] < dims[1]);
    uint min_dim = dims[1 - dir];
    uint ind = index;

    uint shift = 0;
    uint loc;

    // define in which sub-hypercube of side 2^{min_dim} the point is
    for (uint i = dims[0] + dims[1] - 1; i >= 2 * min_dim; --i)
    {
        loc = Bit(ind, i);
        shift += loc * (1 << (i - min_dim));
        ind -= loc * (1 << i);
    }

    // run the cubic inversion algorithm in the sub-hypercube
    Hilbert_index_inverse(min_dim, coords, ind, 0, dir);

    // shift the appropriate coordinate by the necessary value
    coords[dir] += shift;
}

//============================================================================//
// 3D version
//============================================================================//
void General_hilbert_index_inverse(const uint * dims, uint * coords, uint index)
{
    uint entry = 0;
    uint dir = 0;

    uint imin;
    uint imed;
    uint imax;

    uint localh;

    uint tempp[3];

    // compare dimension sizes
    if (dims[0] >= dims[1] && dims[0] >= dims[2]) { imax = 0; }
    else if (dims[1] > dims[0] && dims[1] >= dims[2]) { imax = 1; }
    else { imax = 2; }

    if (dims[(imax + 1) % 3] >= dims[(imax + 2) % 3])
    {
        imed = (imax + 1) % 3;
        imin = (imax + 2) % 3;
    }
    else { imed = (imax + 2) % 3; imin = (imax + 1) % 3; }

    // localize in which sub hypercube the point is
    // do not account for the first 3 * imin bits of index
    localh = (index >> (dims[imin] * 3));

    // run the 2D inversion algorithm on the reduced domain
    General_hilbert_index_inverse(
        dims[imax] - dims[imin], dims[imed] - dims[imin],
        dims[imax], dims[imed], localh
    );

    // now local P stores the position of the cube in the 2D domain
    // we need to run the 3D inversion algorithm on this cube
    // with the correct entry point and direction

    // run the 2D indexgenerator in order to evaluate entry and dir
    localh
        = General_hilbert_index_orientation(
            dims[imax] - dims[imin], dims[imed] - dims[imin],
            dims[imax], dims[imed], entry, dir
        );

    // transform coordinates in the global frame
    dims[imax] *= (1 << dims[imin]);
    dims[imed] *= (1 << dims[imin]);
    dims[imin] = 0;

    // use only first bits of index for the local hypercube
    localh = index & ((1 << (dims[imin] * 3)) - 1);

    // run the cubic inversion algorithm in the local sub hypercube
    Hilbert_index_inverse(
        dims[imin], &tempp[imax], &tempp[imed], &tempp[imin], localh,
        entry, dir
    );

    // add results to the coordinates
    dims[imax] += tempp[imax];
    dims[imed] += tempp[imed];
    dims[imin] += tempp[imin];
}

} // namespace maxwell

// hilbert.cpp
