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

    for (uint j = 1U; (1U << j) <= g; ++j) { i ^= g >> j; }

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
//  Hilbert index2D
//============================================================================//
// calculates the Hilbert index h of a patch of coordinates x,y
// for a simulation box with 2^m patches per side (2^(2*m) patches in total)
uint Hilbert_index_2D(
    const uint m, const uint * coords, uint & entry, uint & dir
)
{
    uint index = 0;

    uint l;
    uint w;

    for (uint i = m - 1; i >= 0; --i)
    {
        // i-th Bit of y at the leftmost position of l,
        // and i-th Bit of x at the rightmost position of l
        l = (Bit(coords[1], i) << 1) + Bit(coords[0], i);
        Transform(entry, dir, l, 2);
        w = Gray_code_inverse(l);

        entry ^= Rotate_left(Entry(w), dir + 1, 2);
        dir = (dir + Direction(w, 2) + 1) & 1;
        index = (index << 2) | w;
    }

    return index;
}

//============================================================================//
//  Hilbert index3D
//============================================================================//
// Hilbert index3D calculates the Hilbert index h of a patch
// of coordinates x,y,z for a simulation box with 2^m patches
// per side (2^(3 * m) patches in total)
uint Hilbert_index(
    const uint m, const uint * coords, const uint entry, const uint dir
)
{
    uint index = 0;

    uint l;
    uint w;

    for (uint i = m - 1; i >= 0; --i)
    {
        l = (Bit(coords[2], i) << 2) + (Bit(coords[1], i) << 1)
            + Bit(coords[0], i);

        Transform(entry, dir, l, 3);
        w = Gray_code_inverse(l);

        entry = entry ^ (Rotate_left(Entry(w), dir + 1, 3));
        dir = (dir + Direction(w, 3) + 1) % 3;
        index = (index << 3) | w;
    }

    return index;
}

//============================================================================//
//  Hilbert index2D inv
//============================================================================//
// calculates the coordinates x,y of the patch
// of Hilbert index h in a simulation box with 2^m patches per side
// (2^(2*m) patches in total)
void Hilbert_index_inverse(
    const uint m,
    uint * coords,
    const uint index,
    const uint entry,
    const uint dir
)
{
    uint l;
    uint w;

    coords[0] = coords[1] = 0;

    for (uint i = m - 1; i >= 0; --i)
    {
        w = (Bit(index, (i << 1) + 1) << 1) + Bit(index, i << 1);
        l = Gray_code(w);

        Transform_inverse(entry, dir, l, 2);
        Set_bit(coords[0], i, Bit(l, 0));
        Set_bit(coords[1], i, Bit(l, 1));

        entry ^= Rotate_left(Entry(w), dir + 1, 2);
        dir = (dir + Direction(w, 2) + 1) & 1U;
    }
}

//============================================================================//
//  Hilbert index3D inv
//============================================================================//
void Hilbert_index_inverse(
    const uint m,
    uint * coords,
    const uint index,
    const uint entry,
    const uint dir
)
{
    uint l;
    uint w;

    coords[0] = coords[1] = coords[2] = 0;

    for (uint i = m - 1; i >= 0; --i)
    {
        w = (Bit(index, 3 * i + 2) << 2) + (Bit(index, 3 * i + 1) << 1)
            + Bit(index, 3 * i);

        l = Gray_code(w);
        Transform_inverse(entry, dir, l, 3);

        Set_bit(coords[0], i, Bit(l, 0));
        Set_bit(coords[1], i, Bit(l, 1));
        Set_bit(coords[2], i, Bit(l, 2));

        entry ^= Rotate_left(Entry(w), dir + 1, 3);
        dir = (dir + Direction(w, 3) + 1) % 3;
    }
}

//============================================================================//
//  The "general" versions of the functions
//  allow a different number of patch in each direction
//============================================================================//

// General Hilbert index2D calculates the  Hilbert index h
// of a patch of coordinates x,y for a simulation box with 2^mi patches per side
// (2^(m0 + m1)) patches in total)
uint General_hilbert_index(
    const uint * dims, const int * coords, uint * entry, uint * dir
)
{
    if (x < 0 || x >= (1 << m0) || y < 0 || y >= (1 << m1))
    {
        return MPI_PROC_NULL;
    }

    uint index = 0;

    uint mmin;
    uint mmax;

    uint l;
    uint localx = (uint)x;
    uint localy = (uint)y;
    uint * target;
    uint * dinit = 0;

    if (m0 >= m1)
    {
        target = &localx;
        mmin = dims[1];
        mmax = dims[0];
    }
    else
    {
        target = &localy;
        *dinit = 1;
        mmin = dims[0];
        mmax = dims[1];
    }

    for (uint i = mmax - 1; i >= mmin; --i)
    {
        l = Bit(*target, i);
        index += l * (1 << (i + mmin));
        *target -= l * (1 << i);
    }

    if (mmin > 0)
    {
        index += Hilbert_index(mmin, localx, localy, entry, dir);
    }

    return index;
}

uint General_hilbert_index(uint m0, uint m1, int x, int y)
{
    if (x < 0 || x >= (1 << m0) || y < 0 || y >= (1 << m1))
    {
        return MPI_PROC_NULL;
    }

    uint h = 0;
    uint l;
    uint localx = uint(x);
    uint localy = uint(y);
    uint einit = 0;
    uint dinit = 0;
    uint * target;

    int mmin;
    int mmax;

    if (m0 >= m1)
    {
        target = &localx;
        mmin = m1;
        mmax = m0;
    }
    else
    {
        target = &localy;
        dinit = 1;
        mmin = m0;
        mmax = m1;
    }

    for (int i = mmax - 1; i >= mmin; --i)
    {
        l = bit(*target, i);
        h += l * (1 << (i + mmin));
        *target -= l * (1 << i);
    }

    if (mmin > 0)
    {
        h += hilbert_index((uint)mmin, localx, localy, &einit, &dinit);
    }

    return h;
}

// General Hilbert index3D calculates the compact Hilbert index h
// of a patch of coordinates x,y,z for a simulation box
// with 2^mi patches per side (2^(m0 + m1 + m2)) patches in total)
uint General_hilbert_index(uint m0, uint m1, uint m2, int x, int y, int z)
{
    if (
        x < 0 || x >= (1 << m0) || y < 0 || y >= (1 << m1) || z < 0
        || z >= (1 << m2)
    )
    {
        return MPI_PROC_NULL;
    }

    uint h = 0;
    uint e = 0;
    uint d = 0;
    uint dimmin;
    uint dimmax;
    uint dimmed;
    uint * einit = &e;
    uint * dinit = &d;

    uint mi[3];
    uint localp[3];
    uint tempp[3];

    //Store positions and dimensions in arrays
    localp[0] = uint(x);
    localp[1] = uint(y);
    localp[2] = uint(z);
    mi[0] = m0;
    mi[1] = m1;
    mi[2] = m2;

    //Compare dimension sizes
    if (m0 >= m1 && m0 >= m2) { dimmax = 0; }
    else if ((m1 > m0) && (m1 >= m2)) { dimmax = 1; }
    else { dimmax = 2; }

    if (mi[(dimmax + 1) % 3] >= mi[(dimmax + 2) % 3])
    {
        dimmed = (dimmax + 1) % 3;
        dimmin = (dimmax + 2) % 3;
    }
    else { dimmed = (dimmax + 2) % 3; dimmin = (dimmax + 1) % 3; }

    // First approach on a flattened 2D grid along dimmax and dimmed
    // The 3D grid is projected along dimmin axis
    // Erase last mi[dimmin] bits. Not relevent for this phase
    tempp[dimmax] = localp[dimmax] >> mi[dimmin];

    //Erase last mi[dimmin] bits. Not relevent for this phase
    tempp[dimmed] = localp[dimmed] >> mi[dimmin];

    h += General_hilbert_index(
        mi[dimmax] - mi[dimmin], mi[dimmed] - mi[dimmin], tempp[dimmax],
        tempp[dimmed], einit, dinit
    ) * (1 << (3 * mi[dimmin]));

    // Now in a local cube of side mi[dimmin]
    // The local Entry point "einit" and initial Direction "dinit" of the local
    // hilbert curve has been determined by the previous call
    // to compacthilbertindex2
    // Relative position in the local cube is given by the last mi[dimmin]
    // bits of the position
    // Only keep the last mi[dimmin] bits
    tempp[dimmax] = localp[dimmax] & ((1 << mi[dimmin]) - 1);
    //Only keep the last mi[dimmin] bits
    tempp[dimmed] = localp[dimmed] & ((1 << mi[dimmin]) - 1);
    //Only keep the last mi[dimmin] bits
    tempp[dimmin] = localp[dimmin] & ((1 << mi[dimmin]) - 1);

    // Add local index to the previously calculated one
    h += Hilbert_index(
        mi[dimmin], tempp[dimmax], tempp[dimmed], tempp[dimmin], *einit, *dinit
    );

    return h;
}

//============================================================================//
// General Hilbert index inverse calculates the coordinates x,y
// of a patch for a given Hilbert index h in a simulation box with 2^mi patches
// per side (2^(m0 + m1) patches in total)
// 2D version
void General_hilbert_index_inverse(
    const uint * dims, uint * coords, const uint index
)
{
    uint entry = 0;
    uint dir = 0;
    uint ind = index;

    uint l;
    uint shift = 0;
    uint * target;

    uint mmin;
    uint mmax;

    // Compare dimensions. Target points at the dimension which must be shifted
    if (m0 >= m1)
    {
        target = x;
        mmin = dims[1];
        mmax = dims[0];
    }
    else
    {
        target = y;
        dir = 1;
        mmin = dims[0];
        mmax = dims[1];
    }

    // First define in which sub-hypercube of side 2^mmin the point is
    for (uint i = mmax + mmin - 1; i >= mmin + mmin; --i)
    {
        l = Bit(ind, i);
        shift += l * (1 << (i - mmin));
        ind -= l * (1 << i);
    }

    // Run the cubic inversion algorithm in the sub hypercube
    Hilbert_index_inverse(mmin, x, y, ind, entry, dir);

    // Shift the appropriate coordinate by the necessary value
    *target += shift;
}

// 3D version
void General_hilbert_index_inverse(const uint * dims, uint * coords, uint index)
{
    uint entry = 0;
    uint dir = 0;

    uint imin;
    uint imed;
    uint imax;

    uint localh;

    uint tempp[3];

    uint * localp[3];

    // store positions and dimensions in arrays
    localp[0] = x;
    localp[1] = y;
    localp[2] = z;

    // compare dimension sizes
    if (m0 >= m1 && m0 >= m2) { imax = 0; }
    else if (m1 > m0 && m1 >= m2) { imax = 1; }
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
        localp[imax], localp[imed], localh
    );

    // now local P stores the position of the cube in the 2D domain
    // we need to run the 3D inversion algorithm on this cube
    // with the correct entry point and direction

    // run the 2D indexgenerator in order to evaluate entry and dir
    localh
        = General_hilbert_index(
            dims[imax] - dims[imin], dims[imed] - dims[imin],
            *localp[imax], *localp[imed], &entry, &dir
        );

    // transform coordinates in the global frame
    *localp[imax] *= (1 << dims[imin]);
    *localp[imed] *= (1 << dims[imin]);
    *localp[imin] = 0;

    // use only first bits of index for the local hypercube
    localh = index & ((1 << (dims[imin] * 3)) - 1);

    // run the cubic inversion algorithm in the local sub hypercube
    Hilbert_index_inverse(
        dims[imin], &tempp[imax], &tempp[imed], &tempp[imin], localh,
        entry, dir
    );

    // add results to the coordinates
    *localp[imax] += tempp[imax];
    *localp[imed] += tempp[imed];
    *localp[imin] += tempp[imin];
}

} // namespace maxwell

// hilbert.cpp
