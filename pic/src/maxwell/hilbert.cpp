// hilbert.cpp

//============================================================================//
// Get k-th Bit of i
uint_t Bit(uint_t i, uint_t k) { return (i >> k) & 1; }

// Set k-th Bit of i to value
void Set_bit(uint_t * i, uint_t k, uint_t value)
{
    *i = (*i & ~(1 << k)) | (value << k);
}

//============================================================================//
//  Bitwise rotation operators
//============================================================================//
// Works only for shift <= dim

// Evaluate left side of the rotated value then right side of it
// Finally retain only the first dim bits.
uint_t rotl(uint_t value, uint_t shift, uint_t dim)
{
    return ((value << shift) | (value >> (dim - shift))) & ((1 << dim) - 1);
}

// Evaluate right side of the rotated value then left side of it
// Finally retain only the first dim bits.
uint_t rotr(uint_t value, uint_t shift, uint_t dim)
{
    return ((value >> shift) | (value << (dim - shift))) & ((1 << dim) - 1);
}

//============================================================================//
// Generates the binary reflected Gray Code
uint_t Gray_code(uint_t i) { return i ^ (i >> 1); }

// Given a non-negative integer g,
// calculates the non-negative integer i such that Gray_code(i) = g.
uint_t Gray_code_inverse(uint_t g)
{
    uint_t i = g;

    for (uint_t j = 1; (1U << j) <= g; ++j) { i ^= (g >> j); }

    return i;
}

//============================================================================//
// TSB = trailing set bit
// It is the number of trailing set bits in the binary representation of i.
// TSB is also the inter sub-hypercube direction, g(i).
uint_t TSB(uint_t i)
{
    uint_t k = 0;

    for ( ; i & 1; ++k) { i >>= 1; }

    return k;
}

// Direction computes the sequence of intra sub-hypercube Direction,
// d(i) for 0 <= i < 2^dim.
uint_t Direction(uint_t i, uint_t dim)
{
    if (!i) { return 0; }
    else if (i & 1) { return TSB(i) % dim; }
    else { return TSB(i - 1) % dim; }
}

// Entry computes the sequence of Entry points, e(i) for 0 <= i < 2^dim.
uint_t Entry(uint_t i)
{
    if (!i) { return 0; }
    else { return Gray_code(((i - 1) >> 1) << 1); }
}

//============================================================================//
// TED is the transformation such that the Gray_code ordering of sub-hypercubes
// in the Hilbert curve defined by e and d will map to the standard binary
// reflected Gray_code.
void TED(uint_t e, uint_t d, uint_t *b, uint_t dim)
{
    *b = rotr(*b ^ e, d + 1, dim);
}

void TED_inverse(uint_t e, uint_t d, uint_t *b, uint_t dim)
{
    *b = rotl(*b, d + 1, dim) ^ e;
}

//============================================================================//
// Hilbert index2D calculates the Hilbert index h of a patch of coordinates x,y
// for a simulation box with 2^m patches per side (2^(2*m) patches in total).
uint_t Hilbert_index(uint_t m, uint_t x, uint_t y, uint_t *einit, uint_t *dinit)
{
    uint_t e = *einit;
    uint_t d = *dinit;
    uint_t h = 0;
    uint_t l;
    uint_t w;

    for (int i = m - 1; i >= 0; --i)
    {
        // i-th Bit of y at the leftmost position of l,
        // and i-th Bit of x at the rightmost position of l.
        l = (Bit(y, i) << 1) + Bit(x, i);
        TED(e, d, &l, 2);
        w = Gray_code_inverse(l);
        e = e ^ (rotl(Entry(w), d + 1, 2));
        d = (d + Direction(w, 2) + 1) % 2;
        h = (h << 2) | w;
    }

    *einit = e;
    *dinit = d;

    return h;
}

// Hilbert index3D calculates the Hilbert index h of a patch
// of coordinates x,y,z for a simulation box with 2^m patches
// per side (2^(3 * m) patches in total).
uint_t Hilbert_index(
    uint_t m, uint_t x, uint_t y, uint_t z, uint_t einit, uint_t dinit
)
{
    uint_t e = einit;
    uint_t d = dinit;
    uint_t h = 0;

    uint_t l;
    uint_t w;

    for (int i = m - 1; i >= 0; --i)
    {
        l = (Bit(z, i) << 2) + (Bit(y, i) << 1) + Bit(x, i);
        TED(e, d, &l, 3);
        w = Gray_code_inverse(l);
        e = e ^ (rotl(Entry(w), d + 1, 3));
        d = (d + Direction(w, 3) + 1) % 3;
        h = (h << 3) | w;
    }

    return h;
}

//============================================================================//
// Hilbert index2D inv  calculates the coordinates x,y of the patch
// of Hilbert index h in a simulation box with 2^m patches per side
// (2^(2*m) patches in total).
void Hilbert_index_inverse(
    uint_t m, uint_t *x, uint_t *y, uint_t h, uint_t einit, uint_t dinit
)
{
    uint_t e = einit;
    uint_t d = dinit;
    uint_t l;
    uint_t w;

    *x = 0;
    *y = 0;

    for (int i = m - 1; i >= 0; --i)
    {
        w = ((Bit(h, (i << 1) + 1)) << 1) + Bit(h, i << 1);
        l = Gray_code(w);

        TED_inverse(e, d, &l, 2);
        Set_bit(x, (uint_t)i, Bit(l, 0));
        Set_bit(y, (uint_t)i, Bit(l, 1));

        e = e ^ (rotl(Entry(w), d + 1, 2));
        d = (d + Direction(w, 2) + 1) % 2;
    }
}

void Hilbert_index_inverse(
    uint_t m,
    uint_t *x,
    uint_t *y,
    uint_t *z,
    uint_t h,
    uint_t einit,
    uint_t dinit
)
{
    uint_t e = einit;
    uint_t d = dinit;
    uint_t l;
    uint_t w;

    *x = 0;
    *y = 0;
    *z = 0;

    for (int i = m - 1; i >= 0; --i)
    {
        w = (Bit(h, 3 * i + 2) << 2) + (Bit(h, 3 * i + 1) << 1) + Bit(h, 3 * i);
        l = Gray_code(w);

        TED_inverse(e, d, &l, 3);
        Set_bit(x, (uint_t)i, Bit(l, 0));
        Set_bit(y, (uint_t)i, Bit(l, 1));
        Set_bit(z, (uint_t)i, Bit(l, 2));

        e = e ^ (rotl(Entry(w), d + 1, 3));
        d = (d + Direction(w, 3) + 1) % 3;
    }
}


//============================================================================//
//  The "general" versions of the functions
//  allow a different number of patch in each direction
//============================================================================//

// General Hilbert index2D calculates the  Hilbert index h
// of a patch of coordinates x,y for a simulation box with 2^mi patches per side
// (2^(m0 + m1)) patches in total).
uint_t generalhilbertindex(
    uint_t m0, uint_t m1, int x, int y, uint_t *einit, uint_t *dinit
)
{

    if (x < 0 || x >= (1 << m0) || y < 0 || y >= (1 << m1))
    {
        return MPI_PROC_NULL;
    }

    uint_t h = 0;
    uint_t mmin;
    uint_t mmax;
    uint_t l;
    uint_t localx = (uint_t)x;
    uint_t localy = (uint_t)y;
    uint_t *target;
    uint_t *dinit = 0;

    if (m0 >= m1)
    {
        target = &localx;
        mmin = m1;
        mmax = m0;
    }
    else
    {
        target = &localy;
        *dinit = 1;
        mmin = m0;
        mmax = m1;
    }

    for (int i = (int)mmax - 1; i >= (int)mmin; --i)
    {
        l = Bit(*target, i);
        h += l * (1 << (i + mmin));
        *target -= l * (1 << i);
    }

    if (mmin > 0) { h += Hilbert_index(mmin, localx, localy, einit, dinit); }

    return h;
}

uint_t generalhilbertindex(uint_t m0, uint_t m1, int x, int y)
{

    if (x < 0 || x >= (1 << m0) || y < 0 || y >= (1 << m1))
    {
        return MPI_PROC_NULL;
    }

    uint_t h, l, localx, localy, *target, einit, dinit;
    int mmin, mmax;
    h = 0;
    dinit = 0;
    einit = 0;
    localx = (uint_t)x;
    localy = (uint_t)y;

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
        l = Bit(*target, i);
        h += l * (1 << (i + mmin));
        *target -= l * (1 << i);
    }

    if (mmin > 0)
    {
        h += Hilbert_index((uint_t)mmin, localx, localy, &einit, &dinit);
    }

    return h;
}

// General Hilbert index3D calculates the compact Hilbert index h
// of a patch of coordinates x,y,z for a simulation box
// with 2^mi patches per side (2^(m0 + m1 + m2)) patches in total).
uint_t generalhilbertindex(uint_t m0, uint_t m1, uint_t m2, int x, int y, int z)
{
    if (
        x < 0 || x >= (1 << m0) || y < 0 || y >= (1 << m1) || z < 0
        || z >= (1 << m2)
    )
    {
        return MPI_PROC_NULL;
    }

    uint_t h, e, d, *einit, *dinit, dimmin, dimmax, dimmed, mi[3], localp[3], tempp[3];
    h = 0;
    e = 0;
    d = 0;
    dinit = &d;
    einit = &e;
    //Store positions and dimensions in arrays
    localp[0] = (uint_t)x;
    localp[1] = (uint_t)y;
    localp[2] = (uint_t)z;
    mi[0] = m0;
    mi[1] = m1;
    mi[2] = m2;

    //Compare dimension sizes
    if (m0 >= m1 && m0 >= m2) { dimmax = 0; }
    else if ((m1 > m0) && (m1 >= m2)) { dimmax = 1; }
    else { dimmax = 2; }

    if (mi[(dimmax + 1) % 3] >= mi[(dimmax + 2) % 3]) {
        dimmed = (dimmax + 1) % 3;
        dimmin = (dimmax + 2) % 3;
    }
    else { dimmed = (dimmax + 2) % 3; dimmin = (dimmax + 1) % 3; }

    // First approach on a flattened 2D grid along dimmax and dimmed.
    // The 3D grid is projected along dimmin axis.
    // Erase last mi[dimmin] bits. Not relevent for this phase.
    tempp[dimmax] = localp[dimmax] >> mi[dimmin];

    //Erase last mi[dimmin] bits. Not relevent for this phase.
    tempp[dimmed] = localp[dimmed] >> mi[dimmin];

    h += generalhilbertindex(
        mi[dimmax] - mi[dimmin], mi[dimmed] - mi[dimmin], tempp[dimmax],
        tempp[dimmed], einit, dinit
    ) * (1 << (3 * mi[dimmin]));

    // Now in a local cube of side mi[dimmin].
    // The local Entry point "einit" and initial Direction "dinit" of the local
    // hilbert curve has been determined by the previous call
    // to compacthilbertindex2.
    // Relative position in the local cube is given by the last mi[dimmin]
    // bits of the position.
    // Only keep the last mi[dimmin] bits.
    tempp[dimmax] = localp[dimmax] & ((1 << mi[dimmin]) - 1);
    //Only keep the last mi[dimmin] bits.
    tempp[dimmed] = localp[dimmed] & ((1 << mi[dimmin]) - 1);
    //Only keep the last mi[dimmin] bits.
    tempp[dimmin] = localp[dimmin] & ((1 << mi[dimmin]) - 1);

    // Add local index to the previously calculated one.
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
void generalhilbertindexinv(
    uint_t m0, uint_t m1, uint_t *x, uint_t *y, uint_t h
)
{
    uint_t einit = 0;
    uint_t dinit = 0;
    uint_t l;
    uint_t localh = h;
    uint_t * target;
    uint_t shift = 0;
    int mmin;
    int mmax;

    // Compare dimensions. Target points at the dimension which must be shifted.
    if (m0 >= m1)
    {
        target = x;
        mmin = m1;
        mmax = m0;
    }
    else
    {
        target = y;
        dinit = 1;
        mmin = m0;
        mmax = m1;
    }

    // First define in which sub-hypercube of side 2^mmin the point is.
    for (int i = mmax + mmin - 1; i >= mmin + mmin; --i)
    {
        l = Bit(localh, i);
        shift += l * (1 << (i - mmin));
        localh -= l * (1 << i);
    }

    // Run the cubic inversion algorithm in the sub hypercube.
    Hilbert_index_inverse(mmin, x, y, localh, einit, dinit);

    // Shift the appropriate coordinate by the necessary value.
    *target += shift;
}

// 3D version
void generalhilbertindexinv(
    uint_t m0, uint_t m1, uint_t m2,  uint_t *x, uint_t *y, uint_t *z, uint_t h
)
{
    uint_t e = 0;
    uint_t d = 0;
    uint_t dimmin;
    uint_t dimmax;
    uint_t dimmed;
    uint_t mi[3];
    uint_t * localp[3];
    uint_t tempp[3];
    uint_t localh;

    // Store positions and dimensions in arrays
    localp[0] = x;
    localp[1] = y;
    localp[2] = z;
    mi[0] = m0;
    mi[1] = m1;
    mi[2] = m2;

    // Compare dimension sizes
    if (m0 >= m1 && m0 >= m2) { dimmax = 0; }
    else if (m1 > m0 && m1 >= m2) { dimmax = 1; }
    else { dimmax = 2; }

    if (mi[(dimmax + 1) % 3] >= mi[(dimmax + 2) % 3])
    {
        dimmed = (dimmax + 1) % 3;
        dimmin = (dimmax + 2) % 3;
    }
    else { dimmed = (dimmax + 2) % 3; dimmin = (dimmax + 1) % 3; }

    // Localize in which sub hypercube the point is
    // Do not account for the first 3 * dimmin bits of h.
    localh = (h >> (mi[dimmin] * 3));

    //Run the 2D inversion algorithm on the reduced domain.
    generalhilbertindexinv(
        mi[dimmax] - mi[dimmin], mi[dimmed] - mi[dimmin], localp[dimmax],
        localp[dimmed], localh
    );

    // Now local P stores the position of the cube in the 2D domain
    // We need to run the 3D inversion algorithm on this cube with the correct
    // entry point and direction.
    // Run the 2D indexgenerator in order to evaluate e and d.
    localh
        = generalhilbertindex(
            mi[dimmax] - mi[dimmin], mi[dimmed] - mi[dimmin], *localp[dimmax],
            *localp[dimmed], &e, &d
        );

    // Transform coordinates in the global frame
    *localp[dimmax] *= (1 << mi[dimmin]);
    *localp[dimmed] *= (1 << mi[dimmin]);
    *localp[dimmin] = 0;

    // Use only first bits of h for the local hypercube
    localh = h & ((1 << (mi[dimmin] * 3)) - 1);

    // Run the cubic inversion algorithm in the local sub hypercube
    Hilbert_index_inverse(
        mi[dimmin], &tempp[dimmax], &tempp[dimmed], &tempp[dimmin], localh, e, d
    );

    // Add results to the coordinates
    *localp[dimmax] += tempp[dimmax];
    *localp[dimmed] += tempp[dimmed];
    *localp[dimmin] += tempp[dimmin];
}

// hilbert.cpp
