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
// set k-th bit of i to value
void Set_bit(uint & i, const uint k, const uint value)
{
    i = (i & ~(1 << k)) | (value << k);
}

//============================================================================//
//  Bitwise cyclic left rotation 
//============================================================================//
// retain only the first dim bits
// works only for (shift <= dim)
uint Rotate_left(const uint value, const uint shift, const uint dim)
{
    return (value << shift | value >> (dim - shift)) & ((1 << dim) - 1);
}

//============================================================================//
//  Bitwise cyclic right rotation 
//============================================================================//
// retain only the first dim bits
// works only for (shift <= dim)
uint Rotate_right(const uint value, const uint shift, const uint dim)
{
    return (value >> shift | value << (dim - shift)) & ((1 << dim) - 1);
}

//============================================================================//
//  Gray code inverse
//============================================================================//
// calculate the non-negative integer i such that Gray_code(i) = g,
// for non-negative integer g
uint Gray_code_inverse(const uint g)
{
    uint i = g;

    for (uint j = 1; (1U << j) <= g; ++j) { i ^= g >> j; }

    return i;
}

//============================================================================//
//  Trailing set bit
//============================================================================//
// compute number of trailing set bits in the binary representation of i,
// which is also the inter sub-hypercube direction
uint Trailing_set_bit(const uint i)
{
    uint k = 0;

    for (uint j = i; j & 1; j >>= 1) { ++k; }

    return k;
}

//============================================================================//
//  Direction
//============================================================================//
// compute intra sub-hypercube direction for (0 <= i < 2^{dim})
uint Direction(const uint i, const uint dim)
{
    if (!i) { return 0; }
    else if (i & 1) { return Trailing_set_bit(i) % dim; }
    else { return Trailing_set_bit(i - 1) % dim; }
}

//============================================================================//
//  Entry
//============================================================================//
// compute entry point for (0 <= i < 2^{dim})
uint Entry(const uint i)
{
    if (i) { return Gray_code(((i - 1) >> 1) << 1); }
    else { return 0; }
}

//============================================================================//
//  Transform
//============================================================================//
// compute transformation such that the Gray code ordering of sub-hypercubes
// in the Hilbert curve defined by "entry" and "dir"
// will map to the standard binary reflected Gray code
void Transform(const uint entry, const uint dir, uint & bin, const uint dim)
{
    bin = Rotate_right(bin ^ entry, dir + 1, dim);
}

void Transform_inverse(
    const uint entry, const uint dir, uint & bin, const uint dim
)
{
    bin = Rotate_left(bin, dir + 1, dim) ^ entry;
}

//============================================================================//
//  Hilbert index 2D with orientation
//============================================================================//
// calculate the Hilbert index of patch of given coordinates for simulation box
// with 2^{dim} patches per side (2^(2 * dim) patches in total)
uint Hilbert_index_orientation(
    const uint dim, const uint * coordinates, uint & entry, uint & dir
)
{
    uint index = 0;

    uint loc;
    uint loc_ind;

    for (uint i = dim - 1; i >= 0; --i)
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
template<>
uint Hilbert_index<2D>(const uint dim, const uint * coordinates)
{
    uint entry = 0;
    uint dir = 0;

    return Hilbert_index_orientation(dim, coordinates, entry, dir);
}

//============================================================================//
//  Hilbert index 3D
//============================================================================//
// calculates the Hilbert index of patch of given coordinates
// for a simulation box with 2^{dim} patches
// per side (2^(3 * dim) patches in total)
template<>
uint Hilbert_index<3D>(
    const uint dim, const uint * coordinates, const uint entry, const uint dir
)
{
    uint index = 0;

    uint loc;
    uint loc_ind;

    for (uint i = dim - 1; i >= 0; --i)
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
// calculates the coordinates of patch of given Hilbert index
// in a simulation box with 2^{dim} patches per side
// (2^(2 * dim) patches in total)
template<>
void Hilbert_index_inverse<2D>(
    const uint dim,
    uint * coordinates,
    const uint index,
    const uint entry,
    const uint dir
)
{
    uint loc;
    uint loc_ind;

    coordinates[0] = coordinates[1] = 0;

    for (uint i = dim - 1; i >= 0; --i)
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
template<>
void Hilbert_index_inverse<3D>(
    const uint dim,
    uint * coordinates,
    const uint index,
    const uint entry,
    const uint dir
)
{
    uint loc;
    uint loc_ind;

    coordinates[0] = coordinates[1] = coordinates[2] = 0;

    for (uint i = dim - 1; i >= 0; --i)
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
// calculates Hilbert index of patch of given coordinates
// for a simulation box with 2^{dims[i]} patches per side
// (2^{dims[0] + dims[1]} patches in total)
uint General_hilbert_index_orientation(
    const uint * dims, const int * coordinates, uint & entry, uint & dir
)
{
    if (
        coordinates[0] < 0 || coordinates[0] >= (1 << dims[0])
        || coordinates[1] < 0 || coordinates[1] >= (1 << dims[1])
    )
    {
        return MPI_PROC_NULL;
    }

    const uint coords[2] = { uint(coordinates[0]), uint(coordinates[1]) }; 

    dir = (dims[0] < dims[1]);

    uint index = 0;
    uint min_dim = dims[1 - dir];

    uint loc;

    for (uint i = dims[dir] - 1; i >= min_dim; --i)
    {
        loc = Get_bit(coords[dir], i);
        index += loc * (1 << (i + min_dim));
        coords[dir] -= loc * (1 << i);
    }

    // calculate entry and direction
    if (min_dim)
    {
        index += Hilbert_index_orientation(min_dim, coords, entry, dir);
    }

    return index;
}

//============================================================================//
//  General Hilbert 2D index
//============================================================================//
uint General_hilbert_index<2D>(const uint * dims, int * coords)
{
    uint entry = 0;
    uint dir = 0;

    return General_hilbert_index_orientation(dims, coords, entry, dir);
}

//============================================================================//
//  General Hilbert 3D index
//============================================================================//
// Calculates the compact Hilbert index of a patch of given coordinates
// for a simulation box with 2^{dims[i]} patches per side
// (2^(dims[0] + dims[1] + dims[2]) patches in total)
template<>
uint General_hilbert_index<3D>(const uint * dims, const int * coordinates)
{
    if (
        coordinates[0] < 0 || coordinates[0] >= (1 << dims[0])
        || coordinates[1] < 0 || coordinates[1] >= (1 << dims[1])
        || coordinates[2] < 0 || coordinates[2] >= (1 << dims[2])
    )
    {
        return MPI_PROC_NULL;
    }

    uint order[3];
                      
    // compare dimension sizes
    Identify_descend_order(dims, order);

    uint coords[3]
        = {
            uint(coordinates[order[0]]),
            uint(coordinates[order[1]]),
            uint(coordinates[order[2]])
        };

    const uint min_dim = dims[order[2]];

    // approach on flattened 2D grid along max and med dims
    // 3D grid is projected along min_dim axis
    // erase last min_dim bits, not relevent for this phase
    const uint flat_dims[2]
        = { dims[order[0]] - min_dim, dims[order[1]] - min_dim };

    const uint flat_coords[2] = { coords[0] >> min_dim, coords[1] >> min_dim };

    uint entry = 0;
    uint dir = 0;

    uint index
        = General_hilbert_index_orientation(flat_dims, flat_coords, entry, dir)
            * (1 << (3 * min_dim));

    // in local cube of side min_dim
    // local entry point "entry" and initial direction "dir" of local hilbert
    // curve has been determined by the previous call of General_hilbert_index
    // relative position in local cube is given by last min_dim bits of position
    // only keep the last min_dim bits
    const uint mask = (1 << min_dim) - 1;

    for (int d = 0; d < 3; ++d) { coords[d] &= mask; }

    // return overall index
    return index + Hilbert_index<3D>(min_dim, coords, entry, dir);
}

//============================================================================//
//  General Hilbert index 2D inverse
//============================================================================//
// General Hilbert index inverse calculates given coordinates of patch
// for given Hilbert index in a simulation box with 2^{dims[i]} patches
// per side (2^(dims[0] + dims[1]) patches in total)
template<>
void General_hilbert_index_inverse<2D>(
    const uint * dims, uint * coordinates, const uint index
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
        loc = Get_bit(ind, i);
        shift += loc * (1 << (i - min_dim));
        ind -= loc * (1 << i);
    }

    // run the cubic inversion algorithm in the sub-hypercube
    Hilbert_index_inverse<2D>(min_dim, coordinates, ind, 0, dir);

    // shift the appropriate coordinate by the necessary value
    coordinates[dir] += shift;
}

//============================================================================//
//  General Hilbert index inverse
//============================================================================//
// 3D version
template<>
void General_hilbert_index_inverse<3D>(
    const uint * dims, uint * coordinates, const uint index
)
{
    uint order[3];
                      
    // compare dimension sizes
    Identify_descend_order(dims, order);

    const uint min_dim = dims[order[2]];

    const uint flat_dims[2]
        = { dims[order[0]] - min_dim, dims[order[1]] - min_dim };

    uint coords[3];

    // localize in which sub-hypercube the point is
    // do not account for the first 3 * min_dim bits of index
    uint loc_ind = index >> (3 * min_dim);

    // run 2D inversion algorithm on reduced domain
    General_hilbert_index_inverse<2D>(flat_dims, coords, loc_ind);

    // now coordinates store the position of the cube in the 2D domain
    // we need to run the 3D inversion algorithm on this cube
    // with the correct entry point and direction
    uint entry = 0;
    uint dir = 0;

    // run the 2D indexgenerator in order to evaluate entry and dir
    loc_ind = General_hilbert_index_orientation(flat_dims, coords, entry, dir);

    // store transformed coordinates in the resulting global frame
    coordinates[order[0]] = coords[0] * (1 << min_dim);
    coordinates[order[1]] = coords[1] * (1 << min_dim);
    coordinates[order[2]] = 0;

    // use only first bits of index for the local hypercube
    loc_ind = index & ((1 << (3 * min_dim)) - 1);

    // run the cubic inversion algorithm in the local sub hypercube
    Hilbert_index_inverse<3D>(min_dim, coords, loc_ind, entry, dir);

    // store the resulting coordinates
    coordinates[order[0]] += coords[0];
    coordinates[order[1]] += coords[1];
    coordinates[order[2]] += coords[2];
}

} // namespace maxwell

// hilbert.cpp
