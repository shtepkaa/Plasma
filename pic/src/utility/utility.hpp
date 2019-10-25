#ifndef UTILITY_HPP
#define UTILITY_HPP

namespace maxwell {

////////////////////////////////////////////////////////////////////////////////
//
//  BLAS-like generic functions
//
////////////////////////////////////////////////////////////////////////////////
//==============================================================================
//  Copy
//==============================================================================
template<typename Type>
/* ??? __host__ __device__*/ void Copy(
    const uint dim,
    const Type * src,
    const uint inc_src,
    Type * dst,
    const uint inc_dst
)
{
/* #pragma omp parallel for */
    for (int s = 0, d = 0; d < dim; s += inc_src, d += inc_dst)
    {
        dst[d] = src[s];
    }
}

//==============================================================================
//  Axpy
//==============================================================================
template<typename Type>
/* ??? __host__ __device__*/ void Axpy(
    const uint dim,
    const Type & mult,
    const Type * src,
    const uint inc_src,
    Type * dst,
    const uint inc_dst
)
{
/* #pragma omp parallel for */
    for (int s = 0, d = 0; d < dim; s += inc_src, d += inc_dst)
    {
        dst[d] += mult * src[s];
    }
}

////////////////////////////////////////////////////////////////////////////////
//
//  Common algorithms generic implementaions
//
////////////////////////////////////////////////////////////////////////////////
//==============================================================================
//  Find_index
//==============================================================================
// Implements a binary search procedure
template<typename ElementType, typename ValueType>
uint Find_index(const Array<ElementType> & arr, const ValueType & val)
{
    // Zero size specialization
    if (!arr.Get_size()) { return -1; }

    uint end = arr.Get_size() - 1;

    // Out of bounds specialization
    if (arr[end] <= val) { return end; }

    // General case specialization
    uint start = 0;
    uint ind = end >> 1;

    while (start < ind)
    {
        if (arr[ind] > val) { end = ind; }
        else { start = ind; }

        ind = start + ((end - start) >> 1);
    }

    return ind;
}

//==============================================================================
//
//==============================================================================
/// TODO /// Move sort algorithms here

////////////////////////////////////////////////////////////////////////////////
//
//  Tuple
//
////////////////////////////////////////////////////////////////////////////////
//==============================================================================
//  Constructor
//==============================================================================
Tuple::Tuple(const Tuple & tup, Type (* const Mut)(const Type &))
{
    Set(tup, Mut);
}

//==============================================================================
//  Subtraction
//==============================================================================
// Calculates the difference between the left argument
// and a Tuple which entries all equal to the right argument
template<uint size, typename Type>
Tuple operator-(Tuple<size, Type> left, const Type & right)
{
    left -= right;
    return left;
}

template<uint size, typename Type>
Tuple operator-(Tuple<size, Type> left, const Tuple<size, Type> & right)
{
    left -= right;
    return left;
}

//==============================================================================
//  Addition
//==============================================================================
template<uint size, typename Type>
Tuple operator+(Tuple<size, Type> left, const Type & right)
{
    left += right;
    return left;
}

template<uint size, typename Type>
Tuple operator+(Tuple<size, Type> left, const Tuple<size, Type> & right)
{
    left += right;
    return left;
}

/// useful /// //==============================================================================
/// useful /// //  Sum
/// useful /// //==============================================================================
/// useful /// // Computes the sum of the first "len" entries of a tuple
/// useful /// template<uint size, typename Type>
/// useful /// Type Sum(const Tuple<size, Type> & tup, const uint len)
/// useful /// {
/// useful ///     Type res = 0;
/// useful /// 
/// useful ///     for (int s = 0; s < std::min(size, len); ++s) { res += tup.data[s]; }
/// useful /// 
/// useful ///     return res;
/// useful /// }

//==============================================================================
//  Product
//==============================================================================
template<uint size, typename Type>
Type Product(const Tuple<size, Type> & tup, const uint count)
{
    Type res = 1;

    for (int s = 0; s < std::min(size, count); ++s) { res *= tup.data[s]; }

    return res;
}

//==============================================================================
//  Initialization
//==============================================================================
template<uint size, typename Type>
void Tuple::Set(const Type * dat)
{
    if (dat) { Copy(size, dat, 1, data, 1); }
    else { Set(0); }
}

template<uint size, typename Type = uint>
void Tuple::Set(const Tuple & tup, Type (* const Mut)(const Type &))
{
    for (int s = 0; s < size; ++s) { data[s] = Mut(tup.data[s]); }
}

//==============================================================================
//  Compound assignment by difference
//==============================================================================
// Subtracts a Tuple which entries all equal to the right argument
// from the left argument
Tuple & Tuple::operator-=(const Type & val)
{
     Axpy(size, -1, &val, 0, data, 1);
     return *this;
}

Tuple & Tuple::operator-=(const Tuple & tup)
{
     Axpy(size, -1, tup.data, 1, data, 1);
     return *this;
}

//==============================================================================
//  Compound assignment by sum
//==============================================================================
Tuple & Tuple::operator+=(const Type & val)
{
     Axpy(size, 1, &val, 0, data, 1);
     return *this;
}

Tuple & Tuple::operator+=(const Tuple & tup)
{
     Axpy(size, 1, tup.data, 1, data, 1);
     return *this;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Array
//
////////////////////////////////////////////////////////////////////////////////
//==============================================================================
//  Data management
//==============================================================================
template<typename Type>
void Array::Reallocate(const uint cap)
{
    if (capacity != cap)
    {
        capacity = cap;

        if (!capacity) { free(data); data = NULL; }
        else
        {
            Type * dat = (Type *)realloc(data, size * sizeof(Type))

            if (dat) { data = dat; }
            else
            {
                /// TODO /// Error handling
                exit(EXIT_FAILURE);
            }
        }

        if (capacity < size) { size = capacity; }
    }
}

//==============================================================================
//  Initialization
//==============================================================================
template<typename Type>
void Array::Set(const uint cap, const Type & val)
{
    Reallocate(cap);
    size = capacity;
    Copy(size, &val, 0, data, 1);
}

template<typename Type>
void Array::Set(const uint cap, const Type * dat)
{
    Reallocate(cap);
    size = capacity;
    if (dat) { Copy(size, dat, 1, data, 1); }
}

template<typename Type>
void Array::Set(const Array & arr)
{
    Reallocate(arr.size);
    size = capacity;
    Copy(size, arr.data, 1, data, 1);
}

//==============================================================================
//  Constructors
//==============================================================================
template<typename Type>
Array::Array(const uint cap, const Type & val): capacity(0), size(0), data(NULL)
{
    Set(cap, val);
}

template<typename Type>
Array::Array(const uint cap, const Type * dat): capacity(0), size(0), data(NULL)
{
    Set(cap, dat);
}

//==============================================================================
//  Append
//==============================================================================
template<typename Type>
void Array::Append()
{
    if (capacity == size) { Reallocate(2 * size); }

    ++size;
}

template<typename Type>
void Array::Append(const Type & val)
{
    if (capacity == size) { Reallocate(2 * size); }

    data[size] = val;
    ++size;
}

//============================================================================//
//  Identify_ghost_index
//============================================================================//
template<Dim dim>
uint8_t Identify_ghost_index(const Tuple<dim, Dir> & directions)
{
    uint8_t res = directions[0];

    for (int d = 1; d < dim; ++d)
    {
        res *= 3;
        res += directions[d];
    }

    return res;
}

//============================================================================//
//  Identify_ghost_directions
//============================================================================//
template<Dim dim>
Tuple<dim, Dir> & Identify_ghost_directions(const uint8_t ind)
{
    Tuple<dim, Dir> res;

    uint8_t tmp = ind;

    for (int d = 0; d < dim - 1; ++d)
    {
        res[d] = tmp % 3;
        tmp /= 3;
    }

    res[dim - 1] = tmp;

    return res;
}

//============================================================================//
//  Reflect_ghost_directions
//============================================================================//
template<Dim dim>
Tuple<dim, Dir> & Reflect_ghost_directions(Tuple<dim, Dir> & directions)
{
    for (int d = 0; d < dim; ++d)
    {
        if (directions[d] == LEFT) { directions[d] = RIGHT; }
        else if (directions[d] == RIGHT) { directions[d] = LEFT; }
    }

    return directions;
}

} // namespace maxwell

#endif // UTILITY_HPP
