#ifndef UTILITY_HPP
#define UTILITY_HPP

namespace maxwell {

////////////////////////////////////////////////////////////////////////////////
//
//  BLAS-like generic functions
//
////////////////////////////////////////////////////////////////////////////////
//============================================================================//
//  Copy
//============================================================================//
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

//============================================================================//
//  Axpy
//============================================================================//
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
//============================================================================//
//  Binary search
//============================================================================//
template<typename Type>
uint Binary_search(const uint size, const Type * arr, const Type & val)
{
    // Zero size specialization
    if (!size) { return -1; }

    uint end = size - 1;

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

/// TODO /// Move sort algorithms here

////////////////////////////////////////////////////////////////////////////////
//
//  Tuple
//
////////////////////////////////////////////////////////////////////////////////
//============================================================================//
//  Constructor
//============================================================================//
Tuple::Tuple(const Tuple & tup, Type (* const Mut)(const Type &))
{
    Set(tup, Mut);
}

//============================================================================//
//  Subtraction
//============================================================================//
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

//============================================================================//
//  Addition
//============================================================================//
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

/// useful /// //============================================================================//
/// useful /// //  Sum
/// useful /// //============================================================================//
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

//============================================================================//
//  Product
//============================================================================//
// Computes the product of first "len" entries of a tuple
template<uint size, typename Type>
Type Product(const Tuple<size, Type> & tup, const uint len)
{
    Type res = 1;

    for (int s = 0; s < std::min(size, len); ++s) { res *= tup.data[s]; }

    return res;
}

//============================================================================//
//  Initialization
//============================================================================//
template<uint size, typename Type>
void Tuple::Set(const Type * dat) { data[0] = 0; Copy(size, dat, 0, data, 1); }

// Constructs a new Tuple
// by the entries of the original Tuple undergone the mutator function
template<uint size, typename Type = uint>
void Tuple::Set(const Tuple & tup, Type (* const Mut)(const Type &))
{
    for (int s = 0; s < size; ++s) { data[s] = Mut(tup.data[s]); }
}

//============================================================================//
//  Compound assignment by difference
//============================================================================//
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

//============================================================================//
//  Compound assignment by sum
//============================================================================//
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
//============================================================================//
//  Data management
//============================================================================//
template<typename Type>
void Array::Reallocate(const uint siz)
{
    if (size != siz)
    {
        size = siz;

        if (!siz) { free(data); data = NULL; }
        else { data = (Type *)realloc(data, size * sizeof(Type)); }
    }
}

//============================================================================//
//  Initialization
//============================================================================//
template<typename Type>
void Array::Set(const uint siz, const Type & val)
{
    Reallocate(siz);
    Copy(size, &val, 0, data, 1);
}

template<typename Type>
void Array::Set(const uint siz, const Type * dat)
{
    Reallocate(siz);
    if (dat) { Copy(size, dat, 1, data, 1); }
}

template<typename Type>
void Array::Set(const Array & arr)
{
    Reallocate(arr.size);
    Copy(size, arr.data, 1, data, 1);
}

//============================================================================//
//  Constructors
//============================================================================//
template<typename Type>
Array::Array(const uint siz, const Type & val): size(0), data(NULL)
{
    Set(siz, val);
}

template<typename Type>
Array::Array(const uint siz, const Type * dat): size(0), data(NULL)
{
    if (dat) { Set(siz, dat); }
    else { Reallocate(siz); }
}

} // namespace maxwell

#endif // UTILITY_HPP
