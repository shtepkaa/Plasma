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

//============================================================================//
//  Binary search
//============================================================================//
// 
template<typename Type>
uint Binary_search(const Array<Type> & arr, const Type & val)
{
    uint ind = (arr.Get_size() - 1) >> 1;

    for (uint range = ind; range; )
    {
        if (range > 1) { range = (range + 1) >> 1; }

        if (arr[ind] > val) { ind -= range; }
        else if (arr[ind + 1] < val) { ind += range; }
        else { break; }
    }

    return ind;
}


unsigned Binary_search(
  6     const unsigned size, const unsigned * arr, const unsigned val
  7 )
  8 {
  9     unsigned start = 0;
 10     unsigned end = size - 1;
 11
 12     unsigned len = size - 1;
 13
 14     unsigned ind = len >> 1;
 15
 16     //while (start < end - 1)
 17     while (len > 1)
 18     {
 19         if (arr[ind] > val)
 20         {
 21             //end = ind;
 22             //ind = start + ((end - start) >> 1);
 23             len = ind - start;
 24             ind = start + (len >> 1);
 25         }
 26         else
 27         {
 28             //start = ind;
 29             //ind = start + ((end - start) >> 1);
 30             len = ind - start;
 31             ind = start + ((end - start) >> 1);
 32         }
 33     }
 34
 35     return ind;
 36 }                      

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
void Array::Reallocate(const uint size)
{
    if (size != new_size)
    {
        size = siz;
        data = (Type *)realloc(data, size * sizeof(Type));
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
