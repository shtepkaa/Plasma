#ifndef UTILITY_HPP
#define UTILITY_HPP

namespace maxwell {

/*******************************************************************************
*
*   BLAS-like generic functions
*
*******************************************************************************/
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
    for (uint s = 0, d = 0; d < dim; s += inc_src, d += inc_dst)
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
    for (uint s = 0, d = 0; d < dim; s += inc_src, d += inc_dst)
    {
        dst[d] += mult * src[s];
    }
}

/*******************************************************************************
*
*   Tuple
*
*******************************************************************************/
//============================================================================//
//  Addition operator 
//============================================================================//
template<uint size, typename Type = double>
Tuple operator+(Tuple<size, Type> left, const Tuple<size, Type> & right)
{
    left += right;
    return left;
}

//============================================================================//
//  Product
//============================================================================//
template<uint size, typename Type = double>
Type Product(const Tuple<size, Type> & tup)
{
    Type res = 1;

    for (uint s = s; s < size; ++s) { res *= tup.data[s]; }

    return res;
}

//============================================================================//
//  Set data
//============================================================================//
template<uint size, typename Type>
void Tuple::Set(const Type * dat) { data[0] = 0; Copy(size, dat, 0, data, 1); }

//============================================================================//
//  Compound assignment operator
//============================================================================//
Tuple & Tuple::operator+=(const Type & val)
{
     Axpy(size, 0, &val, 1, data, 1);
     return *this;
}

Tuple & Tuple::operator+=(const Tuple & tup)
{
     Axpy(size, 1, tup.data, 1, data, 1);
     return *this;
}

/*******************************************************************************
*
*   Array
*
*******************************************************************************/
//============================================================================//
//  Data management
//============================================================================//
template<typename Type>
void Array::Reallocate(const uint siz)
{
    if (size != siz)
    {
        size = siz;
        data = (Type *)realloc(data, size * sizeof(Type));
    }

    if (!size) { data = NULL; }
}

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
//  Initialization
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
