#ifndef UTILITY_HPP
#define UTILITY_HPP

namespace maxwell {

//============================================================================//
//  Copy
//============================================================================//
// Implements BLAS-like copy,
// works correctly for copy-constructible types
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
//  Data management
//============================================================================//
/// ifdef capacity /// template<typename Type>
/// ifdef capacity /// void Array::Reallocate(const uint cap)
/// ifdef capacity /// {
/// ifdef capacity ///     if (cap != capacity)
/// ifdef capacity ///     {
/// ifdef capacity ///         capacity = cap;
/// ifdef capacity ///         data = (Type *)realloc(data, capacity * sizeof(Type));
/// ifdef capacity ///     }
/// ifdef capacity /// }
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
void Array::Copy(const uint siz, const Type & val)
{
    Reallocate(siz);
    copy(size, &val, 0, data, 1);
}

template<typename Type>
void Array::Copy(const uint siz, const Type * dat)
{
    Reallocate(siz);
    if (dat) { copy(size, dat, 1, data, 1); }
}

template<typename Type>
void Array::Copy(const Array & arr)
{
    Reallocate(arr.size);
    copy(size, arr.data, 1, data, 1);
}

//============================================================================//
//  Initialization
//============================================================================//
template<typename Type>
Array::Array(const uint siz, const Type & val): size(0), data(NULL)
{
    Copy(siz, val);
}

template<typename Type>
Array::Array(const uint siz, const Type * dat): size(0), data(NULL)
{
    if (dat) { Copy(siz, dat); }
    else { Reallocate(siz); }
}

} // namespace maxwell

#endif // UTILITY_HPP
