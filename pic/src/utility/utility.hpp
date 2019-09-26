#ifndef UTILITY_HPP
#define UTILITY_HPP

namespace maxwell {

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
Array::Copy(const uint siz, const Type & val)
{
    Reallocate(siz);
    copy(size, &val, 0, data, 1);
}

template<typename Type>
Array::Copy(const uint siz, const Type * dat)
{
    Reallocate(siz);
    copy(size, dat, 1, data, 1);
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
Array::Array(const uint siz, const Type * dat = NULL): size(0), data(NULL)
{
    if (dat) { Copy(siz, dat); }
    else { Reallocate(siz); }
}

} // namespace maxwell

#endif // UTILITY_HPP
