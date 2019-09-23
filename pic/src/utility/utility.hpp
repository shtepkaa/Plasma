#ifndef UTILITY_HPP
#define UTILITY_HPP

namespace maxwell {

/// FIXME ///

template<typename Type>
void Init_array(const uint size, const Type & arg, Type *& res)
{
    res = (Type *)realloc(res, size * sizeof(Type));
    copy(size, arg, 0, res, 1);
}

template<typename Type>
void Init_array(const uint, const Type * args, Type *& res)
{
    res = (Type *)realloc(res, size * sizeof(Type));
    copy(size, args, 1, res, 1);
}

} // namespace maxwell

#endif // UTILITY_HPP
