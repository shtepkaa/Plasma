#ifndef UTILITY_H
#define UTILITY_H

namespace maxwell {

template<typename Type>
void Init_array(const uint, const Type &, Type *&);

template<typename Type>
void Init_array(const uint, const Type *, Type *&);

} // namespace maxwell

#endif // UTILITY_H
