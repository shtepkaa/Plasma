#ifndef UTILITY_H
#define UTILITY_H

#include "types.h"

namespace maxwell {

/*******************************************************************************
*
*   BLAS-like generic functions
*
*******************************************************************************/
//============================================================================//
//  Copy
//============================================================================//
// Implements BLAS-like Copy,
// works correctly for copy-constructible types
template<typename Type>
void Copy(const uint, const Type *, const uint, Type *, const uint);

//============================================================================//
//  Axpy
//============================================================================//
// Implements BLAS-like Axpy,
// works correctly for copy-constructible types,
// allowing compound assignment and multiplication operators
template<typename Type>
void Axpy(
    const uint, const Type &, const Type *, const uint, Type *, const uint
);

/*******************************************************************************
*
*   Sub-hypercube counting functions
*
*******************************************************************************/
//============================================================================//
//  Power
//============================================================================//
// Calculates integer power for integer
uint Power(const uint, const uint);

//============================================================================//
//  Sub_hypercube_count
//============================================================================//
/// FIXME /// probably switch to c++11 constexpr function
// Counts sub-hypercubes of given dimension in cube 
uint Sub_hypercube_count(const uint, const uint);

//============================================================================//
//  Total_hypercube_count
//============================================================================//
/// FIXME /// probably switch to c++11 constexpr function
// Counts total amount of sub-hypercubes in cube
uint Total_hypercube_count(const uint dim) { return Power(3U, dim); }

/*******************************************************************************
*
*   Tuple
*
*******************************************************************************/
// Contains tuple of given type 
// works correctly for Copy-constructible types
// allowing conversion from 0 and 1
template<uint size, typename Type = uint>
struct Tuple
{
    friend Tuple operator+(Tuple, const Tuple &);
    friend Type Product(const Tuple &);

    Type data[size];

    // initialization
    void Set(const Type & val) { Copy(size, &val, 0, data, 1); }
    void Set(const Type * = NULL);
    void Set(const Tuple & tup) { Copy(size, tup.data, 1, data, 1); }

    Tuple(const Type & val) { Set(val); }  
    Tuple(const Type * dat = NULL) { Set(dat); }
    Tuple(const Tuple & tup) { Set(tup); }

    Tuple & operator=(const Tuple & tup) { Set(tup); return *this; }

    // conversion
    inline operator const Type *() const { return data; } 

    // compound addition
    Tuple & operator+=(const Type &);
    Tuple & operator+=(const Tuple &);

    // element mutate / access
    inline Type & operator[](const uint ind) { return data[ind]; }
    inline const Type & operator[](const uint ind) const { return data[ind]; }
};

/*******************************************************************************
*
*   Array
*
*******************************************************************************/
// Contains host array of given type
template<typename Type>
class Array
{
    private:

        uint size;

        Type * data;

    public:

        // data management
        void Reallocate(const uint = 0);

        // initialization
        void Set(const uint, const Type &);
        void Set(const uint, const Type * = NULL);
        void Set(const Array &);

        Array(): size(0), data(NULL) {}
        Array(const uint, const Type &);
        Array(const uint, const Type * = NULL);
        Array(const Array & arr): size(0), data(NULL) { Set(arr); }

        Array & operator=(const Array & arr) { Set(arr); return *this; }

        // deallocation
        ~Array() { Reallocate(); }

        // field access
        inline uint Get_size() const { return size; }

        // element mutate / access
        inline Type & operator[](const uint in) { return data[in]; }
        inline const Type & operator[](const uint in) const { return data[in]; }

        // data access
        inline const Type * Export_data() const { return data; }
};

} // namespace maxwell

#endif // UTILITY_H
