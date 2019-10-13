#ifndef UTILITY_H
#define UTILITY_H

#include "types.h"

namespace maxwell {

////////////////////////////////////////////////////////////////////////////////
//
//  BLAS-like generic functions
//
////////////////////////////////////////////////////////////////////////////////
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

////////////////////////////////////////////////////////////////////////////////
//
//  Sub-hypercube counting functions
//
////////////////////////////////////////////////////////////////////////////////
//============================================================================//
//  Binary_logarithm
//============================================================================//
// Calculates exponent for a given integer power of two
uint Binary_logarithm(const uint);

//============================================================================//
//  Power
//============================================================================//
// Calculates integer power for integer base
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
uint Total_hypercube_count(const uint dim) { return Power(3, dim); }

////////////////////////////////////////////////////////////////////////////////
//
//  Binary search
//
////////////////////////////////////////////////////////////////////////////////
//============================================================================//
//  Binary search
//============================================================================//
// Counts total amount of sub-hypercubes in cube
template<typename Type>
uint Binary_search(const Array<Type> & arr, const Type & val);

////////////////////////////////////////////////////////////////////////////////
//
//  Tuple
//
////////////////////////////////////////////////////////////////////////////////
// Contains tuple of given type 
// works correctly for Copy-constructible types
// allowing conversion from 0 and 1
template<uint size, typename Type = uint>
struct Tuple
{
    //========================================================================//
    //  Friends
    //========================================================================//
    // Subtraction
    friend Tuple operator-(Tuple, const Type &);
    friend Tuple operator-(Tuple, const Tuple &);

    // Addition
    friend Tuple operator+(Tuple, const Type &);
    friend Tuple operator+(Tuple, const Tuple &);

    // Reduction
    /// friend Type Sum(const Tuple &, const uint = size);
    friend Type Product(const Tuple &, const uint = size);

    //========================================================================//
    //  Data
    //========================================================================//
    // Data array
    Type data[size];

    //========================================================================//
    //  Functionality
    //========================================================================//
    // Initialization
    void Set(const Type & val) { Copy(size, &val, 0, data, 1); }
    void Set(const Type * = NULL);
    void Set(const Tuple & tup) { Copy(size, tup.data, 1, data, 1); }
    void Set(const Tuple &, Type (* const)(const Type &));

    // Constructors
    Tuple(const Type & val) { Set(val); }  
    Tuple(const Type * dat = NULL) { Set(dat); }
    Tuple(const Tuple & tup) { Set(tup); }
    Tuple(const Tuple &, Type (* const)(const Type &));

    // Assignment
    Tuple & operator=(const Tuple & tup) { Set(tup); return *this; }

    // Conversion
    inline operator Type * const() { return data; } 
    inline operator Type const * const() const { return data; } 

    // Compound assignment by difference
    Tuple & operator-=(const Type &);
    Tuple & operator-=(const Tuple &);
    
    // Compound assignment by sum
    Tuple & operator+=(const Type &);
    Tuple & operator+=(const Tuple &);

    // Element mutate / access
    inline Type & operator[](const uint ind) { return data[ind]; }
    inline const Type & operator[](const uint ind) const { return data[ind]; }
};

////////////////////////////////////////////////////////////////////////////////
//
//  Array
//
////////////////////////////////////////////////////////////////////////////////
// Contains host array of given type
template<typename Type>
class Array
{
    private:

        uint size;

        Type * data;

    public:

        //====================================================================//
        //  Functionality
        //====================================================================//
        // Data management
        void Reallocate(const uint = 0);

        // Initialization
        void Set(const uint, const Type &);
        void Set(const uint, const Type * = NULL);
        void Set(const Array &);

        // Constructors
        Array(): size(0), data(NULL) {}
        Array(const uint, const Type &);
        Array(const uint, const Type * = NULL);
        Array(const Array & arr): size(0), data(NULL) { Set(arr); }

        // Assignment
        Array & operator=(const Array & arr) { Set(arr); return *this; }

        // Deallocation
        ~Array() { Reallocate(); }

        // Get size
        inline uint Get_size() const { return size; }

        // Element mutate / access
        inline Type & operator[](const uint ind) { return data[ind]; }
        inline const Type & operator[](const uint i) const { return data[i]; }

        // Data access
        inline const Type * Export_data() const { return data; }
};

} // namespace maxwell

#endif // UTILITY_H
