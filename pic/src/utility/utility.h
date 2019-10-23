#ifndef UTILITY_H
#define UTILITY_H

#include "types.h"

namespace maxwell {

////////////////////////////////////////////////////////////////////////////////
//
//  BLAS-like generic functions
//
////////////////////////////////////////////////////////////////////////////////
//==============================================================================
//  Copy
//==============================================================================
// Implements BLAS-like array copy,
// works correctly for copy-constructible types
template<typename Type>
void Copy(const uint, const Type *, const uint, Type *, const uint);

//==============================================================================
//  Axpy
//==============================================================================
// Implements BLAS-like Axpy,
// works correctly for copy-constructible types,
// allowing compound assignment and multiplication operators
//==============================================================================
template<typename Type>
void Axpy(
    const uint, const Type &, const Type *, const uint, Type *, const uint
);

////////////////////////////////////////////////////////////////////////////////
//
//  Sub-hypercube counting functions
//
////////////////////////////////////////////////////////////////////////////////
//==============================================================================
//  Binary_logarithm
//==============================================================================
// Calculates exponent for a given integer power of two
uint Binary_logarithm(const uint);

//==============================================================================
//  Power
//==============================================================================
// Calculates integer power for integer base
//==============================================================================
uint Power(const uint, const uint);

//==============================================================================
//  Sub_hypercube_count
//==============================================================================
/// FIXME /// probably switch to c++11 constexpr function
// Counts sub-hypercubes of given dimension in cube
//==============================================================================
uint Sub_hypercube_count(const uint, const uint);

//==============================================================================
//  Total_hypercube_count
//==============================================================================
/// FIXME /// probably switch to c++11 constexpr function
// Counts total amount of sub-hypercubes in cube
// -----------------------------------------------------------------------------
uint Total_hypercube_count(const uint dim) { return Power(3, dim); }

////////////////////////////////////////////////////////////////////////////////
//
//  Common algorithms generic implementaions
//
////////////////////////////////////////////////////////////////////////////////
//==============================================================================
//  Binary search
//==============================================================================
// Identifies the left bound of interval (including left, excluding right),
// which contains a given value
// -----------------------------------------------------------------------------
// Returns uint(-1) if array is empty
// -----------------------------------------------------------------------------
// Input parameters:
//     arr -- the array of interval bounds
//     val -- a reference value
// -----------------------------------------------------------------------------
template<typename Type>
uint Binary_search(const Array<Type> & arr, const Type & val);

////////////////////////////////////////////////////////////////////////////////
//
//  Tuple
//
////////////////////////////////////////////////////////////////////////////////
// Contains tuple of given type
// -----------------------------------------------------------------------------
// works correctly for Copy-constructible types
// allowing conversion from 0 and 1
// -----------------------------------------------------------------------------
// Template parameters:
//     size -- number of components
//     Type -- supported arithmetic type
////////////////////////////////////////////////////////////////////////////////
template<uint size, typename Type = uint>
struct Tuple
{
    //==========================================================================
    //  Friend functions
    //==========================================================================
    // Perform binary vector operations with two tuples
    // one of which may be constructed from a single given component
    friend Tuple operator-(Tuple, const Type &);
    friend Tuple operator-(Tuple, const Tuple &);

    friend Tuple operator+(Tuple, const Type &);
    friend Tuple operator+(Tuple, const Tuple &);

    // Reduce the tuple by given operation
    /// useful /// friend Type Sum(const Tuple &, const uint = size);
    friend Type Product(const Tuple &, const uint = size);

    //==========================================================================
    //  Data
    //==========================================================================
    Type data[size];

    //==========================================================================
    //  Data management
    //==========================================================================
    // Initializes the tuple with predefined components
    void Set(const Type & val) { Copy(size, &val, 0, data, 1); }
    void Set(const Type * = NULL);
    void Set(const Tuple & tup) { Copy(size, tup.data, 1, data, 1); }
    void Set(const Tuple &, Type (* const)(const Type &));

    Tuple(const Type & val) { Set(val); }
    Tuple(const Type * dat = NULL) { Set(dat); }
    Tuple(const Tuple & tup) { Set(tup); }
    Tuple(const Tuple &, Type (* const)(const Type &));

    //==========================================================================
    //  Assignment methods
    //==========================================================================
    Tuple & operator=(const Tuple & tup) { Set(tup); return *this; }

    // Compound assign with another tuple
    // or a tuple, constructed from a single given component
    Tuple & operator-=(const Type &);
    Tuple & operator-=(const Tuple &);

    Tuple & operator+=(const Type &);
    Tuple & operator+=(const Tuple &);

    //==========================================================================
    //  Access / mutate methods
    //==========================================================================
    // Returns a pointer to data
    inline operator Type * const() { return data; }
    inline operator const Type * const() const { return data; }

    // Gets / sets component
    inline Type & operator[](const uint ind) { return data[ind]; }
    inline const Type & operator[](const uint ind) const { return data[ind]; }
};

////////////////////////////////////////////////////////////////////////////////
//
//  Array
//
////////////////////////////////////////////////////////////////////////////////
// Contains host array of given type of variadic size
// -----------------------------------------------------------------------------
// Template parameter:
//     Type -- supported arithmetic type
////////////////////////////////////////////////////////////////////////////////
template<typename Type>
class Array
{
    private:

        // total number of allocated elements
        uint capacity;

        uint size;

        Type * data;

    public:

        //======================================================================
        //  Data management
        //======================================================================
        // Changes capacity
        void Reallocate(const uint = 0);

        // Removes trailing elements from the end
        void Truncate() { Reallocate(size); }

        // Initialize the array with predefined elements
        void Set(const uint, const Type &);
        void Set(const uint, const Type * = NULL);
        void Set(const Array &);

        // Construct the array
        Array(): capacity(0), size(0), data(NULL) {}
        Array(const uint, const Type &);
        Array(const uint, const Type * = NULL);
        Array(const Array & arr): capacity(0), size(0), data(NULL) { Set(arr); }

        Array & operator=(const Array & arr) { Set(arr); return *this; }

        ~Array() { Reallocate(); }

        //======================================================================
        //  Access / mutate methods
        //======================================================================
        // Returns capacity
        inline uint Get_capacity() const { return capacity; }

        // Returns size
        inline uint Get_size() const { return size; }

        // Returns raw pointer to data
        inline operator Type * const() { return data; }
        inline operator const Type * const() const { return data; }

        // Get / set element
        inline Type & operator[](const uint ind) { return data[ind]; }
        inline const Type & operator[](const uint i) const { return data[i]; }

        // Get / set the last element
        inline Type & End() { return data[size - 1]; }
        inline const Type & End() const { return data[size - 1]; }

        // Creates a default element at the end, enlarges the array if needed
        void Append();

        // Copies the given element to the end, enlarges the array if needed
        void Append(const Type &);
};

} // namespace maxwell

#endif // UTILITY_H
