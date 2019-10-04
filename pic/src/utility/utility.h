#ifndef UTILITY_H
#define UTILITY_H

#include "types.h"

namespace maxwell {

//============================================================================//
//  Copy
//============================================================================//
// Implements BLAS-like Copy,
// works correctly for Copy-constructible types
template<typename Type>
void Copy(const uint, const Type *, const uint, Type *, const uint);

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

//============================================================================//
//  Total_hypercube_count
//============================================================================//
template<Dim dim>
uint & Shift(const uint);

//============================================================================//
//  Tuple
//============================================================================//
// Contains tuple of given type 
template<uint size, typename Type = uint>
struct Tuple
{
    Type data[size];

    // set
    void Set(const Type & val) { Copy(size, &val, 0, data, 1); }
    void Set(const Type * dat = NULL) { data[0] = 0; Copy(size, data, 0, data, 1); }
    void Set(const Tuple & tup) { Copy(); }

    // initialization
    Tuple() { data[0] = 0; Copy(size, data, 0, data, 1); }
    Tuple(const Type & val) { Set(val); }  
    Tuple(const Type * dat = NULL) { Set(dat); }
    Tuple(const Tuple & tup) { Set(tup); }

    Tuple & operator=(const Tuple & tup) { Copy(tup); return *this; }

    // field access
    inline uint Get_size() const { return size; }

    // element mutate
    inline Type & operator[](const uint ind) { return data[ind]; }

    // element access
    inline const Type & operator[](const uint ind) const { return data[ind]; }
};

//============================================================================//
//  Array
//============================================================================//
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

        // set
        void Copy(const uint, const Type &);
        void Copy(const uint, const Type * = NULL);
        void Copy(const Array &);

        // initialization
        Array(): size(0), data(NULL) {}
        Array(const uint, const Type &);
        Array(const uint, const Type * = NULL);
        Array(const Array & arr): size(0), data(NULL) { Copy(arr); }

        Array & operator=(const Array & arr) { Copy(arr); return *this; }

        // deallocation
        ~Array() { Reallocate(); }

        // field access
        inline uint Get_size() const { return size; }

        // element mutate
        inline Type & operator[](const uint i) { return data[i]; }

        // element access
        inline const Type & operator[](const uint i) const { return data[i]; }

        inline const Type * Export_data() const { return data; }
};

} // namespace maxwell

#endif // UTILITY_H
