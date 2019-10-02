#ifndef UTILITY_H
#define UTILITY_H

namespace maxwell {

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

//============================================================================//
//  DeviceArray
//============================================================================//
// Contains device array of given type
template<typename Type>
class DeviceArray
{
    private:

        uint size;

        Type * data;

    public:

        // data management
        void Reallocate(const uint = 0);

        // set
        void Copy(const uint, const Type &);
        __device__ void Copy(const uint, const Type * = NULL);
        void Copy(const Array &);

        // initialization
        DeviceArray(): size(0), data(NULL) {}
        DeviceArray(const uint, const Type &);
        DeviceArray(const uint, const Type * = NULL);
        DeviceArray(const Array & arr): size(0), data(NULL) { Copy(arr); }

        Array & operator=(const Array & arr) { Copy(arr); return *this; }

        // deallocation
        ~Array() { Reallocate(); }

        // field access
        inline uint Get_size() const { return size; }

        inline Type * Export_data() const { return data; }
};

} // namespace maxwell

#endif // UTILITY_H
