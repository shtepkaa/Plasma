#ifndef UTILITY_H
#define UTILITY_H

namespace maxwell {

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
uint Total_hypercube_count(const uint);

//============================================================================//
//  Array
//============================================================================//
// Container type
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
            /// ifdef capacity /// uint Get_capacity() const { return capacity; }

        // element mutate
        Type & operator[](const uint);

        // element access
        const Type & operator[](const uint) const;

        inline const Type * Export_data() const { return data; }
};

/// //============================================================================//
/// //  DeviceArray
/// //============================================================================//
/// // Container type
/// template<typename Type>
/// class DeviceArray
/// {
///     private:
/// 
///         uint size;
/// 
///         Type * data;
/// 
///     public:
/// 
///         // data management
///         void Reallocate(const uint = 0);
/// 
///         void Copy(const uint, const Type &);
///         void Copy(const uint, const Type * = NULL);
///         void Copy(const Array &);
/// 
///         // initialization
///         Array(): size(0), data(NULL) {}
///         Array(const uint, const Type &);
///         Array(const uint, const Type * = NULL);
///         Array(const Array & arr): size(0), data(NULL) { Copy(arr); }
/// 
///         Array & operator=(const Array & arr) { Copy(arr); return *this; }
/// 
///         // deallocation
///         ~Array() { Reallocate(); }
/// 
///         // field access
///         inline uint Get_size() const { return size; }
/// 
///         // element access
///         inline Type & operator[](const uint i) { return data[i]; }
///         inline const Type & operator[](const uint i) const { return data[i]; }
/// 
///         inline const Type * Export_data() const { return data; }
/// };

} // namespace maxwell

#endif // UTILITY_H
