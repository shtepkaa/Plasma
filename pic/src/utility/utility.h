#ifndef UTILITY_H
#define UTILITY_H

#ifdef __CUDACC__
#define HOST_DEVICE  __host__ __device__
#else
#define HOST_DEVICE
#endif

namespace maxwell {

//============================================================================//
//  Array
//============================================================================//
// Container type
template<typename Type, Device>
class Array
{
    private:

        uint size;

        Type * data;

    public:

        // data management
        void Reallocate(const uint = 0);

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

        // element access
        HOST_DEVICE Type & operator[](const uint ind) { return data[ind]; }
        HOST_DEVICE const Type & operator[](const uint ind) const { return data[ind]; }

        inline const Type * Export_data() const { return data; }
};

} // namespace maxwell

#endif // UTILITY_H
