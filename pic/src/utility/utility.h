#ifndef UTILITY_H
#define UTILITY_H

namespace maxwell {

/// template<typename Type>
/// void Init_array(const uint, const Type &, Type *&);
/// 
/// template<typename Type>
/// void Init_array(const uint, const Type *, Type *&);

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
        uint Get_size() const { return size; }
        /// ifdef capacity /// uint Get_capacity() const { return capacity; }

        // element access
        Type & operator[](const uint ind) { return data[ind]; }
        const Type & operator[](const uint ind) const { return data[ind]; }

        const Type * Export_data() const { return data; }
};

} // namespace maxwell

#endif // UTILITY_H
