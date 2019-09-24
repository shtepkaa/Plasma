#ifndef UTILITY_H
#define UTILITY_H

namespace maxwell {

template<typename Type>
void Init_array(const uint, const Type &, Type *&);

template<typename Type>
void Init_array(const uint, const Type *, Type *&);

template<typename Type>
class Array
{
    private:

        uint size;
        uint capacity;

        Type * data;

    public:

        Array();
        Array(const uint);
        Array(const uint, const Type &);
        Array(const uint, const Type *);

        ~Array();

        Type & operator[](const uint ind) { return data[ind]; }
        const Type & operator[](const uint ind) const { return data[ind]; }
}

} // namespace maxwell

#endif // UTILITY_H
