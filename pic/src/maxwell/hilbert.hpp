#ifndef HILBERT_HPP
#define HILBERT_HPP

#include <iostream>

namespace maxwell {

/// FIXME ///

template<Dim d>
uint Hilbert_index(const uint, const uint *, const uint, const uint)
{
    if (d == 1D)
    {
        std::cerr << "Error: Hilbert index is not implemented in 1D"
            << std::endl;

        exit(EXIT_FAILURE);
    }
}

template<Dim d>
void Hilbert_index_inverse(
    const uint, uint *, const uint, const uint, const uint
)
{
    if (d == 1D)
    {
        std::cerr << "Error: Hilbert index inverse is not implemented in 1D"
            << std::endl;

        exit(EXIT_FAILURE);
    }
}

//============================================================================//
//  General Hilbert index
//============================================================================//
template<Dim d>
uint General_Hilbert_index(const uint *, const int *)
{
    if (d == 1D)
    {
        std::cerr << "Error: General Hilbert index is not implemented in 1D"
            << std::endl;

        exit(EXIT_FAILURE);
    }
}

template<Dim d>
void General_Hilbert_index_inverse(const uint *, uint *, const uint)
{
    if (d == 1D)
    {
        std::cerr << "Error: General Hilbert index inverse is not implemented "
            << "in 1D" << std::endl;

        exit(EXIT_FAILURE);
    }
}

} // namespace maxwell

#endif // HILBERT_HPP
