#ifndef HILBERT_HPP
#define HILBERT_HPP

#include <iostream>

namespace maxwell {

template<Dim d>
uint Hilbert_index(const uint, const uint *, const uint, const uint)
{
    if (d == 1D)
    {
        std::cerr << "Error: Hilbert index is not implemented for 1D"
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
        std::cerr << "Error: Hilbert index inverse is not implemented for 1D"
            << std::endl;

        exit(EXIT_FAILURE);
    }
}

//============================================================================//
//  General Hilbert index
//============================================================================//
template<Dim d>
uint General_hilbert_index(const uint *, const int *)
{
    if (d == 1D)
    {
        std::cerr << "Error: Hilbert index is not implemented for 1D"
            << std::endl;

        exit(EXIT_FAILURE);
    }
}

template<Dim d>
void General_hilbert_index_inverse(const uint *, uint *, const uint);
} // namespace maxwell

#endif // HILBERT_HPP
