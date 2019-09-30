// src/utility/utility.cpp

#include <cmath>

namespace maxwell {

//============================================================================//
//  Sub_hypercube_count
//============================================================================//
// recursive
uint Sub_hypercube_count(const int sub_dim, const int dim)
{
    if (sub_dim < 0 || dim < 0 || sub_dim >= dim) { return 0U; }
    else if (sub_dim == dim == 0) { return 1U; }
    else
    {
        return (Sub_hypercube_count(sub_dim, dim - 1) << 1)
            + Sub_hypercube_count(sub_dim - 1, dim - 1);
    }
}

//============================================================================//
//  Total_hypercube_count
//============================================================================//
uint Total_hypercube_count(const uint dim)
{
    return std::ceil(std::pow(3., dim));
}

} // namespace maxwell

// src/utility/utility.cpp
