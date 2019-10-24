#ifndef TYPES_H
#define TYPES_H

namespace maxwell {

typedef unsigned int uint;

// Dimensionality
enum Dim: uint { OneDim = 1, TwoDim, ThreeDim };

// Patch order
enum Order: uint { Cartesian = 0, Hilbertian };

// Direction
enum Dir: uint8_t { Left = 0, Center, Right };


} // namespace maxwell

#endif // TYPES_H
