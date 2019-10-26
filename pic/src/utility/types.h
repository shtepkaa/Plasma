#ifndef TYPES_H
#define TYPES_H

namespace maxwell {

typedef unsigned int uint;

// Utility constants
enum Util: uint { UNDEFINED = ~0 };

// Dimensionality
enum Dim: uint { ONE_DIM = 1, TWO_DIM, THREE_DIM };

// Patch order
enum Order: uint { CARTESIAN = 0, HILBERTIAN };

// Direction
enum Dir: uint8_t { LEFT = 0, CENTER, RIGHT };


} // namespace maxwell

#endif // TYPES_H
