// src/.../domain.cpp

#include "domain.h"

namespace maxwell {

////////////////////////////////////////////////////////////////////////////////
//
//  BufferMarking
//
////////////////////////////////////////////////////////////////////////////////
BufferMarking::BufferMarking(): 
    size(0),
    offset(~0),
    local_patch_index(~0),
    extern_patch_index(~0),
    local_ghost_index(~0),
    extern_ghost_index(~0)
{}

BufferMarking::BufferMarking(
    const uint buffer_size,
    const uint off,
    const uint local_patch_ind,
    const uint extern_patch_ind,
    const uint8_t local_ghost_ind,
    const uint8_t extern_ghost_ind
):
    size(buffer_size),
    offset(off),
    local_patch_index(local_patch_ind),
    extern_patch_index(extern_patch_ind),
    local_ghost_index(local_ghost_ind),
    extern_ghost_index(extern_ghost_ind)
{}

} // namespace maxwell

// src/.../domain.cpp
