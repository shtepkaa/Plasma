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
    send_patch_index(~0),
    recv_patch_index(~0),
    send_ghost_index(~0),
    recv_ghost_index(~0)
{}

BufferMarking::BufferMarking(
    const uint buffer_size,
    const uint off,
    const uint send_patch_ind,
    const uint recv_patch_ind,
    const uint8_t send_ghost_ind,
    const uint8_t recv_ghost_ind
):
    size(buffer_size),
    offset(off),
    send_patch_index(send_patch_ind),
    recv_patch_index(recv_patch_ind),
    send_ghost_index(send_ghost_ind),
    recv_ghost_index(recv_ghost_ind)
{}

} // namespace maxwell

// src/.../domain.cpp
