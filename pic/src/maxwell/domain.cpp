// src/.../domain.cpp

#include "domain.h"

namespace maxwell {

////////////////////////////////////////////////////////////////////////////////
//
//  TransferMarking
//
////////////////////////////////////////////////////////////////////////////////
TransferMarking::TransferMarking(): 
    size(0),
    offset(~0),
    sending_patch_index(~0),
    receiving_patch_index(~0),
    sending_ghost_index(~0),
    receiving_ghost_index(~0)
{}

TransferMarking::TransferMarking(
    const uint buffer_size,
    const uint buffer_offset,
    const uint sending_patch_ind,
    const uint receiving_patch_ind,
    const uint8_t sending_ghost_ind,
    const uint8_t receiving_ghost_ind
):
    size(buffer_size),
    offset(buffer_offset),
    sending_patch_index(sending_patch_ind),
    receiving_patch_index(receiving_patch_ind),
    sending_ghost_index(sending_ghost_ind),
    receiving_ghost_index(receiving_ghost_ind)
{}

} // namespace maxwell

// src/.../domain.cpp
