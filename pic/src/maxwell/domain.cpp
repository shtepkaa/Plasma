// src/.../domain.cpp

#include "domain.h"

namespace maxwell {

////////////////////////////////////////////////////////////////////////////////
//
//  BufferMarking
//
////////////////////////////////////////////////////////////////////////////////
BufferMarking::BufferMarking(
    const uint send_ind,
    const uint recv_ind,
    const uint buffer_size,
    const uint off
):
    send_index(send_ind), recv_index(recv_ind), size(buffer_size), offset(off)
{}

} // namespace maxwell

// src/.../domain.cpp
