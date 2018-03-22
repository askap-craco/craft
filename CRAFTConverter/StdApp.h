///////////////////////////////////////////////////////////////////////////////
//
//  Copyright © 2016, Wayne Arcus, All Rights Reserved.
//
//  Filename:
//
//  Author(s):  Wayne Arcus
//
//  Purpose:    Include file for standard or project specific includes files
//              or definitions that are used frequently, but changed infrequently.
//
//  Notes:
//
//  References:
//
///////////////////////////////////////////////////////////////////////////////

#ifndef STDAPP_H_INCLUDED
#define STDAPP_H_INCLUDED

///////////////////////////////////////////////////////////////////////////////
// Macros.

#ifdef __GNUC__
#  define UNUSED(x) UNUSED_ ## x __attribute__((__unused__))
#else
#  define UNUSED(x) UNUSED_ ## x
#endif

#ifdef _DEBUG
#define UNREFERENCED_PARAMETER(x)       ((void)(x))
#define DBG_UNREFERENCED_PARAMETER(x)   ((void)(x))
#else
#define UNREFERENCED_PARAMETER(x)       ((void)(x))
#define DBG_UNREFERENCED_PARAMETER(x)
#endif

///////////////////////////////////////////////////////////////////////////////
// Headers.

#include <cstdint>                          // For fixed integer types (e.g., uint8_t).
#include <assert.h>                         // Access to assert().
#include <string>                           // Access to std::string.
#include <deque>                            // Access to std::deque.

///////////////////////////////////////////////////////////////////////////////
// Alias and namespaces (or part thereof) used.

using byte_t        = std::uint8_t;          // Common alias for a byte.
using WordDeque_t   = std::deque<uint32_t>;  // Common alias for a byte deque.
using std::string;                           // Using string throughout.

#include "Codec.h"                          // Base-class definition for a codec - keep after
                                            // the above alias for byte_t.

#endif  // STDAPP_H_INCLUDED


