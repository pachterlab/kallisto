// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Facade header for the basic/smart_pointer sub module.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_BASIC_SMART_POINTER_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_BASIC_SMART_POINTER_H_

// --------------------------------------------------------------------------
// Dependencies
// --------------------------------------------------------------------------

#include <seqan/platform.h>
#include <seqan/basic/basic_fundamental.h>
#include <seqan/basic/basic_allocator.h>  // TODO(holtgrew): Drop dependency, only for pointer holder which does not work anyway.
#include <seqan/basic/basic_type.h>  // TODO(holtgrew): Really required? Deprecated header?

// --------------------------------------------------------------------------
// Sub Module Headers
// --------------------------------------------------------------------------

#include <seqan/basic/holder_base.h>
#include <seqan/basic/holder_simple.h>
#include <seqan/basic/holder_tristate.h>

#include <seqan/basic/volatile_ptr.h>

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BASIC_BASIC_SMART_POINTER_H_
