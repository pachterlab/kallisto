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
// Misc simple functionality built upon atomic primitives.
// ==========================================================================

// SEQAN_NO_GENERATED_FORWARDS: No forwards are generated for this file.

#ifndef SEQAN_PARALLEL_PARALLEL_ATOMIC_MISC_H_
#define SEQAN_PARALLEL_PARALLEL_ATOMIC_MISC_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function atomicMax()
// ----------------------------------------------------------------------------

/*!
 * @fn atomicMax
 * @headerfile <seqan/parallel.h>
 * @brief Lock-free implemenattion of <tt>x = max(x, y)</tt>.
 *
 * @signature void atomicMax(x, y);
 *
 * @param[in,out] x Integer to set to <tt>max(x, y)</tt>.
 * @param[in]     y Other integer.
 *
 * This is equivalent to
 *
 * @code{.cpp}
 * atomic {
 *     x = max(x, y);
 * }
 * @endcode
 */

template <typename T>
inline void
atomicMax(T volatile & x, T y)
{
    T val = x;
    while (val < y) {
        T m = _max(val, y);
        val = atomicCas(x, val, m);
    }
}

// ----------------------------------------------------------------------------
// Function atomicMin()
// ----------------------------------------------------------------------------

/*!
 * @fn atomicMin
 * @headerfile <seqan/parallel.h>
 * @brief Lock-free implemenattion of <tt>x = min(x, y)</tt>.
 *
 * @signature void atomicMin(x, y);
 *
 * @param[in,out] x Integer to set to <tt>min(x, y)</tt>.
 * @param[in]     y Other integer.
 *
 * This is equivalent to
 *
 * @code{.cpp}
 * atomic {
 *     x = min(x, y);
 * }
 * @endcode
 */

template <typename T>
inline void
atomicMin(T volatile & x, T y)
{
    T val = x;
    while (val > y) {
        T m = _min(val, y);
        val = atomicCas(x, val, m);
    }
}

} // namespace seqan

#endif  // #define SEQAN_PARALLEL_PARALLEL_ATOMIC_MISC_H_
