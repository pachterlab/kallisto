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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Logic types and operations for metaprogramming.
//
// We can perform logical computations both type and value based.  The type
// based metafunctions are the "main" ones, the value based meta functions
// have the postfix "C" (e.g. "Or" vs "OrC"), similar to how Boost's MPL does
// this.
// ==========================================================================

// TODO(holtgrew): Group all these metafunctions?
// TODO(holtgrew): Make a comment on the C prefix and the auto-evaluation to TParam::Type in the group?

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_METAPROGRAMMING_LOGIC_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_METAPROGRAMMING_LOGIC_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tags True and False
// ----------------------------------------------------------------------------

/*!
 * @defgroup LogicalValuesTags Logical Values
 * @brief Tags for representing true and false.
 *
 * @section Examples
 *
 * Print the values of the tags/metafunctions <tt>True</tt> and <tt>False</tt>.
 *
 * @snippet demos/basic/metaprogramming_logic.cpp tags true false
 *
 * @section Inheriting from True and False
 *
 * The two tags True and False have the special property that they can also be used as metafunctions and both have a
 * <tt>VALUE</tt> as well as a <tt>TYPE</tt>.  This property makes it very convenient to define metafunctions by inheriting from the <tt>True</tt> or <tt>False</tt>.
 *
 * @snippet demos/basic/metaprogramming_logic.cpp inheriting from true false
 *
 * The metafunction <tt>IsInt32</tt> can now be used as follows.
 *
 * @snippet demos/basic/metaprogramming_logic.cpp using isint32
 */

/*!
 * @tag LogicalValuesTags#True
 * @headerfile <seqan/basic.h>
 * @brief Representation for True.
 *
 * @signature struct True;
 * @signature True::Type;
 * @signature bool True::VALUE = true;
 */

/*!
 * @tag LogicalValuesTags#False
 * @headerfile <seqan/basic.h>
 * @brief Representation for False.
 *
 * @signature struct False;
 * @signature False::Type;
 * @signature bool False::VALUE = false;
 */

struct True
{
    typedef True Type;
    static const bool VALUE = true;

    operator bool () const
    {
        return true;
    }

    template <typename TValue>
    bool operator() (TValue const &) const
    {
        return true;
    }
};

struct False
{
    typedef False Type;
    static const bool VALUE = false;

    operator bool () const
    {
        return false;
    }

    template <typename TValue>
    bool operator() (TValue const &) const
    {
        return false;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Eval
// ----------------------------------------------------------------------------

/*!
 * @defgroup LogicMetaprogramming Logic Metaprogramming
 * @brief Logic Metaprogramming operations.
 *
 * This group contains metafunctions for logical operations.
 *
 * For each boolean operation, there is one metafunction called <tt>Operation</tt> and one called <tt>OperationC</tt>
 * (i.e. having a suffix <tt>C</tt>.  The first one works on boolean tag types <tt>True</tt> and <tt>False</tt>.  The second one takes <tt>bool</tt> constant parameters.
 *
 * The metafunctions allow a shortcut using the SFNAE (substitution failure is not an error) feature of the C++
 * programming language.  When passing metafunctions returning <tt>True</tt> and <tt>False</tt> in their <tt>Type</tt>
 * member typedef, you can ommit the <tt>::Type</tt>.
 *
 * Here is an example for this:
 *
 * @snippet demos/basic/metaprogramming_logic.cpp shortcut to type feature
 *
 * @see LogicalValuesTags#True
 * @see LogicalValuesTags#False
 */

/*!
 * @mfn LogicMetaprogramming#Eval
 * @headerfile <seqan/basic.h>
 * @brief Conversion from <tt>bool</tt> to tags <tt>True</tt> and <tt>False</tt>.
 *
 * @signature Eval<VALUE>::Type
 *
 * @tparam VALUE A <tt>bool</tt> to convert to <tt>True</tt>/<tt>False</tt>.
 *
 * @return Type The resulting tag, one of <tt>True</tt> and <tt>False</tt>.
 *
 * @section Examples
 *
 * We define the following two helper functions.
 *
 * @snippet demos/basic/metaprogramming_logic.cpp true false print helpers
 *
 * Now, we can write the following code and achieve the following output:
 *
 * @snippet demos/basic/metaprogramming_logic.cpp print bool type eval
 */

template <bool B>
struct Eval : False {};

template <>
struct Eval<true> : True {};

// ----------------------------------------------------------------------------
// Metafunction Not
// ----------------------------------------------------------------------------

/*!
 * @mfn LogicMetaprogramming#Not
 * @headerfile <seqan/basic.h>
 * @brief Metaprograming boolean "not" operator.
 *
 * @signature Not<TBool>::Type;
 * @signature Not<TBool>::VALUE;
 *
 * @tparam TBool The tag to negate.
 *
 * @return Type  The inverted TBool.
 * @return VALUE Shortcut for <tt>Not&lt;TBool&gt;::Type::VALUE</tt>.
 *
 * @section Example
 *
 * We define the following two helper functions.
 *
 * @snippet demos/basic/metaprogramming_logic.cpp true false print helpers
 *
 * Now, we can write the following code using Not.
 *
 * @snippet demos/basic/metaprogramming_logic.cpp print bool type not
 */

template <typename TBool1>
struct Not : Not<typename TBool1::Type> {};

template <>
struct Not<False> : True {};
template <>
struct Not<True> : False {};

// ----------------------------------------------------------------------------
// Metafunction NotC
// ----------------------------------------------------------------------------

/*!
 * @mfn LogicMetaprogramming#NotC
 * @headerfile <seqan/basic.h>
 * @brief Metaprograming boolean "not" operator for values.
 *
 * @signature NotC<BOOL>::Type;
 * @signature NotC<BOOL>::VALUE;
 *
 * @tparam BOOL The <tt>bool</tt> value to negate.
 *
 * @return Type  The corresponding Tag for <tt>!BOOL</tt> (<tt>True</tt>/<tt>False</tt>).
 * @return VALUE Shortcut for <tt>NotC&lt;BOOL&gt;::Type::VALUE</tt>.
 *
 * @section Example
 *
 * We define the following two helper functions.
 *
 * @snippet demos/basic/metaprogramming_logic.cpp true false print helpers
 *
 * Now, we can write the following code using NotC.
 *
 * @snippet demos/basic/metaprogramming_logic.cpp print bool type notc
 */

template <bool B>
struct NotC;

template <>
struct NotC<false> : True {};
template <>
struct NotC<true> : False {};

// ----------------------------------------------------------------------------
// Metafunction Or
// ----------------------------------------------------------------------------

/*!
 * @mfn LogicMetaprogramming#Or
 * @headerfile <seqan/basic.h>
 * @brief Metaprograming "or" operator.
 *
 * @signature Or<TLhs, TRhs>::Type;
 * @signature Or<TLhs, TRhs>::VALUE;
 *
 * @tparam TLhs The left hand side logical value tag.
 * @tparam TRhs The right hand side logical value tag.
 *
 * @return Type  The logical value tag result for the or operation.
 * @return VALUE Shortcut for <tt>Or<TLhs, TRhs>::Type::VALUE</tt>.
 *
 * @section Example
 *
 * We define the following two helper functions.
 *
 * @snippet demos/basic/metaprogramming_logic.cpp true false print helpers
 *
 * Now, we can write the following code using Or.
 *
 * @snippet demos/basic/metaprogramming_logic.cpp print bool type or
 */

template <typename TBool1, typename TBool2>
struct Or : Or<typename TBool1::Type, typename TBool2::Type> {};

template <>
struct Or<False, False> : False {};
template <>
struct Or<False, True> : True {};
template <>
struct Or<True, False> : True {};
template <>
struct Or<True, True> : True {};

// ----------------------------------------------------------------------------
// Metafunction OrC
// ----------------------------------------------------------------------------

/*!
 * @mfn LogicMetaprogramming#OrC
 * @headerfile <seqan/basic.h>
 * @brief Metaprograming boolean "or" operator.
 *
 * @signature OrC<LHS, RHS>::Type;
 * @signature OrC<LHS, RHS>::VALUE;
 *
 * @tparam LHS Left hand side <tt>bool</tt> constant.
 * @tparam RHS Right hand side <tt>bool</tt> constant.
 *
 * @return Type  The logical value tag result for the or operation.
 * @return VALUE Shortcut for <tt>OrC<LHS, RHS>::Type::VALUE</tt>.
 *
 * @section Example
 *
 * We define the following two helper functions.
 *
 * @snippet demos/basic/metaprogramming_logic.cpp true false print helpers
 *
 * Now, we can write the following code using OrC.
 *
 * @snippet demos/basic/metaprogramming_logic.cpp print bool type orc
 */

template <bool B1, bool B2>
struct OrC;

template <>
struct OrC<false, false> : False {};
template <>
struct OrC<false, true> : True {};
template <>
struct OrC<true, false> : True {};
template <>
struct OrC<true, true> : True {};

// ----------------------------------------------------------------------------
// Metafunction And
// ----------------------------------------------------------------------------

/*!
 * @mfn LogicMetaprogramming#And
 * @headerfile <seqan/basic.h>
 * @brief Metaprograming "and" operatand.
 *
 * @signature And<TLhs, TRhs>::Type;
 * @signature And<TLhs, TRhs>::VALUE;
 *
 * @tparam TLhs The left hand side logical value tag.
 * @tparam TRhs The right hand side logical value tag.
 *
 * @return Type  The logical value tag result fand the and operation.
 * @return VALUE Shandtcut fand <tt>And<TLhs, TRhs>::Type::VALUE</tt>.
 *
 * @section Example
 *
 * We define the following two helper functions.
 *
 * @snippet demos/basic/metaprogramming_logic.cpp true false print helpers
 *
 * Now, we can write the following code using And.
 *
 * @snippet demos/basic/metaprogramming_logic.cpp print bool type and
 */

template <typename TBool1, typename TBool2>
struct And : And<typename TBool1::Type, typename TBool2::Type> {};

template <>
struct And<False, False> : False {};
template <>
struct And<False, True> : False {};
template <>
struct And<True, False> : False {};
template <>
struct And<True, True> : True {};

// ----------------------------------------------------------------------------
// Metafunction AndC
// ----------------------------------------------------------------------------

/*!
 * @mfn LogicMetaprogramming#AndC
 * @headerfile <seqan/basic.h>
 * @brief Metaprograming boolean "and" operator.
 *
 * @signature AndC<LHS, RHS>::Type;
 * @signature AndC<LHS, RHS>::VALUE;
 *
 * @tparam LHS Left hand side <tt>bool</tt> constant.
 * @tparam RHS Right hand side <tt>bool</tt> constant.
 *
 * @return Type  The logical value tag result for the or operation.
 * @return VALUE Shortcut for <tt>AndC<LHS, RHS>::Type::VALUE</tt>.
 *
 * @section Example
 *
 * We define the following two helper functions.
 *
 * @snippet demos/basic/metaprogramming_logic.cpp true false print helpers
 *
 * Now, we can write the following code using AndC.
 *
 * @snippet demos/basic/metaprogramming_logic.cpp print bool type andc
 */

template <bool B1, bool B2>
struct AndC;

template <>
struct AndC<false, false> : False {};
template <>
struct AndC<false, true> : False {};
template <>
struct AndC<true, false> : False {};
template <>
struct AndC<true, true> : True {};

// ----------------------------------------------------------------------------
// Metafunction If
// ----------------------------------------------------------------------------

/*!
 * @mfn LogicMetaprogramming#If
 * @headerfile <seqan/basic.h>
 * @brief Metaprogramming implication.
 *
 * @signature If<TCondition, TResultTrue, TResultFalse>::Type
 *
 * @tparam TCondition  The condition.
 * @tparam TResultTrue Result if TCondition evaluates to <tt>True</tt>.
 * @tparam TResultFalse Result if TCondition evaluates to <tt>False</tt>.
 *
 * @return Type The resulting type.
 *
 * @section Example
 *
 * We define the following two helper functions.
 *
 * @snippet demos/basic/metaprogramming_logic.cpp true false print helpers
 *
 * Now, we can write the following code using If.
 *
 * @snippet demos/basic/metaprogramming_logic.cpp print bool type if
 */

template <typename TCondition, typename T1, typename T2>
struct If : If<typename TCondition::Type, T1, T2>{};

template <typename T1, typename T2>
struct If<True, T1, T2>
{
    typedef T1 Type;
};

template <typename T1, typename T2>
struct If<False, T1, T2>
{
    typedef T2 Type;
};



template <typename TArgT, typename TArgF>
inline TArgT SEQAN_FORWARD_ARG
ifSwitch(True, TArgT SEQAN_FORWARD_ARG argTrue, TArgF const &)
{
    return SEQAN_FORWARD(TArgT, argTrue);
}

template <typename TArgT, typename TArgF>
inline TArgF SEQAN_FORWARD_ARG
ifSwitch(False, TArgT const &, TArgF SEQAN_FORWARD_ARG argFalse)
{
    return SEQAN_FORWARD(TArgF, argFalse);
}

#ifndef SEQAN_CXX11_STANDARD

template <typename TArgT, typename TArgF>
inline TArgT const &
ifSwitch(True, TArgT const & argTrue, TArgF const &)
{
    return SEQAN_FORWARD(TArgT, argTrue);
}

template <typename TArgT, typename TArgF>
inline TArgF const &
ifSwitch(False, TArgT const &, TArgF const & argFalse)
{
    return SEQAN_FORWARD(TArgF, argFalse);
}

#endif

// ----------------------------------------------------------------------------
// Metafunction IfC
// ----------------------------------------------------------------------------

/*!
 * @mfn LogicMetaprogramming#IfC
 * @headerfile <seqan/basic.h>
 * @brief Metaprogramming boolean, implication.
 *
 * @signature If<CONDITION, TResultTrue, TResultFalse>::Type
 *
 * @tparam CONDITION    The condition, <tt>bool</tt>.
 * @tparam TResultTrue  Result if TCondition evaluates to <tt>True</tt>.
 * @tparam TResultFalse Result if TCondition evaluates to <tt>False</tt>.
 *
 * @return Type The resulting type.
 *
 * @section Example
 *
 * We define the following two helper functions.
 *
 * @snippet demos/basic/metaprogramming_logic.cpp true false print helpers
 *
 * Now, we can write the following code using If.
 *
 * @snippet demos/basic/metaprogramming_logic.cpp print bool type if
 */

template <bool FLAG, typename T1, typename T2>
struct IfC
{
    typedef T1 Type;
};

template <typename T1, typename T2>
struct IfC<false, T1, T2>
{
    typedef T2 Type;
};

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BASIC_METAPROGRAMMING_LOGIC_H_
