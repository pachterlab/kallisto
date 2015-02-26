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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

// TODO(holtgrew): Split into modified_string_mod_view.h and modified_iterator_mod_view.h.
// TODO(holtgrew): Move out convert()

#ifndef SEQAN_MODIFIER_MODIFIER_VIEW_H_
#define SEQAN_MODIFIER_MODIFIER_VIEW_H_

namespace seqan
{

// ==========================================================================
// Forwards
// ==========================================================================

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class ModView
// --------------------------------------------------------------------------

/*!
 * @class ModViewModifiedIterator
 * @extends ModifiedIterator
 * @headerfile <seqan/modifier.h>
 *
 * @brief Transforms the character of a host using a custom functor.
 *
 * @signature template <typename THost, typename TFunctor>
 *            class ModifiedIterator<THost, ModView<TFunctor> >;
 *
 * @tparam THost    The host iterator.
 * @tparam TFunctor A unary functor type.
 */

/*!
 * @class ModViewModifiedString
 * @extends ModifiedString
 * @headerfile <seqan/modifier.h>
 *
 * @brief Transforms the character of a host using a custom functor.
 *
 * @signature template <typename THost, typename TFunctor>
 *            class ModifiedString<THost, ModView<TFunctor> >;
 *
 * @tparam THost    The host iterator.
 * @tparam TFunctor A unary functor type.
 */

template <typename TFunctor>
struct ModView {};

template <typename TFunctor>
struct ModViewCargo
{
    TFunctor    func;

    ModViewCargo() : func()
    {}
};

template <typename THost, typename TFunctor>
class ModifiedIterator<THost, ModView<TFunctor> >
{
public:
    typedef typename Cargo<ModifiedIterator>::Type TCargo_;

    THost _host;
    TCargo_ _cargo;

    mutable typename Value<ModifiedIterator>::Type tmp_value;

    ModifiedIterator() : _host(), tmp_value()
    {}

    template <typename TOtherHost>
    ModifiedIterator(ModifiedIterator<TOtherHost, ModView<TFunctor> > & origin) :
        _host(origin._host), _cargo(origin._cargo), tmp_value()
    {}

    explicit
    ModifiedIterator(THost const & host) :
        _host(host), tmp_value()
    {}

    ModifiedIterator(THost const & host, TFunctor const & functor):
        _host(host), tmp_value()
    {
        cargo(*this).func = functor;
    }
};

// --------------------------------------------------------------------------
// Class ModifiedString
// --------------------------------------------------------------------------

template <typename THost, typename TFunctor>
class ModifiedString<THost, ModView<TFunctor> >
{
public:
    typedef typename Pointer_<THost>::Type       THostPointer_;
    typedef typename Cargo<ModifiedString>::Type TCargo_;

    mutable THostPointer_ _host;
    TCargo_ _cargo;

    mutable typename Value<ModifiedString>::Type tmp_value;

    // Default constructor.
    ModifiedString() : _host(), tmp_value()
    {}

    // Construct with the actual host.
    explicit
    ModifiedString(typename Parameter_<THost>::Type host):
        _host(_toPointer(host)), tmp_value()
    {}

    // Construct with the functor.
    explicit
    ModifiedString(TFunctor const & functor):
        _host(), tmp_value()
    {
        cargo(*this).func = functor;
    }

    // Constructor for creating a ModifiedString with const host from a non-const host.
    template <typename THost_>
    explicit
    ModifiedString(THost_ & host,
                   SEQAN_CTOR_ENABLE_IF(IsConstructible<THost, THost_>)) :
            _host(_toPointer(host)), tmp_value()
    {
        ignoreUnusedVariableWarning(dummy);
    }

    // Construct with the actual host; variant with functor.
    ModifiedString(typename Parameter_<THost>::Type host, TFunctor const & functor) :
            _host(_toPointer(host)), tmp_value()
    {
        cargo(*this).func = functor;
    }

    // Constructor for creating a ModifiedString with const host with a non-const host; variant with functor.
    template <typename THost_>
    explicit
    ModifiedString(THost_ & host,
                   TFunctor const & functor,
                   SEQAN_CTOR_ENABLE_IF(IsConstructible<THost, THost_>)) :
            _host(_toPointer(host)), tmp_value()
    {
        ignoreUnusedVariableWarning(dummy);
        cargo(*this).func = functor;
    }

#ifdef SEQAN_CXX11_STANDARD

    // Constructor for innermost type; hand down to _host which is a ModifiedString itself.
    template <typename THost_>
    explicit
    ModifiedString(THost_ && host,
                   SEQAN_CTOR_ENABLE_IF(IsAnInnerHost<
                                            typename RemoveReference<THost>::Type,
                                            typename RemoveReference<THost_>::Type >)) :
            _host(std::forward<THost_>(host)), tmp_value()
    {
        ignoreUnusedVariableWarning(dummy);
    }

    // Constructor for innermost type; hand down to _host which is a ModifiedString itself.  Variant with functor.
    template <typename THost_>
    explicit
    ModifiedString(THost_ && host,
                   TFunctor const & functor,
                   SEQAN_CTOR_ENABLE_IF(IsAnInnerHost<
                                            typename RemoveReference<THost>::Type,
                                            typename RemoveReference<THost_>::Type >)) :
            _host(std::forward<THost_>(host)), tmp_value()
    {
        ignoreUnusedVariableWarning(dummy);
        cargo(*this).func = functor;
    }

#else

    // Constructor for innermost type; hand down to _host which is a ModifiedString itself.  Non-const variant.
    template <typename THost_>
    explicit
    ModifiedString(THost_ & host,
                   SEQAN_CTOR_ENABLE_IF(IsAnInnerHost<THost, THost_>)) :
            _host(host), tmp_value()
    {
        ignoreUnusedVariableWarning(dummy);
    }

    // Constructor for innermost type; hand down to _host which is a ModifiedString itself.  Non-const variant.
    template <typename THost_>
    explicit
    ModifiedString(THost_ const & host,
                   SEQAN_CTOR_ENABLE_IF(IsAnInnerHost<THost, THost_ const>)) :
            _host(host), tmp_value()
    {
        ignoreUnusedVariableWarning(dummy);
    }

    // Constructor for innermost type; hand down to _host which is a ModifiedString itself.  Non-const variant with
    // functor.
    template <typename THost_>
    explicit
    ModifiedString(THost_ & host,
                   TFunctor const & functor,
                   SEQAN_CTOR_ENABLE_IF(IsAnInnerHost<THost, THost_>)) :
            _host(host), tmp_value()
    {
        ignoreUnusedVariableWarning(dummy);
        cargo(*this).func = functor;
    }

    // Constructor for innermost type; hand down to _host which is a ModifiedString itself.  Non-const variant with
    // functor.
    template <typename THost_>
    explicit
    ModifiedString(THost_ const & host,
                   TFunctor const & functor,
                   SEQAN_CTOR_ENABLE_IF(IsAnInnerHost<THost, THost_ const>)) :
            _host(host), tmp_value()
    {
        ignoreUnusedVariableWarning(dummy);
        cargo(*this).func = functor;
    }

#endif

    template <typename TPos>
    inline typename Reference<ModifiedString>::Type
    operator[](TPos pos)
    {
        return value(*this, pos);
    }

    template <typename TPos>
    inline typename Reference<ModifiedString const>::Type
    operator[](TPos pos) const
    {
        return value(*this, pos);
    }
};

// ==========================================================================
// Metafunctions
// ==========================================================================

// --------------------------------------------------------------------------
// Metafunction Cargo                                      [ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost, typename TFunctor>
struct Cargo<ModifiedIterator<THost, ModView<TFunctor> > >
{
    typedef ModViewCargo<TFunctor>    Type;
};

// --------------------------------------------------------------------------
// Metafunction Value                                      [ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost, typename TFunctor>
struct Value<ModifiedIterator<THost, ModView<TFunctor> > >
{
    typedef typename TFunctor::result_type            TResult_;
    typedef typename RemoveConst_<TResult_>::Type   Type;
};

// --------------------------------------------------------------------------
// Metafunction GetValue                                   [ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost, typename TFunctor>
struct GetValue<ModifiedIterator<THost, ModView<TFunctor> > > : Value<ModifiedIterator<THost, ModView<TFunctor> > >
{};

// --------------------------------------------------------------------------
// Metafunction Reference                                  [ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost, typename TFunctor>
struct Reference<ModifiedIterator<THost, ModView<TFunctor> > >
{
    typedef typename Value<ModifiedIterator<THost, ModView<TFunctor> > >::Type & Type;
};

// --------------------------------------------------------------------------
// Metafunction Cargo                                        [ModifiedString]
// --------------------------------------------------------------------------

template <typename THost, typename TFunctor>
struct Cargo< ModifiedString<THost, ModView<TFunctor> > >
{
    typedef ModViewCargo<TFunctor>    Type;
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function value()                                        [ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost, typename TFunctor>
inline typename Reference<ModifiedIterator<THost, ModView<TFunctor> > >::Type
value(ModifiedIterator<THost, ModView<TFunctor> > & me)
{
    me.tmp_value = cargo(me).func(getValue(host(me)));
    return me.tmp_value;
}

template <typename THost, typename TFunctor>
inline typename Reference<ModifiedIterator<THost, ModView<TFunctor> > const>::Type
value(ModifiedIterator<THost, ModView<TFunctor> > const & me)
{
    me.tmp_value = cargo(me).func(getValue(host(me)));
    return me.tmp_value;
}

// --------------------------------------------------------------------------
// Function getValue()                                     [ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost, typename TFunctor>
inline typename GetValue<ModifiedIterator<THost, ModView<TFunctor> > >::Type
getValue(ModifiedIterator<THost, ModView<TFunctor> > & me)
{
    return cargo(me).func(getValue(host(me)));
}

template <typename THost, typename TFunctor>
inline typename GetValue<ModifiedIterator<THost, ModView<TFunctor> > const>::Type
getValue(ModifiedIterator<THost, ModView<TFunctor> > const & me)
{
    return cargo(me).func(getValue(host(me)));
}

// --------------------------------------------------------------------------
// Function value()                                          [ModifiedString]
// --------------------------------------------------------------------------

template <typename THost, typename TFunctor, typename TPos>
inline typename Reference<ModifiedString<THost, ModView<TFunctor> > >::Type
value(ModifiedString<THost, ModView<TFunctor> > & me, TPos pos)
{
    me.tmp_value = cargo(me).func(getValue(host(me), pos));
    return me.tmp_value;
}

template <typename THost, typename TFunctor, typename TPos>
inline typename Reference<ModifiedString<THost, ModView<TFunctor> > const>::Type
value(ModifiedString<THost, ModView<TFunctor> > const & me, TPos pos)
{
    me.tmp_value = cargo(me).func(getValue(host(me), pos));
    return me.tmp_value;
}

// --------------------------------------------------------------------------
// Function getValue()                                       [ModifiedString]
// --------------------------------------------------------------------------

template <typename THost, typename TFunctor, typename TPos>
inline typename GetValue<ModifiedString<THost, ModView<TFunctor> > >::Type
getValue(ModifiedString<THost, ModView<TFunctor> > & me, TPos pos)
{
    return cargo(me).func(getValue(host(me), pos));
}

template <typename THost, typename TFunctor, typename TPos>
inline typename GetValue<ModifiedString<THost, ModView<TFunctor> > const>::Type
getValue(ModifiedString<THost, ModView<TFunctor> > const & me, TPos pos)
{
    return cargo(me).func(getValue(host(me), pos));
}

// --------------------------------------------------------------------------
// Function assignModViewFunctor()
// --------------------------------------------------------------------------

template <typename THost, typename TFunctor>
inline void
assignModViewFunctor(ModifiedString<THost, ModView<TFunctor> > & me, TFunctor const & functor)
{
    cargo(me).func = functor;
}

// --------------------------------------------------------------------------
// Function convert()
// --------------------------------------------------------------------------

template < typename TSequence, typename TFunctor >
inline void
convert(TSequence & sequence, TFunctor const &F)
{
#if defined (_OPENMP) && defined (SEQAN_PARALLEL)
    // OpenMP does not support for loop with iterators. Therefore use index variables.
    typedef typename Position<TSequence>::Type    TPos;
    typedef typename MakeSigned_<TPos>::Type    TSignedPos;

    #pragma omp parallel for if(length(sequence) > 1000000)
    for(TSignedPos p = 0; p < (TSignedPos)length(sequence); ++p)
        sequence[p] = F(sequence[p]);

#else
    typedef typename Iterator<TSequence, Standard>::Type    TIter;

    TIter it = begin(sequence, Standard());
    TIter itEnd = end(sequence, Standard());
    for(; it != itEnd; ++it)
        *it = F(*it);
#endif
}

template < typename TSequence, typename TFunctor >
inline void
convert(TSequence const & sequence, TFunctor const &F)
{
#if defined (_OPENMP) && defined (SEQAN_PARALLEL)
    // OpenMP does not support for loop with iterators. Therefore use index variables.
    typedef typename Position<TSequence>::Type    TPos;
    typedef typename MakeSigned_<TPos>::Type    TSignedPos;

    #pragma omp parallel for if(length(sequence) > 1000000)
    for(TSignedPos p = 0; p < (TSignedPos)length(sequence); ++p)
        sequence[p] = F(sequence[p]);

#else
    typedef typename Iterator<TSequence const, Standard>::Type    TIter;

    TIter it = begin(sequence, Standard());
    TIter itEnd = end(sequence, Standard());
    for(; it != itEnd; ++it)
        *it = F(*it);
#endif
}

}  // namespace seqan

#endif  // SEQAN_MODIFIER_MODIFIER_VIEW_H_
