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

// TODO(holtgrew): We should probably specialize Automaton for this.

#ifndef SEQAN_HEADER_GRAPH_IMPL_ORACLE_H
#define SEQAN_HEADER_GRAPH_IMPL_ORACLE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Graph - Oracle
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TPropertyMap, typename TChar>
inline void
_addLetterToOracle(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
                   TPropertyMap& supplyState,
                   TChar const c)
{
    typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
    TVertexDescriptor newState = addVertex(g);
    TVertexDescriptor pred = newState - 1;
    addEdge(g, pred, newState, c);
    TVertexDescriptor k = getProperty(supplyState, pred);
    while ((k!=nilVal) &&
            (getTarget(&g.data_vertex[k].data_edge[ordValue(TAlphabet(c))])==nilVal))
    {
        addEdge(g,k,newState,c);
        k = getProperty(supplyState, k);
    }
    TVertexDescriptor s;
    if (k==nilVal) s=0;
    else s = getTarget(&g.data_vertex[k].data_edge[ordValue(TAlphabet(c))]);
    assignProperty(supplyState, newState, s);
}

//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Automaton#createOracle
 * @brief Creates a factor oracle into an Automaton object.
 *
 * @signature void createOracle(g, text);
 *
 * @param[out] g    The oracle to create.
 * @param[in]  text A @link String @endlink to create the oracle with.
 *
 * @see Automaton#createOracleOnReverse
 */

template<typename TAlphabet, typename TCargo, typename TSpec, typename TText>
inline void
createOracle(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
             TText const text)
{
    SEQAN_CHECKPOINT
    typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Size<TText>::Type TSize;
    TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
    TSize len = length(text);
    String<TVertexDescriptor> supplyState;
    resize(supplyState, len+1);
    TVertexDescriptor v1 = addVertex(g);
    assignRoot(g,v1);
    assignProperty(supplyState, v1, nilVal);
    for(TSize i = 0; i<len; ++i) _addLetterToOracle(g, supplyState, getValue(text,i));
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Automaton#createOracleOnReverse
 * @brief Creates a factor oracle on the reversed string.
 *
 * @signature void createOracleOnReverse(g, text);
 *
 * @param[out] g    The oracle to create.
 * @param[in]  text A @link String @endlink to create the oracle with.
 *
 * @see Automaton#createOracle
 */

template<typename TAlphabet, typename TCargo, typename TSpec, typename TText>
inline void
createOracleOnReverse(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
                      TText const text)
{
    SEQAN_CHECKPOINT
    typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Size<TText>::Type TSize;
    TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
    TSize len = length(text);
    String<TVertexDescriptor> supplyState;
    resize(supplyState, len+1);
    TVertexDescriptor v1 = addVertex(g);
    assignRoot(g,v1);
    assignProperty(supplyState, v1, nilVal);
    for(TSize i = len-1; i>0; --i) _addLetterToOracle(g, supplyState, getValue(text,i));
    _addLetterToOracle(g, supplyState, getValue(text,0));
}


//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TTerminalStateMap, typename TKeywords>
inline void
createSetOracle(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
                TTerminalStateMap& terminalStateMap,
                TKeywords const& keywords)
{
    typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Position<TKeywords>::Type TPos;
    typedef typename Value<TKeywords>::Type TKeyword;
    typedef typename Iterator<TKeyword const, Standard>::Type TIterator;

    createTrie(g, terminalStateMap, keywords);

    String<TVertexDescriptor> supplyState;
    resizeVertexMap(supplyState, g);
    String<bool> visited;
    resizeVertexMap(visited, g);
    arrayFill(begin(visited), end(visited), false);

    TVertexDescriptor nil_ = getNil<TVertexDescriptor>();
    assignProperty(supplyState, root(g), nil_);

    TVertexDescriptor _root = getRoot(g);

    TPos len = length(keywords);
    String<TVertexDescriptor> _here_v;
    resize(_here_v, len, _root);
    String<TIterator> _here_it;
    resize(_here_it, len);
    for (TPos i = 0; i < len; ++i)
    {
        _here_it[i] = begin(keywords[i], Standard());
    }
    TPos _active_count = len;
    while (_active_count)
    {
        for (TPos i = 0; i < len; ++i)
        {
            TIterator & it = _here_it[i];
            TIterator it_end = end(keywords[i], Standard());
            TVertexDescriptor & _parent = _here_v[i];

            if (it != it_end)
            {
                TVertexDescriptor _current = getSuccessor(g, _parent, *it);

                if (!getProperty(visited, _current))
                {
                    assignProperty(visited, _current, true);

                    TVertexDescriptor _down = getProperty(supplyState, _parent);
                    TVertexDescriptor _supply = _root;
                    while (_down != nil_)
                    {
                        TVertexDescriptor _next = getSuccessor(g, _down, *it);
                        if (_next != nil_)
                        {
                            _supply = _next;
                            break;
                        }

                        addEdge(g, _down, _current, *it);
                        _down = getProperty(supplyState, _down);
                    }
                    assignProperty(supplyState, _current, _supply);
                }
                _parent = _current;
                ++it;
                if (it == it_end)
                {
                    --_active_count;
                }
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
