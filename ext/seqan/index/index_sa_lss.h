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
// This file contains the suffix array algorithm implementation
// by Larsson and Sadakane with the following modifications.
//
// MODIFICATIONS:
// - a context is used for reentrance
// - a generic TValue is used instead of int
// - functions are surrounded by SeqAn's namespace to omit namespace pollution
// ==========================================================================

/* qsufsort.c
   Copyright 1999, N. Jesper Larsson, all rights reserved.

   This file contains an implementation of the algorithm presented in "Faster
   Suffix Sorting" by N. Jesper Larsson (jesper@cs.lth.se) and Kunihiko
   Sadakane (sada@is.s.u-tokyo.ac.jp).

   This software may be used freely for any purpose. However, when distributed,
   the original source must be clearly stated, and, when the source code is
   distributed, the copyright notice must be retained and any alterations in
   the code must be clearly marked. No warranty is given regarding the quality
   of this software.*/

#ifndef SEQAN_HEADER_INDEX_SA_LSS_H
#define SEQAN_HEADER_INDEX_SA_LSS_H

namespace SEQAN_NAMESPACE_MAIN
{

template <typename TValue>
struct ContextLss_
{
    TValue *I,                   /* group array, ultimately suffix array.*/
    *V,                          /* inverse array, ultimately inverse of I.*/
    r,                           /* number of symbols aggregated by transform.*/
    h;                           /* length of already-sorted prefixes.*/

    // MODIFIED: renamed defines according to SeqAn's naming conventions
    #define SEQAN_LSSKEY(p)          (V[*(p)+(h)])
    #define SEQAN_LSSSWAP(p, q)      (tmp=*(p), *(p)=*(q), *(q)=tmp)
    #define SEQAN_LSSMED3(a, b, c)   (SEQAN_LSSKEY(a)<SEQAN_LSSKEY(b) ?                        \
            (SEQAN_LSSKEY(b)<SEQAN_LSSKEY(c) ? (b) : SEQAN_LSSKEY(a)<SEQAN_LSSKEY(c) ? (c) : (a))       \
            : (SEQAN_LSSKEY(b)>SEQAN_LSSKEY(c) ? (b) : SEQAN_LSSKEY(a)>SEQAN_LSSKEY(c) ? (c) : (a)))

    /* Subroutine for select_sort_split and sort_split. Sets group numbers for a
       group whose lowest position in I is pl and highest position is pm.*/

    inline void update_group(TValue *pl, TValue *pm)
    {
       TValue g;

       g=pm-I;                      /* group number.*/
       V[*pl]=g;                    /* update group number of first position.*/
       if (pl==pm)
          *pl=-1;                   /* one element, sorted group.*/
       else
          do                        /* more than one element, unsorted group.*/
             V[*++pl]=g;            /* update group numbers.*/
          while (pl<pm);
    }

    /* Quadratic sorting method to use for small subarrays. To be able to update
       group numbers consistently, a variant of selection sorting is used.*/

    inline void select_sort_split(TValue *p, TValue n) {
       TValue *pa, *pb, *pi, *pn;
       TValue f, v, tmp;

       pa=p;                        /* pa is start of group being picked out.*/
       pn=p+n-1;                    /* pn is last position of subarray.*/
       while (pa<pn) {
          for (pi=pb=pa+1, f=SEQAN_LSSKEY(pa); pi<=pn; ++pi)
             if ((v=SEQAN_LSSKEY(pi))<f) {
                f=v;                /* f is smallest key found.*/
                SEQAN_LSSSWAP(pi, pa);       /* place smallest element at beginning.*/
                pb=pa+1;            /* pb is position for elements equal to f.*/
             } else if (v==f) {     /* if equal to smallest key.*/
                SEQAN_LSSSWAP(pi, pb);       /* place next to other smallest elements.*/
                ++pb;
             }
          update_group(pa, pb-1);   /* update group values for new group.*/
          pa=pb;                    /* continue sorting rest of the subarray.*/
       }
       if (pa==pn) {                /* check if last part is single element.*/
          V[*pa]=pa-I;
          *pa=-1;                   /* sorted group.*/
       }
    }

    /* Subroutine for sort_split, algorithm by Bentley & McIlroy.*/

    inline TValue choose_pivot(TValue *p, TValue n) {
       TValue *pl, *pm, *pn;
       TValue s;

       pm=p+(n>>1);                 /* small arrays, middle element.*/
       if (n>7) {
          pl=p;
          pn=p+n-1;
          if (n>40) {               /* big arrays, pseudomedian of 9.*/
             s=n>>3;
             pl=SEQAN_LSSMED3(pl, pl+s, pl+s+s);
             pm=SEQAN_LSSMED3(pm-s, pm, pm+s);
             pn=SEQAN_LSSMED3(pn-s-s, pn-s, pn);
          }
          pm=SEQAN_LSSMED3(pl, pm, pn);      /* midsize arrays, median of 3.*/
       }
       return SEQAN_LSSKEY(pm);
    }

    /* Sorting routine called for each unsorted group. Sorts the array of integers
       (suffix numbers) of length n starting at p. The algorithm is a ternary-split
       quicksort taken from Bentley & McIlroy, "Engineering a Sort Function",
       Software -- Practice and Experience 23(11), 1249-1265 (November 1993). This
       function is based on Program 7.*/

    inline void sort_split(TValue *p, TValue n)
    {
       TValue *pa, *pb, *pc, *pd, *pl, *pm, *pn;
       TValue f, v, s, t, tmp;

       if (n<7) {                   /* multi-selection sort smallest arrays.*/
          select_sort_split(p, n);
          return;
       }

       v=choose_pivot(p, n);
       pa=pb=p;
       pc=pd=p+n-1;
       while (1) {                  /* split-end partition.*/
          while (pb<=pc && (f=SEQAN_LSSKEY(pb))<=v) {
             if (f==v) {
                SEQAN_LSSSWAP(pa, pb);
                ++pa;
             }
             ++pb;
          }
          while (pc>=pb && (f=SEQAN_LSSKEY(pc))>=v) {
             if (f==v) {
                SEQAN_LSSSWAP(pc, pd);
                --pd;
             }
             --pc;
          }
          if (pb>pc)
             break;
          SEQAN_LSSSWAP(pb, pc);
          ++pb;
          --pc;
       }
       pn=p+n;
       if ((s=pa-p)>(t=pb-pa))
          s=t;
       for (pl=p, pm=pb-s; s; --s, ++pl, ++pm)
          SEQAN_LSSSWAP(pl, pm);
       if ((s=pd-pc)>(t=pn-pd-1))
          s=t;
       for (pl=pb, pm=pn-s; s; --s, ++pl, ++pm)
          SEQAN_LSSSWAP(pl, pm);

       s=pb-pa;
       t=pd-pc;
       if (s>0)
          sort_split(p, s);
       update_group(p+s, p+n-t-1);
       if (t>0)
          sort_split(p+n-t, t);
    }

    /* Bucketsort for first iteration.

       Input: x[0...n-1] holds integers in the range 1...k-1, all of which appear
       at least once. x[n] is 0. (This is the corresponding output of transform.) k
       must be at most n+1. p is array of size n+1 whose contents are disregarded.

       Output: x is V and p is I after the initial sorting stage of the refined
       suffix sorting algorithm.*/

    inline void bucketsort(TValue *x, TValue *p, TValue n, TValue k)
    {
       TValue *pi, i, c, d, g;

       for (pi=p; pi<p+k; ++pi)
          *pi=-1;                   /* mark linked lists empty.*/
       for (i=0; i<=n; ++i) {
          x[i]=p[c=x[i]];           /* insert in linked list.*/
          p[c]=i;
       }
       for (pi=p+k-1, i=n; pi>=p; --pi) {
          d=x[c=*pi];               /* c is position, d is next in list.*/
          x[c]=g=i;                 /* last position equals group number.*/
          if (d>=0) {               /* if more than one element in group.*/
             p[i--]=c;              /* p is permutation for the sorted x.*/
             do {
                d=x[c=d];           /* next in linked list.*/
                x[c]=g;             /* group number in x.*/
                p[i--]=c;           /* permutation in p.*/
             } while (d>=0);
          } else
             p[i--]=-1;             /* one element, sorted group.*/
       }
    }

    /* Transforms the alphabet of x by attempting to aggregate several symbols into
       one, while preserving the suffix order of x. The alphabet may also be
       compacted, so that x on output comprises all integers of the new alphabet
       with no skipped numbers.

       Input: x is an array of size n+1 whose first n elements are positive
       integers in the range l...k-1. p is array of size n+1, used for temporary
       storage. q controls aggregation and compaction by defining the maximum value
       for any symbol during transformation: q must be at least k-l; if q<=n,
       compaction is guaranteed; if k-l>n, compaction is never done; if q is
       INT_MAX, the maximum number of symbols are aggregated into one.

       Output: Returns an integer j in the range 1...q representing the size of the
       new alphabet. If j<=n+1, the alphabet is compacted. The global variable r is
       set to the number of old symbols grouped into one. Only x[n] is 0.*/

    inline TValue transform(TValue *x, TValue *p, TValue n, TValue k, TValue l, TValue q)
    {
       TValue b, c, d, e, i, j, m, s;
       TValue *pi, *pj;

       for (s=0, i=k-l; i; i>>=1)
          ++s;                      /* s is number of bits in old symbol.*/
       e=MaxValue<TValue>::VALUE>>s; /* e is for overflow checking.*/
       for (b=d=r=0; r<n && d<=e && (c=d<<s|(k-l))<=q; ++r) {
          b=b<<s|(x[r]-l+1);        /* b is start of x in chunk alphabet.*/
          d=c;                      /* d is max symbol in chunk alphabet.*/
       }
       m=(1<<(r-1)*s)-1;            /* m masks off top old symbol from chunk.*/
       x[n]=l-1;                    /* emulate zero terminator.*/
       if (d<=n) {                  /* if bucketing possible, compact alphabet.*/
          for (pi=p; pi<=p+d; ++pi)
             *pi=0;                 /* zero transformation table.*/
          for (pi=x+r, c=b; pi<=x+n; ++pi) {
             p[c]=1;                /* mark used chunk symbol.*/
             c=(c&m)<<s|(*pi-l+1);  /* shift in next old symbol in chunk.*/
          }
          for (i=1; i<r; ++i) {     /* handle last r-1 positions.*/
             p[c]=1;                /* mark used chunk symbol.*/
             c=(c&m)<<s;            /* shift in next old symbol in chunk.*/
          }
          for (pi=p, j=1; pi<=p+d; ++pi)
             if (*pi)
                *pi=j++;            /* j is new alphabet size.*/
          for (pi=x, pj=x+r, c=b; pj<=x+n; ++pi, ++pj) {
             *pi=p[c];              /* transform to new alphabet.*/
             c=(c&m)<<s|(*pj-l+1);  /* shift in next old symbol in chunk.*/
          }
          while (pi<x+n) {          /* handle last r-1 positions.*/
             *pi++=p[c];            /* transform to new alphabet.*/
             c=(c&m)<<s;            /* shift right-end zero in chunk.*/
          }
       } else {                     /* bucketing not possible, don't compact.*/
          for (pi=x, pj=x+r, c=b; pj<=x+n; ++pi, ++pj) {
             *pi=c;                 /* transform to new alphabet.*/
             c=(c&m)<<s|(*pj-l+1);  /* shift in next old symbol in chunk.*/
          }
          while (pi<x+n) {          /* handle last r-1 positions.*/
             *pi++=c;               /* transform to new alphabet.*/
             c=(c&m)<<s;            /* shift right-end zero in chunk.*/
          }
          j=d+1;                    /* new alphabet size.*/
       }
       x[n]=0;                      /* end-of-string symbol is zero.*/
       return j;                    /* return new alphabet size.*/
    }

    /* Makes suffix array p of x. x becomes inverse of p. p and x are both of size
       n+1. Contents of x[0...n-1] are integers in the range l...k-1. Original
       contents of x[n] is disregarded, the n-th symbol being regarded as
       end-of-string smaller than all other symbols.*/

    void suffixsort(TValue *x, TValue *p, TValue n, TValue k, TValue l)
    {
       TValue *pi, *pk;
       TValue i, j, s, sl;

       V=x;                         /* set global values.*/
       I=p;

       if (n>=k-l) {                /* if bucketing possible,*/
          j=transform(V, I, n, k, l, n);
          bucketsort(V, I, n, j);   /* bucketsort on first r positions.*/
       } else {
          transform(V, I, n, k, l, MaxValue<TValue>::VALUE);
          for (i=0; i<=n; ++i)
             I[i]=i;                /* initialize I with suffix numbers.*/
          h=0;
          sort_split(I, n+1);       /* quicksort on first r positions.*/
       }
       h=r;                         /* number of symbols aggregated by transform.*/

       while (*I>=-n) {
          pi=I;                     /* pi is first position of group.*/
          sl=0;                     /* sl is negated length of sorted groups.*/
          do {
             if ((s=*pi)<0) {
                pi-=s;              /* skip over sorted group.*/
                sl+=s;              /* add negated length to sl.*/
             } else {
                if (sl) {
                   *(pi+sl)=sl;     /* combine sorted groups before pi.*/
                   sl=0;
                }
                pk=I+V[s]+1;        /* pk-1 is last position of unsorted group.*/
                sort_split(pi, pk-pi);
                pi=pk;              /* next group.*/
             }
          } while (pi<=I+n);
          if (sl)                   /* if the array ends with a sorted group.*/
             *(pi+sl)=sl;           /* combine sorted groups at end of I.*/
          h=2*h;                    /* double sorted-depth.*/
       }

       for (i=0; i<=n; ++i)         /* reconstruct suffix array from inverse.*/
          I[V[i]]=i;
    }
};



//////////////////////////////////////////////////////////////////////////////
// SeqAn interface

    struct LarssonSadakane {};

    // WARNING:
    // 1. the content of s will be destroyed
    // 2. s and SA have to have size n+1
    // 3. s[n] must be unique and less than any other character in s
    //
    // better use LarssonSadakane as a pipe (look down)

    template < typename TSA,
               typename TText >
    void createSuffixArray(
        TSA &SA,
        TText const &s,
        LarssonSadakane const &,
        unsigned K)
    {
        typedef typename Value<TSA>::Type            TValue;
        typedef typename MakeSigned_<TValue>::Type    TSValue;    // LarssonSadakane expects signed values
        ContextLss_<TSValue> c;
        c.suffixsort(
            (TSValue*)begin(s, Standard()),        // text
            (TSValue*)begin(SA, Standard()),    // SA
            length(s) - 1,                        // n
            K,                                    // text[i] <  K
            0);                                    // text[i] >= 0
    }

    //////////////////////////////////////////////////////////////////////////////
    // qsufsort pipe
    template < typename TInput >
    struct Pipe< TInput, LarssonSadakane >
    {
        typedef typename SAValue<TInput>::Type    TValue;
        typedef String<TValue, Alloc<> >        TText;
        typedef String<TValue, Alloc<> >        TSA;
        typedef Pipe<TSA, Source<> >            TSource;

        TSA        sa;
        TSource    in;

        Pipe(TInput &_textIn):
            in(sa)
        {
            typedef typename Iterator<TText, Standard>::Type TIter;

            TValue len = length(_textIn);
            TText text;
            resize(text, len + 1, Exact());

            TIter    it = begin(text);
            TIter    itEnd = begin(text) + len;
            TValue    maxChar = 0;

            beginRead(_textIn);
            for(; it != itEnd; ++it, ++_textIn) {
                TValue val = (*it = 1 + (TValue)*_textIn);
                if (val > maxChar)
                    maxChar = val;
            }
            endRead(_textIn);

            resize(sa, len + 1, Exact());
            createSuffixArray(sa, text, LarssonSadakane(), maxChar + 1);
        }

        inline typename Value<TSource>::Type const & operator*() {
            return *in;
        }

        inline Pipe& operator++() {
            ++in;
            return *this;
        }
    };

    template < typename TInput >
    inline bool control(Pipe< TInput, LarssonSadakane > &me, ControlBeginRead const &command) {
        control(me.in, command);
        ++me.in;
        return true;
    }

    template < typename TInput >
    inline typename Size< Pipe< TInput, LarssonSadakane > >::Type
    length(Pipe< TInput, LarssonSadakane > const &me) {
        return length(me.in) - 1;
    }

}

#endif
