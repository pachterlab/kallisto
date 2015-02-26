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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Implementation of the biological SimpleType specializations Dna, Dna5,
// DnaQ, Dna5Q, Rna, Rna5, Iupac, and AminoAcid.  The conversion tables are
// in alphabet_residue_tabs.h.
//
// This header's structure is an exception to the standard.  Because
// splitting into one header for each specialization is a bit too much, we
// define all types in one header.  We define the classes, metafunctions and
// functions in one section, one subsection for reach type to make the whole
// thing more readable.  Conversion through assignment is defined in the
// Function section.
// ==========================================================================

// TODO(holtgrew): Add RnaQ and Rna5Q? Can we create a tag/type for Dna and Rna that is then differentiated with one additional tag?

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_BASIC_ALPHABET_RESIDUE_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_BASIC_ALPHABET_RESIDUE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename T> struct BaseAlphabet;

// ============================================================================
// Classes, Metafunctions, Functions
// ============================================================================

// Also see comment at the top for more information on this exceptional
// structure.
//
// We define the SimpleType specializations, the metafunctions ValueSize and
// BitsPerValue and the functions unknownValueImpl() for each specialization in
// this section.

// ----------------------------------------------------------------------------
// Specialization Dna
// ----------------------------------------------------------------------------

/*!
 * @class Dna
 * @extends SimpleType
 * @headerfile <seqan/basic.h>
 * @brief Alphabet for DNA.
 *
 * @signature typedef SimpleType<unsigned char, Dna_> Dna;
 *
 * The ValueSize of <tt>Dna</tt> is 4.  The nucleotides are enumerated this way: <tt>'A' = 0, 'C' = 1, 'G' = 2, 'T' =
 * 3</tt>.
 *
 * Objects of type <tt>Dna</tt> can be converted to various other types and vice versa.  An object that has a value not
 * in <tt>{'A', 'C', 'G', 'T'}</tt> is converted to <tt>'A'</tt>.
 *
 * @see Dna5
 * @see DnaString
 * @see DnaIterator
 */

struct Dna_ {};
typedef SimpleType<unsigned char, Dna_> Dna;

template <>
struct ValueSize<Dna>
{
    typedef __uint8 Type;
    static const Type VALUE = 4;
};

template <>
struct BitsPerValue< Dna >
{
    typedef __uint8 Type;
    static const Type VALUE = 2;
};

// ----------------------------------------------------------------------------
// Specialization Dna5
// ----------------------------------------------------------------------------

/*!
 * @class Dna5
 * @extends SimpleType
 * @headerfile <seqan/basic.h>
 * @brief Alphabet for DNA including 'N' character.
 *
 * @signature typedef SimpleType<unsigned char, Dna5_> Dna5;
 *
 * The @link FiniteOrderedAlphabetConcept#ValueSize @endlink of <tt>Dna5</tt> is 5.  The nucleotides are enumerated this
 * way: <tt>'A' = 0, 'C' = 1, 'G' = 2, 'T' = 3</tt>.  The 'N' character ("unkown nucleotide") is encoded by 4.
 *
 * Objects of type <tt>Dna5</tt> can be converted to various other types and vice versa.  An object that has a value not
 * in <tt>{'A', 'C', 'G', 'T'}</tt> is converted to <tt>'N'</tt>.
 *
 * @see Dna5Iterator
 * @see Dna5String
 * @see Dna
 */

struct Dna5_ {};
typedef SimpleType<unsigned char, Dna5_> Dna5;

template <>
struct ValueSize<Dna5>
{
    typedef __uint8 Type;
    static const Type VALUE = 5;
};

template <>
struct BitsPerValue<Dna5>
{
    typedef __uint8 Type;
    static const Type VALUE = 3;
};

inline Dna5
unknownValueImpl(Dna5 *)
{
    static const Dna5 _result = Dna5('N');
    return _result;
}

// ----------------------------------------------------------------------------
// Specialization DnaQ
// ----------------------------------------------------------------------------

/*!
 * @class DnaQ
 * @extends SimpleType
 * @headerfile <seqan/basic.h>
 * @implements AlphabetWithQualitiesConcept
 * @brief Alphabet for DNA plus PHRED quality.
 *
 * @signature typedef SimpleType<unsigned char, DnaQ_> DnaQ;
 *
 * The ValueSize of <tt>DnaQ</tt> is 4.  The nucleotides are enumerated this way: <tt>'A' = 0, 'C' = 1, 'G' = 2, 'T' =
 * 3</tt>.
 *
 * Objects of type <tt>DnaQ</tt> can be converted to various other types and vice versa.
 *
 * Note that the default quality value is set to 60.
 *
 * @see Dna5Q
 * @see Dna
 */

#ifndef SEQAN_DEFAULT_QUALITY
#define SEQAN_DEFAULT_QUALITY 40
#endif

struct DnaQ_ {};
typedef SimpleType <unsigned char, DnaQ_> DnaQ;

template <> struct ValueSize<DnaQ>
{
    typedef __uint8 Type;
    static const ValueSize<DnaQ>::Type VALUE = 4;  // Considering nucleotides.
};

template <> struct InternalValueSize_<DnaQ>
{
    enum { VALUE = 252 };  // Considering nucleotides x Quality 0..62.
};

template <> struct BitsPerValue<DnaQ>
{
    enum { VALUE = 8 };
    typedef __uint8 Type;
};

template <> struct HasQualities<DnaQ>
{
    typedef True Type;
    static const bool VALUE = true;
};

template <>
struct BaseAlphabet<DnaQ>
{
    typedef Dna Type;
};

template <>
struct QualityValueSize<DnaQ>
{
    enum { VALUE = 63 }; // 64 - 1 (N)
};

template <typename TValue>
inline int getQualityValue(TValue const &)
{
    return 0;
}

inline int getQualityValue(DnaQ const & c)
{
    return c.value >> 2;
}

inline
void assignQualityValue(DnaQ & c, int q)
{
    if (q < 0) q = 0;
    if (q >= QualityValueSize<DnaQ>::VALUE)
        q = QualityValueSize<DnaQ>::VALUE - 1;
    c.value = (c.value & 3) | (q << 2);
}

inline
void assignQualityValue(DnaQ & c, char q)
{
    int q1 = static_cast<int>(q - '!');
    if (q1 < 0) q1 = 0;
    if (q1 >= QualityValueSize<DnaQ>::VALUE)
        q1 = QualityValueSize<DnaQ>::VALUE - 1;
    assignQualityValue(c, q1);
}

inline
void assignQualityValue(char & q, DnaQ c)
{
    q = '!' + getQualityValue(c);
}


// ----------------------------------------------------------------------------
// Specialization Dna5Q
// ----------------------------------------------------------------------------

/*!
 * @class Dna5Q
 * @extends SimpleType
 * @headerfile <seqan/basic.h>
 * @implements AlphabetWithQualitiesConcept
 * @brief Alphabet for DNA plus PHRED quality including 'N' character.
 *
 * @signature typedef SimpleType<unsigned char, Dna5Q_> Dna5Q;
 *
 * The ValueSize of <tt>Dna5Q</tt> is 5.  The nucleotides are enumerated this way: <tt>'A' = 0, 'C' = 1, 'G' = 2, 'T' =
 * 3</tt>. The 'N' character ("unknown nucleotide") is encoded by 4.
 *
 * Objects of type <tt>Dna5</tt> can be converted to various other types and vice versa.
 *
 * Note that the default quality value is set to 40.
 *
 * @see Dna5
 * @see DnaQ
 */

struct Dna5Q_ {};
typedef SimpleType <unsigned char, Dna5Q_> Dna5Q;

static const unsigned char Dna5QValueN_ = 252;                              // value representing N

template <> struct ValueSize<Dna5Q>
{
    typedef __uint8 Type;
    static const Type VALUE = 5;  // Considering nucleotides + N.
};

template <> struct InternalValueSize_<Dna5Q>
{
    enum { VALUE = 253 };  // Considering (nucleotides x Quality 0..62) + N.
};

template <> struct BitsPerValue<Dna5Q>
{
    enum { VALUE = 8 };
    typedef __uint8 Type;
};

template <> struct HasQualities<Dna5Q>
{
    typedef True Type;
    static const bool VALUE = true;
};

template <>
struct BaseAlphabet<Dna5Q>
{
    typedef Dna5 Type;
};

template <> struct
QualityValueSize<Dna5Q>
{
    enum { VALUE = 63 };
};

inline Dna5Q
unknownValueImpl(Dna5Q *)
{
    static const Dna5Q _result = Dna5Q('N');
    return _result;
}

inline int getQualityValue(Dna5Q const &c)
{
    // We use a lookup table to extract the qualities from DNA5Q.  The lookup
    // table based code is equivalent to the following line:
    // return (c.value == Dna5QValueN_)? 0: c.value >> 2;

    static const unsigned table[] = {
         0,  0,  0,  0,  1,  1,  1,  1,  2,  2,  2,  2,  3,  3,  3,  3,  4,
         4,  4,  4,  5,  5,  5,  5,  6,  6,  6,  6,  7,  7,  7,  7,  8,  8,
         8,  8,  9,  9,  9,  9, 10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12,
        12, 13, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15, 15, 16, 16, 16, 16,
        17, 17, 17, 17, 18, 18, 18, 18, 19, 19, 19, 19, 20, 20, 20, 20, 21,
        21, 21, 21, 22, 22, 22, 22, 23, 23, 23, 23, 24, 24, 24, 24, 25, 25,
        25, 25, 26, 26, 26, 26, 27, 27, 27, 27, 28, 28, 28, 28, 29, 29, 29,
        29, 30, 30, 30, 30, 31, 31, 31, 31, 32, 32, 32, 32, 33, 33, 33, 33,
        34, 34, 34, 34, 35, 35, 35, 35, 36, 36, 36, 36, 37, 37, 37, 37, 38,
        38, 38, 38, 39, 39, 39, 39, 40, 40, 40, 40, 41, 41, 41, 41, 42, 42,
        42, 42, 43, 43, 43, 43, 44, 44, 44, 44, 45, 45, 45, 45, 46, 46, 46,
        46, 47, 47, 47, 47, 48, 48, 48, 48, 49, 49, 49, 49, 50, 50, 50, 50,
        51, 51, 51, 51, 52, 52, 52, 52, 53, 53, 53, 53, 54, 54, 54, 54, 55,
        55, 55, 55, 56, 56, 56, 56, 57, 57, 57, 57, 58, 58, 58, 58, 59, 59,
        59, 59, 60, 60, 60, 60, 61, 61, 61, 61, 62, 62, 62, 62,
        0,  0,  0,  0};
    return table[c.value];
}

inline
void assignQualityValue(Dna5Q &c, int q)
{
    if (q < 0) q = 0;
    if (q >= QualityValueSize<Dna5Q>::VALUE)
        q = QualityValueSize<Dna5Q>::VALUE - 1;
    if (c.value != Dna5QValueN_)
        c.value = (c.value & 3) | (q << 2);
}

inline
void assignQualityValue(Dna5Q &c, char q)
{
    int q1 = static_cast<int>(q - '!');
    if (q1 < 0) q1 = 0;
    if (q1 >= QualityValueSize<Dna5Q>::VALUE)
        q1 = QualityValueSize<Dna5Q>::VALUE - 1;
    assignQualityValue(c, q1);
}

inline
void assignQualityValue(char & q, Dna5Q c)
{
    q = '!' + getQualityValue(c);
}

// ----------------------------------------------------------------------------
// Specialization Rna
// ----------------------------------------------------------------------------

/*!
 * @class Rna
 * @extends SimpleType
 * @headerfile <seqan/basic.h>
 * @brief Alphabet for RNA.
 *
 * @signature typedef SimpleType<unsigned char, Rna_> Rna;
 *
 * The ValueSize of <tt>Rna</tt> is 4.  The nucleotides are enumerated this way: <tt>'A' = 0, 'C' = 1, 'G' = 2, 'U' =
 * 3</tt>.
 *
 * Objects of type <tt>Rna</tt> can be converted to various other types and vice versa.  An object that has a value not
 * in <tt>{'A', 'C', 'G', 'U'}</tt> is converted to <tt>'A'</tt>.
 *
 * <tt>Rna</tt> is typedef for <tt>SimpleType<char,Rna_></tt>, while <tt>Rna_</tt> is a helper specialization tag class.
 *
 * @see Rna5
 * @see RnaString
 * @see RnaIterator
 */

struct Rna_ {};
typedef SimpleType<unsigned char, Rna_> Rna;

template <>
struct ValueSize<Rna>
{
    typedef __uint8 Type;
    static const Type VALUE = 4;
};

template <>
struct BitsPerValue<Rna>
{
    typedef __uint8 Type;
    static const Type VALUE = 2;
};

// ----------------------------------------------------------------------------
// Specialization Rna5
// ----------------------------------------------------------------------------

/*!
 * @class Rna5
 * @extends SimpleType
 * @headerfile <seqan/basic.h>
 * @brief Alphabet for RNA including 'N' character.
 *
 * @signature typedef SimpleType<unsigned char, Rna5_> Rna5;
 *
 * The ValueSize of <tt>Rna5</tt> is 5.  The nucleotides are enumerated this way: <tt>'A' = 0, 'C' = 1, 'G' = 2, 'U' =
 * 3</tt>.  The 'N' character ("unkown nucleotide") is encoded by 4.
 *
 * Objects of type <tt>Rna5</tt> can be converted to various other types and vice versa.  An object that has a value not
 * in <tt>{'A', 'C', 'G', 'U'}</tt> is converted to <tt>'N'</tt>.
 *
 * @see Rna5Iterator
 * @see Rna5String
 * @see Rna
 */

struct Rna5_ {};
typedef SimpleType<unsigned char, Rna5_> Rna5;

template <>
struct ValueSize<Rna5>
{
    typedef __uint8 Type;
    static const Type VALUE = 5;
};

template <> struct BitsPerValue<Rna5>
{
    typedef __uint8 Type;
    static const Type VALUE = 3;
};

inline Rna5
unknownValueImpl(Rna5 *)
{
    static const Rna5 _result = Rna5('N');
    return _result;
}

// ----------------------------------------------------------------------------
// Specialization Iupac
// ----------------------------------------------------------------------------

// TODO(holtgrew): We should support retrieval of nucleotides represented by a IUPAC char.

/*!
 * @class Iupac
 * @extends SimpleType
 * @headerfile <seqan/basic.h>
 * @brief Iupac code for DNA.
 *
 * @signature typedef SimpleType<unsigned char, Iupac_> Iupac;
 *
 * The ValueSize of <tt>Iupac</tt> is 16.  The nucleotides are enumerated from 0 to 19 in this order: 'U'=0, 'T', 'A',
 * 'W', 'C', 'Y', 'M', 'H', 'G', 'K', 'R', 'D', 'S', 'B', 'V', 'N'=15.
 *
 * Objects of type <tt>Iupac</tt> can be converted to various other types and vice versa.  Unknown values are converted
 * to <tt>'N'</tt>.
 *
 * @see IupacString
 * @see IupacIterator
 */

struct Iupac_ {};
typedef SimpleType<unsigned char, Iupac_> Iupac;

template <> struct ValueSize<Iupac>
{
    typedef __uint8 Type;
    static const Type VALUE = 16;
};

template <> struct BitsPerValue<Iupac>
{
    typedef __uint8 Type;
    static const Type VALUE = 4;
};

inline Iupac
unknownValueImpl(Iupac *)
{
    static const Iupac _result = Iupac('N');
    return _result;
}

// ----------------------------------------------------------------------------
// Specialization AminoAcid
// ----------------------------------------------------------------------------

/*!
 * @class AminoAcid
 * @extends SimpleType
 * @headerfile <seqan/basic.h>
 * @brief IUPAC code for amino acids.
 * @signature typedef SingleType<unsigned char, AminoAcid_> AminoAcid;
 *
 * The ValueSize of <tt>AminoAcid</tt> is 26.
 *
 * The amino acid symbols are as follows, i.e. they are sorted alphabetically
 * up until the last two symbols:
 *
 * 'A' = 0, 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'Y', 'Z', 'X'. '*'
 *
 * Of these 'B' is a wildcard for (Aspartic Acid, Asparagine),
 * 'J' for (Leucine, Isoleucine), 'Z' for (Glutamic Acid, Glutamine) and
 * 'X' for "any amino acid".
 *
 * 'U' refers to the rare Selenocysteine and '*' to the terminator tRNA.
 *
 * Objects of type <tt>AminoAcid</tt> can be converted to <tt>char</tt> and vice versa.  Unknown values are converted to
 * <tt>'X'</tt>.
 *
 * @see FiniteOrderedAlphabetConcept#ValueSize
 * @see PeptideIterator
 * @see Peptide
 */

struct AminoAcid_ {};
typedef SimpleType<unsigned char, AminoAcid_> AminoAcid;

template <> struct ValueSize<AminoAcid>
{
    typedef __uint8 Type;
    static const Type VALUE = 26;
};

template <> struct BitsPerValue<AminoAcid>
{
    typedef __uint8 Type;
    static const Type VALUE = 5;
};

inline AminoAcid
unknownValueImpl(AminoAcid *)
{
    static const AminoAcid _result = AminoAcid('X');
    return _result;
}

// ----------------------------------------------------------------------------
// Specialization Finite
// ----------------------------------------------------------------------------

/*!
 * @class Finite
 * @extends SimpleType
 * @headerfile <seqan/basic.h>
 *
 * @brief A finite alphabet of a fixed size.
 *
 * @signature template <typename TValue, unsigned SIZE>
 *            class SimpleType<TValue, Finite<SIZE> >;
 *
 * @tparam TValue The type that is use to store the values.
 * @tparam SIZE   The ValueSize of the alphabet.
 */

template <unsigned SIZE>
struct Finite;

template <typename TValue, unsigned SIZE>
struct ValueSize<SimpleType<TValue, Finite<SIZE> > >
{
    typedef __uint8 Type;
    static const Type VALUE = SIZE;
};

template <typename TValue, unsigned SIZE>
struct BitsPerValue<SimpleType<TValue, Finite<SIZE> > >
{
    typedef __uint8 Type;
    static const Type VALUE = Log2<SIZE>::VALUE;
};

// ============================================================================
// Assignment / Conversion Functions
// ============================================================================

// ----------------------------------------------------------------------------
// char
// ----------------------------------------------------------------------------

inline void assign(char & c_target,
                   Dna const & source)
{
    c_target = TranslateTableDna5ToChar_<>::VALUE[source.value];
}

inline void assign(char & c_target,
                   Dna5 const & source)
{
    c_target = TranslateTableDna5ToChar_<>::VALUE[source.value];
}

inline void assign(char & c_target,
                   Rna const & source)
{
    c_target = TranslateTableRna5ToChar_<>::VALUE[source.value];
}

inline void assign(char & c_target,
                   Rna5 const & source)
{
    c_target = TranslateTableRna5ToChar_<>::VALUE[source.value];
}

inline void assign(char & c_target, Iupac const & source)
{
    c_target = TranslateTableIupacToChar_<>::VALUE[source.value];
}

inline void assign(char & c_target, AminoAcid const & source)
{
    c_target = TranslateTableAAToChar_<>::VALUE[source.value];
}

// ----------------------------------------------------------------------------
// Dna
// ----------------------------------------------------------------------------

template <>
struct CompareTypeImpl<Dna, __uint8>
{
    typedef Dna Type;
};

inline void assign(Dna & target, __uint8 c_source)
{
    target.value = TranslateTableByteToDna_<>::VALUE[c_source];
}

template <>
struct CompareTypeImpl<Dna, char>
{
    typedef Dna Type;
};

inline void assign(Dna & target, char c_source)
{
#ifdef __CUDA_ARCH__
    target.value = _translateCharToDna<void>((unsigned char)c_source);
#else
    target.value = TranslateTableCharToDna_<>::VALUE[(unsigned char) c_source];
#endif
}

template <>
struct CompareTypeImpl<Dna, Dna5>
{
    typedef Dna Type;
};

inline void assign(Dna & target, Dna5 const & c_source)
{
    target.value = c_source.value & 0x03;
}

template <>
struct CompareTypeImpl<Dna, Iupac>
{
    typedef Dna Type;
};

inline void assign(Dna & target, Iupac const & source)
{
    target.value = TranslateTableIupacToDna_<>::VALUE[source.value];
}

// ----------------------------------------------------------------------------
// Dna5
// ----------------------------------------------------------------------------

template <>
struct CompareTypeImpl<Dna5, __uint8>
{
    typedef Dna5 Type;
};

inline void assign(Dna5 & target, __uint8 c_source)
{
    target.value = TranslateTableByteToDna5_<>::VALUE[c_source];
}

template <>
struct CompareTypeImpl<Dna5, char>
{
    typedef Dna5 Type;
};

inline void assign(Dna5 & target, char c_source)
{
#ifdef __CUDA_ARCH__
    target.value = _translateCharToDna5<void>((unsigned char)c_source);
#else
    target.value = TranslateTableCharToDna5_<>::VALUE[(unsigned char) c_source];
#endif
}

template <>
struct CompareTypeImpl<Dna5, Iupac>
{
    typedef Dna5 Type;
};

inline void assign(Dna5 & target, Iupac const & source)
{
    target.value = TranslateTableIupacToDna5_<>::VALUE[source.value];
}

template <>
struct CompareTypeImpl<Dna5, Dna>
{
    typedef Dna Type;
};

inline void assign(Dna5 & target, Dna const & c_source)
{
    target.value = c_source.value;
}

// ----------------------------------------------------------------------------
// Rna
// ----------------------------------------------------------------------------

template <>
struct CompareTypeImpl<Rna, __uint8>
{
    typedef Rna Type;
};

inline void assign(Rna & target, __uint8 c_source)
{
    target.value = TranslateTableByteToDna_<>::VALUE[c_source];
}

template <>
struct CompareTypeImpl<Rna, char>
{
    typedef Rna Type;
};

inline void assign(Rna & target, char c_source)
{
    target.value = TranslateTableCharToDna_<>::VALUE[(unsigned char)c_source];
}

template <>
struct CompareTypeImpl<Rna, Rna5>
{
    typedef Rna Type;
};

inline void assign(Rna & target, Rna5 const & c_source)
{
    target.value = c_source.value & 0x03;
}

// ---------------------------------------------------------------------------
// Rna5
// ---------------------------------------------------------------------------

template <>
struct CompareTypeImpl<Rna5, __uint8>
{
    typedef Rna5 Type;
};

inline void assign(Rna5 & target, __uint8 c_source)
{
    target.value = TranslateTableByteToDna5_<>::VALUE[c_source];
}

template <>
struct CompareTypeImpl<Rna5, char>
{
    typedef Rna5 Type;
};

inline void assign(Rna5 & target, char c_source)
{
    target.value = TranslateTableCharToDna5_<>::VALUE[(unsigned char)c_source];
}

template <>
struct CompareTypeImpl<Rna5, Rna>
{
    typedef Dna Type;
};

inline void assign(Rna5 & target, Rna const & c_source)
{
    target.value = c_source.value;
}

// ---------------------------------------------------------------------------
// Iupac
// ---------------------------------------------------------------------------

template <>
struct CompareTypeImpl<Iupac, __uint8>
{
    typedef Iupac Type;
};

inline void assign(Iupac & target, __uint8 c_source)
{
    target.value = TranslateTableByteToIupac_<>::VALUE[c_source];
}

template <>
struct CompareTypeImpl<Iupac, char>
{
    typedef Iupac Type;
};

inline void assign(Iupac & target, char c_source)
{
    target.value = TranslateTableCharToIupac_<>::VALUE[(unsigned char) c_source];
}

inline void assign(Iupac & target, Dna const & source)
{
    target.value = TranslateTableDna5ToIupac_<>::VALUE[source.value];
}

inline void assign(Iupac & target, Dna5 const & source)
{
    target.value = TranslateTableDna5ToIupac_<>::VALUE[source.value];
}

// ---------------------------------------------------------------------------
// Amino Acid
// ---------------------------------------------------------------------------

template <>
struct CompareTypeImpl<AminoAcid, __uint8>
{
    typedef AminoAcid Type;
};

inline void assign(AminoAcid & target, __uint8 c_source)
{
    target.value = TranslateTableByteToAA_<>::VALUE[c_source];
}

template <>
struct CompareTypeImpl<AminoAcid, char>
{
    typedef AminoAcid Type;
};

inline void assign(AminoAcid & target, char c_source)
{
    target.value = TranslateTableCharToAA_<>::VALUE[(unsigned char) c_source];
}

// ---------------------------------------------------------------------------
// DnaQ
// ---------------------------------------------------------------------------

// template <typename TValue, typename TValue2>
// struct CompareTypeImpl<SimpleType<TValue,DnaQ_>, SimpleType<TValue2,Dna_> >
// {
//  typedef SimpleType<TValue2,Dna_> Type;
// };
//
// template <typename TValue, typename TValue2>
// struct CompareTypeImpl<SimpleType<TValue,Dna_>, SimpleType<TValue2,DnaQ_> >
// {
//  typedef SimpleType<TValue,Dna_> Type;
// };

template <>
struct CompareTypeImpl<DnaQ, DnaQ>
{
    typedef Dna Type;
};

template <>
struct CompareTypeImpl<DnaQ, Dna>
{
    typedef Dna Type;
};

inline void assign(DnaQ & target, Dna const & source)
{
    target.value = source.value | (SEQAN_DEFAULT_QUALITY << 2);
}

template <>
struct CompareTypeImpl<Dna, DnaQ>
{
    typedef Dna Type;
};

inline void assign(Dna & target, DnaQ const & source)
{
    target.value = source.value & 3;
}

template <>
struct CompareTypeImpl<DnaQ, Iupac>
{
    typedef Iupac Type;
};

inline void assign(DnaQ & target, Iupac const & source)
{
    assign(target, (Dna) source);
}

template <>
struct CompareTypeImpl<Iupac, DnaQ>
{
    typedef Iupac Type;
};

inline void assign(Iupac & target, DnaQ const & source)
{
    assign(target, (Dna) source);
}

template <>
struct CompareTypeImpl<DnaQ, Dna5>
{
    typedef Dna Type;
};

inline void assign(DnaQ & target, Dna5 const & source)
{
    assign(target, (Dna) source);
}

template <>
struct CompareTypeImpl<DnaQ, __uint8>
{
    typedef Dna Type;
};

inline void assign(DnaQ & target, __uint8 c_source)
{
    assign(target, (Dna) c_source);
}

template <>
struct CompareTypeImpl<DnaQ, char>
{
    typedef Dna Type;
};

inline void assign(DnaQ & target, char c_source)
{
    assign(target, (Dna) c_source);
}

inline void
assign(DnaQ & target, DnaQ const & source)
{
    target.value = source.value;
}

template <typename TSource>
inline void
assign(DnaQ & target, TSource const & source)
{
    target.value = (Dna)source;
}

inline void
assign(__int64 & c_target,
       DnaQ & source)
{
    c_target = Dna(source);
}

inline void
assign(__int64 & c_target,
       DnaQ const & source)
{
    c_target = Dna(source);
}

// __uint64

inline void
assign(__uint64 & c_target,
       DnaQ & source)
{
    c_target = Dna(source);
}

inline void
assign(__uint64 & c_target,
       DnaQ const & source)
{
    c_target = Dna(source);
}

// int

inline void
assign(int & c_target,
       DnaQ & source)
{
    c_target = Dna(source);
}

inline void
assign(int & c_target,
       DnaQ const & source)
{
    c_target = Dna(source);
}

// unsigned int

inline void
assign(unsigned int & c_target,
       DnaQ & source)
{
    c_target = Dna(source);
}

inline void
assign(unsigned int & c_target,
       DnaQ const & source)
{
    c_target = Dna(source);
}

// short

inline void
assign(short & c_target,
       DnaQ & source)
{
    c_target = Dna(source);
}

inline void
assign(short & c_target,
       DnaQ const & source)
{
    c_target = Dna(source);
}

// unsigned short

inline void
assign(unsigned short & c_target,
       DnaQ & source)
{
    c_target = Dna(source);
}

inline void
assign(unsigned short & c_target,
       DnaQ const & source)
{
    c_target = Dna(source);
}

// char

inline void
assign(char & c_target,
       DnaQ & source)
{
    c_target = Dna(source);
}

inline void
assign(char & c_target,
       DnaQ const & source)
{
    c_target = Dna(source);
}

// signed char

inline void
assign(signed char & c_target,
       DnaQ & source)
{
    c_target = Dna(source);
}

inline void
assign(signed char & c_target,
       DnaQ const & source)
{
    c_target = Dna(source);
}

// unsigned char

inline void
assign(unsigned char & c_target,
       DnaQ & source)
{
    c_target = Dna(source);
}

inline void
assign(unsigned char & c_target,
       DnaQ const & source)
{
    c_target = Dna(source);
}

// ---------------------------------------------------------------------------
// Dna5Q
// ---------------------------------------------------------------------------

// template <typename TValue, typename TValue2>
// struct CompareTypeImpl<SimpleType<TValue,Dna5Q_>, SimpleType<TValue2,Dna5_> >
// {
//  typedef SimpleType<TValue2,Dna5_> Type;
// };
//
// template <typename TValue, typename TValue2>
// struct CompareTypeImpl<SimpleType<TValue,Dna5_>, SimpleType<TValue2,Dna5Q_> >
// {
//  typedef SimpleType<TValue,Dna5_> Type;
// };


template <>
struct CompareTypeImpl<Dna5Q, Dna5Q>
{
    typedef Dna5 Type;
};

template <>
struct CompareTypeImpl<DnaQ, Dna5Q>
{
    typedef Dna Type;
};

inline void assign(DnaQ & target, Dna5Q const & source)
{
    // We perform the converstion from DNA5 to DNA5 with qualities by a simple
    // table lookup.  The lookup below is equivalent to the following line:
    //
    // target.value = (source.value == Dna5QValueN_)? 0: source.value;

    static const unsigned table[] = {
          0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
         16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
         32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
         48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
         64,  65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,  78,  79,
         80,  81,  82,  83,  84,  85,  86,  87,  88,  89,  90,  91,  92,  93,  94,  95,
         96,  97,  98,  99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111,
        112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127,
        128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
        144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
        160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
        176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
        192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
        208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
        224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
        240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 0,   0,   0,   0
    };
    target.value = table[source.value];
}

template <>
struct CompareTypeImpl<Dna5Q, DnaQ>
{
    typedef Dna Type;
};

inline void assign(Dna5Q & target, DnaQ const & source)
{
    target.value = source.value;
}


template <>
struct CompareTypeImpl<Dna5, Dna5Q>
{
    typedef Dna5 Type;
};

inline void assign(Dna5 & target, Dna5Q const & source)
{
        SEQAN_CHECKPOINT;;

    // We perform the conversion from DNA5 to DNA5 with qualities by a simple
    // table lookup.  The lookup below is equivalent to the following line:
    //
    // target.value = (source.value == Dna5QValueN_)? 4: source.value & 3;

    static const unsigned table[] = {
        0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 4, 4, 4, 4 // <-- note the 4
    };
    target.value = table[source.value];
}

template <>
struct CompareTypeImpl<Dna5Q, Dna5>
{
    typedef Dna5 Type;
};

inline void assign(Dna5Q & target, Dna5 const & source)
{

    // We perform the conversion from DNA5 with qualities to DNA5 by a simple
    // table lookup.  The lookup below is equivalent to the following line:
    //
    // target.value = (source.value == 4)? Dna5QValueN_ : source.value | (40 << 2);

    static const unsigned table[] = {
        (SEQAN_DEFAULT_QUALITY << 2) + 0, (SEQAN_DEFAULT_QUALITY << 2) + 1,
        (SEQAN_DEFAULT_QUALITY << 2) + 2, (SEQAN_DEFAULT_QUALITY << 2) + 3, Dna5QValueN_
    };
    target.value = table[source.value];
}

template <>
struct CompareTypeImpl<Dna5Q, Dna>
{
    typedef Dna Type;
};

inline void assign(Dna5Q & target, Dna const & source)
{
    assign(target, (DnaQ) source);
}

template <>
struct CompareTypeImpl<Dna, Dna5Q>
{
    typedef Dna Type;
};

inline void assign(Dna & target, Dna5Q const & source)
{
    assign(target, (Dna5)source);
}

template <>
struct CompareTypeImpl<Dna5, DnaQ>
{
    typedef Dna5 Type;
};

inline void assign(Dna5 & target, DnaQ const & source)
{
    assign(target, (Dna5Q)source);
}

template <>
struct CompareTypeImpl<Dna5Q, __uint8>
{
    typedef Dna5 Type;
};

inline void assign(Dna5Q & target, __uint8 c_source)
{
    assign(target, (Dna5)c_source);
}

template <>
struct CompareTypeImpl<Dna5Q, char>
{
    typedef Dna5 Type;
};

inline void assign(Dna5Q & target, char c_source)
{
    assign(target, (Dna5)c_source);
}

template <>
struct CompareTypeImpl<Dna5Q, Iupac>
{
    typedef Iupac Type;
};

inline void assign(Dna5Q & target, Iupac const & source)
{
    assign(target, (Dna5)source);
}

template <>
struct CompareTypeImpl<Iupac, Dna5Q>
{
    typedef Iupac Type;
};

inline void assign(Iupac & target, Dna5Q const & source)
{
    assign(target, (Dna5)source);
}

inline void
assign(Dna5Q & target, Dna5Q const & source)
{
    target.value = source.value;
}

template <typename TSource>
inline void
assign(Dna5Q & target, TSource const & source)
{
    assign(target, (Dna5)source);
}

// __int64

inline void
assign(__int64 & c_target,
       Dna5Q & source)
{
    c_target = Dna5(source);
}

inline void
assign(__int64 & c_target,
       Dna5Q const & source)
{
    c_target = Dna5(source);
}

// __uint64

inline void
assign(__uint64 & c_target,
       Dna5Q & source)
{
    c_target = Dna5(source);
}

inline void
assign(__uint64 & c_target,
       Dna5Q const & source)
{
    c_target = Dna5(source);
}

// int

inline void
assign(int & c_target,
       Dna5Q & source)
{
    c_target = Dna5(source);
}

inline void
assign(int & c_target,
       Dna5Q const & source)
{
    c_target = Dna5(source);
}

// unsigned int

inline void
assign(unsigned int & c_target,
       Dna5Q & source)
{
    c_target = Dna5(source);
}

inline void
assign(unsigned int & c_target,
       Dna5Q const & source)
{
    c_target = Dna5(source);
}


//short

inline void
assign(short & c_target,
       Dna5Q & source)
{
    c_target = Dna5(source);
}

inline void
assign(short & c_target,
       Dna5Q const & source)
{
    c_target = Dna5(source);
}

//unsigned short

inline void
assign(unsigned short & c_target,
       Dna5Q & source)
{
    c_target = Dna5(source);
}

inline void
assign(unsigned short & c_target,
       Dna5Q const & source)
{
    c_target = Dna5(source);
}

// char

inline void
assign(char & c_target,
       Dna5Q & source)
{
    c_target = Dna5(source);
}

inline void
assign(char & c_target,
       Dna5Q const & source)
{
    c_target = Dna5(source);
}

// signed char

inline void
assign(signed char & c_target,
       Dna5Q & source)
{
    c_target = Dna5(source);
}

inline void
assign(signed char & c_target,
       Dna5Q const & source)
{
    c_target = Dna5(source);
}

// unsigned char

inline void
assign(unsigned char & c_target,
       Dna5Q & source)
{
    c_target = Dna5(source);
}

inline void
assign(unsigned char & c_target,
       Dna5Q const & source)
{
    c_target = Dna5(source);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BASIC_BASIC_ALPHABET_RESIDUE_H_
