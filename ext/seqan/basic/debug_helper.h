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
// Author: Andreas Gogol-Döring <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Helper code for debugging and testing.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_DEBUG_HELPER_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_DEBUG_HELPER_H_

#include <cstdio>
#include <fstream>

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

// TODO(holtgrew): Document, make public.

// compare two files, do not translate linebreaks
inline bool
_compareBinaryFiles(const char * file1, const char * file2)
{
//IOREV see above
    bool ret = false;

    FILE * fl1 = fopen(file1, "rb");
    if (!fl1) return ret;

    FILE * fl2 = fopen(file2, "rb");
    if (!fl2)
    {
        fclose(fl1);

        return ret;
    }

    while (!feof(fl1) && !feof(fl2))
    {
        if (fgetc(fl1) != fgetc(fl2)) goto End;
    }

    ret = feof(fl1) && feof(fl2);

End:
    fclose(fl2);
    fclose(fl1);

    return ret;

}
//____________________________________________________________________________

//one line break is either \r, \n, or \r\n.
//a single line break is skipped.
//the second line break is transformed into \n
inline void
_compareTextFilesReadChar(FILE * fl, char & c, int & num_lb, bool & is_eof)
{
//IOREV see above
    num_lb = 0;
    is_eof = false;

    c = fgetc(fl);
    while ((c == '\r') || (c == '\n'))
    {
        ++num_lb;
        if (c == '\r')
        {
            c = fgetc(fl);
            if (feof(fl)) is_eof = true;
            else
            {
                if (c == '\n')
                {
                    c = fgetc(fl);
                    if (feof(fl)) is_eof = true;
                }
            }
        }
        else if (c == '\n')
        {
            c = fgetc(fl);
            if (feof(fl)) is_eof = true;
        }
    }
}

// compare two files, translate linebreaks
inline bool
_compareTextFiles(const char * file1, const char * file2)
{
//IOREV see above
    FILE * fl1 = fopen(file1, "rb");
    if (!fl1) return false;

    FILE * fl2 = fopen(file2, "rb");
    if (!fl2)
    {
        fclose(fl1);
        return false;
    }

    bool ret = false;

    int num_lb1, num_lb2;
    bool is_eof1, is_eof2;
    char c1, c2;

    while (!feof(fl1) && !feof(fl2))
    {
        _compareTextFilesReadChar(fl1, c1, num_lb1, is_eof1);
        _compareTextFilesReadChar(fl2, c2, num_lb2, is_eof2);

        if (num_lb1 != num_lb2)
        {
            goto End;
        }
        if (is_eof1 ^ is_eof2)
        {
            goto End;
        }
        if (c1 != c2)
        {
            goto End;
        }
    }

    ret = feof(fl1) && feof(fl2);

End:
    fclose(fl2);
    fclose(fl1);

    return ret;

}


// compare two files, translate linebreaks
// more helpful output in case of a difference
inline bool
_compareTextFilesAlt(const char * file1, const char * file2)
{
    std::ifstream fl1(file1);
    std::ifstream fl2(file2);

    if (!fl1.good())
    {
        std::cerr << "Couldn't open file \"" << file1 << "\"" << std::endl;
        return false;
    }
    if (!fl2.good())
    {
        std::cerr << "Couldn't open file \"" << file2 << "\"" << std::endl;
        return false;
    }

    std::string line1;
    std::string line2;

    for (__uint64 lineNo = 1; !fl1.eof() && !fl2.eof(); ++lineNo)
    {
        getline(fl1, line1);
        getline(fl2, line2);

        // Normalize line endings.
        if (line1.size() > 2u && line1[line1.size() - 2] == '\r' && line1[line1.size() - 2] == '\n')
        {
            line1[line1.size() - 2] = '\n';
            line1.erase(line1.size() - 1);
        }
        if (line2.size() > 2u && line2[line2.size() - 2] == '\r' && line2[line2.size() - 2] == '\n')
        {
            line2[line2.size() - 2] = '\n';
            line2.erase(line2.size() - 1);
        }

        if (line1 != line2)
        {
            std::cerr << "The following files are different:" << std::endl;
            std::cerr << '\t' << file1 << std::endl;
            std::cerr << '\t' << file2 << std::endl;
            std::cerr << "Line " << lineNo << " of the text files differ:" << std::endl;
            std::cerr << line1 << std::endl;
            std::cerr << line2 << std::endl;
            return false;
        }
    }
    if (fl1.eof() != fl2.eof())
    {
        std::cerr << "File sizes differ" << std::endl;
        return false;
    }
    return true;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BASIC_DEBUG_HELPER_H_
