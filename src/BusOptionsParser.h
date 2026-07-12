#ifndef KALLISTO_BUS_OPTIONS_PARSER_H
#define KALLISTO_BUS_OPTIONS_PARSER_H

#include <string>
#include <vector>

#include "common.h"

namespace bus_parser {

void ListSingleCellTechnologies();

void ParseOptionsBus(int argc, char **argv, ProgramOptions &opt);

bool ParseTechnology(const std::string &techstr, BUSOptions &busopt,
                     std::vector<std::string> &errorList);

// Apply preset technology (e.g. "10XV3") or a custom -x string to opt.busOptions.
// Sets busopt.{bc,umi,seq,nfiles,paired,...} and writes the technology's
// default strand into strand_out (caller decides whether to apply it).
// Returns true on success; on failure, error messages are pushed onto errorList.
// Does NOT touch fields that depend on batch_mode/bam/tagsequence/etc.
bool ApplyBusPresetTechnology(ProgramOptions &opt,
                              ProgramOptions::StrandType &strand_out,
                              std::vector<std::string> &errorList);

}  // namespace bus_parser

#endif  // KALLISTO_BUS_OPTIONS_PARSER_H
