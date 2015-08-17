#ifndef KALLISTO_PLAINTEXT_WRITER_H
#define KALLISTO_PLAINTEXT_WRITER_H

#include <assert.h>
#include <cstdlib>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

void plaintext_writer(
    const std::string& out_name,
    const std::vector<std::string>& targ_ids,
    const std::vector<double>& alpha,
    const std::vector<double>& eff_lens,
    const std::vector<int>& lens
    );

std::string to_json(const std::string& id, const std::string& val, bool quote,
    bool comma = true, int level = 1);

void plaintext_aux(
    const std::string& out_fname,
    const std::string& n_targs,
    const std::string& n_bootstrap,
    const std::string& n_processed,
    const std::string& version,
    const std::string& index_v,
    const std::string& start_time,
    const std::string& call);

#endif
