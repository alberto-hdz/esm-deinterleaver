#pragma once

#include <string>
#include <vector>

#include "esm/pipeline.hpp"

namespace esm {

// One row of the emitter library: the parameter windows a known emitter type
// is expected to show up in. This is the same idea as the threat tables a radar
// warning receiver uses, but every entry in my library is made up.
struct LibraryEntry {
    std::string name;
    std::string role;
    double rf_min_mhz = 0.0;
    double rf_max_mhz = 0.0;
    double pw_min_us = 0.0;
    double pw_max_us = 0.0;
    PriType pri_type = PriType::Constant;
    double pri_min_us = 0.0;   // mean PRI window (for a stagger: frame / levels)
    double pri_max_us = 0.0;
    int stagger_levels = 0;
    double jitter_max = 0.0;
};

// Reads the library CSV (columns by header name, like the PDW reader).
// Throws std::runtime_error with the line number on a bad row.
std::vector<LibraryEntry> read_library_csv(const std::string& path);

struct IdParams {
    double rf_margin_mhz = 10.0;  // how far outside a window a value can be and still partly count
    double pw_margin_frac = 0.10;
    double pw_margin_min_us = 0.1;
    double pri_margin_frac = 0.02;
    double min_score = 0.8;       // below this the emitter is reported as UNKNOWN
    double ambiguous_gap = 0.05;  // runner-up this close to the best match = ambiguous
};

struct IdResult {
    std::string name = "UNKNOWN";
    double score = 0.0;
    std::string runner_up;        // second-best entry, for the ambiguity check
    double runner_up_score = 0.0;
    bool ambiguous = false;
};

// Score one library entry against a measured emitter, 0 to 1.
double match_score(const Emitter& e, const LibraryEntry& lib, const IdParams& p);

// Best library match for an emitter, or UNKNOWN if nothing scores high enough.
IdResult identify(const Emitter& e, const std::vector<LibraryEntry>& library, const IdParams& p);

}  // namespace esm
