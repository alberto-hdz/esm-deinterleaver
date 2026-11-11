#include "esm/emitter_library.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace esm {
namespace {

std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream ss(line);
    std::string field;
    while (std::getline(ss, field, ',')) {
        if (!field.empty() && field.back() == '\r') field.pop_back();
        fields.push_back(field);
    }
    return fields;
}

PriType parse_pri_type(const std::string& s, std::size_t line_no) {
    if (s == "constant") return PriType::Constant;
    if (s == "stagger") return PriType::Stagger;
    if (s == "jitter") return PriType::Jitter;
    throw std::runtime_error("line " + std::to_string(line_no) + ": unknown pri_type '" + s + "'");
}

// 1 inside [lo, hi], falling off linearly to 0 at `margin` outside the window.
// A hard in/out check would throw away an emitter measured 1 MHz past the edge
// of its band, which is exactly the kind of thing measurement noise does.
double window_fit(double x, double lo, double hi, double margin) {
    if (x >= lo && x <= hi) return 1.0;
    const double dist = x < lo ? lo - x : x - hi;
    return std::max(0.0, 1.0 - dist / margin);
}

}  // namespace

std::vector<LibraryEntry> read_library_csv(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::string line;
    if (!std::getline(in, line)) throw std::runtime_error(path + " is empty");

    std::map<std::string, std::size_t> col;
    const auto header = split_csv_line(line);
    for (std::size_t c = 0; c < header.size(); ++c) col[header[c]] = c;
    for (const char* need : {"name", "role", "rf_min_mhz", "rf_max_mhz", "pw_min_us", "pw_max_us",
                             "pri_type", "pri_min_us", "pri_max_us", "stagger_levels", "jitter_max"}) {
        if (!col.count(need)) throw std::runtime_error(path + ": missing column " + need);
    }

    std::vector<LibraryEntry> lib;
    std::size_t line_no = 1;
    while (std::getline(in, line)) {
        ++line_no;
        if (line.empty() || line == "\r") continue;
        const auto f = split_csv_line(line);
        if (f.size() != header.size()) {
            throw std::runtime_error("line " + std::to_string(line_no) + ": wrong number of fields");
        }
        auto num = [&](const char* name) {
            try {
                return std::stod(f[col[name]]);
            } catch (const std::exception&) {
                throw std::runtime_error("line " + std::to_string(line_no) + ": bad " + name);
            }
        };
        LibraryEntry e;
        e.name = f[col["name"]];
        e.role = f[col["role"]];
        e.rf_min_mhz = num("rf_min_mhz");
        e.rf_max_mhz = num("rf_max_mhz");
        e.pw_min_us = num("pw_min_us");
        e.pw_max_us = num("pw_max_us");
        e.pri_type = parse_pri_type(f[col["pri_type"]], line_no);
        e.pri_min_us = num("pri_min_us");
        e.pri_max_us = num("pri_max_us");
        e.stagger_levels = static_cast<int>(num("stagger_levels"));
        e.jitter_max = num("jitter_max");
        if (e.rf_min_mhz > e.rf_max_mhz || e.pw_min_us > e.pw_max_us || e.pri_min_us > e.pri_max_us) {
            throw std::runtime_error("line " + std::to_string(line_no) + ": a min is bigger than its max");
        }
        lib.push_back(e);
    }
    return lib;
}

double match_score(const Emitter& e, const LibraryEntry& lib, const IdParams& p) {
    const double rf = window_fit(e.freq_mhz, lib.rf_min_mhz, lib.rf_max_mhz, p.rf_margin_mhz);
    const double pw_margin = std::max(p.pw_margin_min_us, p.pw_margin_frac * lib.pw_max_us);
    const double pw = window_fit(e.pw_us, lib.pw_min_us, lib.pw_max_us, pw_margin);

    // The library stores mean PRI, so for a stagger I compare frame / levels.
    double pri = 0.5;        // no PRI measured: can't confirm or rule out
    double type_factor = 0.6;
    if (e.pri_us > 0.0) {
        const bool stag = e.pri_type == PriType::Stagger && !e.stagger_us.empty();
        const double mean_pri = stag ? e.pri_us / static_cast<double>(e.stagger_us.size()) : e.pri_us;
        pri = window_fit(mean_pri, lib.pri_min_us, lib.pri_max_us, p.pri_margin_frac * lib.pri_max_us);

        if (e.pri_type == lib.pri_type) {
            type_factor = 1.0;
            // a 3-slot stagger shouldn't fully match a 4-slot library entry
            if (stag && lib.stagger_levels > 0 &&
                static_cast<int>(e.stagger_us.size()) != lib.stagger_levels) {
                type_factor = 0.7;
            }
        } else if (e.pri_type != PriType::Stagger && lib.pri_type != PriType::Stagger) {
            type_factor = 0.6;  // constant vs jitter: 1% jitter is easy to measure as constant
        } else {
            type_factor = 0.3;  // stagger vs anything else is a real mismatch
        }
    }
    return type_factor * (rf + pw + pri) / 3.0;
}

IdResult identify(const Emitter& e, const std::vector<LibraryEntry>& library, const IdParams& p) {
    IdResult r;
    double best = -1.0;
    double second = -1.0;
    std::string best_name;
    std::string second_name;
    for (const auto& lib : library) {
        const double s = match_score(e, lib, p);
        if (s > best) {
            second = best;
            second_name = best_name;
            best = s;
            best_name = lib.name;
        } else if (s > second) {
            second = s;
            second_name = lib.name;
        }
    }
    r.score = std::max(best, 0.0);
    r.runner_up = second_name;
    r.runner_up_score = std::max(second, 0.0);
    if (best >= p.min_score) {
        r.name = best_name;
        r.ambiguous = second >= p.min_score && best - second < p.ambiguous_gap;
    }
    return r;
}

}  // namespace esm
