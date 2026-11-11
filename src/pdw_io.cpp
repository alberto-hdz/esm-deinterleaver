#include "esm/pdw.hpp"

#include <array>
#include <cstddef>
#include <fstream>
#include <iomanip>
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
        // strip a trailing \r in case the file came from Windows
        if (!field.empty() && field.back() == '\r') field.pop_back();
        fields.push_back(field);
    }
    return fields;
}

double parse_double(const std::string& s, std::size_t line_no, const std::string& col) {
    std::size_t used = 0;
    double v = 0.0;
    try {
        v = std::stod(s, &used);
    } catch (const std::exception&) {
        used = 0;
    }
    if (used == 0 || used != s.size()) {
        throw std::runtime_error("line " + std::to_string(line_no) + ": bad value '" + s +
                                 "' in column " + col);
    }
    return v;
}

}  // namespace

std::vector<Pdw> read_pdw_csv(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);

    std::string line;
    if (!std::getline(in, line)) throw std::runtime_error(path + " is empty");
    const auto header = split_csv_line(line);

    // Find where each column lives. -1 means not present.
    const std::array<std::string, 6> names = {"toa_us", "freq_mhz", "pw_us",
                                              "amp_db", "aoa_deg", "truth_id"};
    std::array<int, 6> idx{};
    idx.fill(-1);
    for (std::size_t c = 0; c < header.size(); ++c) {
        for (std::size_t k = 0; k < names.size(); ++k) {
            if (header[c] == names[k]) idx[k] = static_cast<int>(c);
        }
    }
    for (std::size_t k = 0; k < 5; ++k) {  // truth_id (k = 5) is optional
        if (idx[k] < 0) throw std::runtime_error(path + ": missing column " + names[k]);
    }

    std::vector<Pdw> pdws;
    std::size_t line_no = 1;
    while (std::getline(in, line)) {
        ++line_no;
        if (line.empty() || line == "\r") continue;
        const auto f = split_csv_line(line);
        if (f.size() != header.size()) {
            throw std::runtime_error("line " + std::to_string(line_no) + ": expected " +
                                     std::to_string(header.size()) + " fields, got " +
                                     std::to_string(f.size()));
        }
        auto get = [&](std::size_t k) {
            return parse_double(f[static_cast<std::size_t>(idx[k])], line_no, names[k]);
        };
        Pdw p;
        p.toa_us = get(0);
        p.freq_mhz = get(1);
        p.pw_us = get(2);
        p.amp_db = get(3);
        p.aoa_deg = get(4);
        if (idx[5] >= 0) p.truth_id = static_cast<int>(get(5));
        pdws.push_back(p);
    }
    return pdws;
}

void write_pdw_csv(const std::string& path, const std::vector<Pdw>& pdws) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path);
    out << "toa_us,freq_mhz,pw_us,amp_db,aoa_deg,truth_id\n";
    out << std::fixed;
    for (const auto& p : pdws) {
        out << std::setprecision(4) << p.toa_us << ',' << std::setprecision(3) << p.freq_mhz
            << ',' << std::setprecision(4) << p.pw_us << ',' << std::setprecision(2)
            << p.amp_db << ',' << p.aoa_deg << ',' << p.truth_id << '\n';
    }
}

}  // namespace esm
