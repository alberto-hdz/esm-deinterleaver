#pragma once

#include <string>
#include <vector>

namespace esm {

// One Pulse Descriptor Word: everything the receiver measured about one pulse.
struct Pdw {
    double toa_us = 0.0;    // time of arrival (microseconds)
    double freq_mhz = 0.0;  // carrier frequency (MHz)
    double pw_us = 0.0;     // pulse width (microseconds)
    double amp_db = 0.0;    // amplitude (dB)
    double aoa_deg = 0.0;   // angle of arrival, 0-360 (degrees)
    int truth_id = -1;      // simulator label, only used for scoring, never by the algorithms
};

// I read columns by header name so the column order in the file doesn't matter.
// truth_id is optional (a real receiver wouldn't have it). Throws std::runtime_error
// with the line number if the file is missing, a column is missing, or a value is bad.
std::vector<Pdw> read_pdw_csv(const std::string& path);

void write_pdw_csv(const std::string& path, const std::vector<Pdw>& pdws);

}  // namespace esm
