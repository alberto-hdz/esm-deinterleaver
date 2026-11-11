// esm_cli: for now it just loads a PDW file and prints a summary.
// I'll grow this into the full pipeline (cluster -> PRI analysis -> emitter ID).
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <exception>
#include <string>

#include "esm/pdw.hpp"

int main(int argc, char** argv) {
    std::string pdw_path;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--pdw" && i + 1 < argc) {
            pdw_path = argv[++i];
        } else {
            std::fprintf(stderr, "usage: esm_cli --pdw <pdw.csv>\n");
            return 2;
        }
    }
    if (pdw_path.empty()) {
        std::fprintf(stderr, "usage: esm_cli --pdw <pdw.csv>\n");
        return 2;
    }

    try {
        const auto t0 = std::chrono::steady_clock::now();
        const auto pdws = esm::read_pdw_csv(pdw_path);
        const auto t1 = std::chrono::steady_clock::now();
        if (pdws.empty()) {
            std::printf("no pulses in %s\n", pdw_path.c_str());
            return 0;
        }

        const auto [fmin, fmax] = std::minmax_element(
            pdws.begin(), pdws.end(),
            [](const esm::Pdw& a, const esm::Pdw& b) { return a.freq_mhz < b.freq_mhz; });
        const double span_ms = (pdws.back().toa_us - pdws.front().toa_us) / 1000.0;
        const double load_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

        std::printf("loaded   %zu PDWs from %s in %.1f ms\n", pdws.size(), pdw_path.c_str(), load_ms);
        std::printf("window   %.1f ms (%.0f pulses/s)\n", span_ms,
                    span_ms > 0 ? static_cast<double>(pdws.size()) / (span_ms / 1000.0) : 0.0);
        std::printf("RF range %.1f - %.1f MHz\n", fmin->freq_mhz, fmax->freq_mhz);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
    return 0;
}
