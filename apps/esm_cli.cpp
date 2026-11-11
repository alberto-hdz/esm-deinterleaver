// esm_cli: load a PDW file, cluster it, and print what it found.
// I'll keep growing this into the full pipeline (cluster -> PRI analysis -> emitter ID).
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <exception>
#include <map>
#include <string>
#include <vector>

#include "esm/clusterer.hpp"
#include "esm/pdw.hpp"

namespace {

void usage() {
    std::fprintf(stderr,
                 "usage: esm_cli --pdw <pdw.csv> [--no-aoa] [--show-truth]\n"
                 "  --no-aoa      cluster on RF and PW only\n"
                 "  --show-truth  debug: show which simulator emitters landed in each cluster\n");
}

double ms_since(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

}  // namespace

int main(int argc, char** argv) {
    std::string pdw_path;
    bool show_truth = false;
    esm::ClusterParams cparams;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--pdw" && i + 1 < argc) {
            pdw_path = argv[++i];
        } else if (arg == "--no-aoa") {
            cparams.use_aoa = false;
        } else if (arg == "--show-truth") {
            show_truth = true;
        } else {
            usage();
            return 2;
        }
    }
    if (pdw_path.empty()) {
        usage();
        return 2;
    }

    try {
        auto t0 = std::chrono::steady_clock::now();
        const auto pdws = esm::read_pdw_csv(pdw_path);
        const double load_ms = ms_since(t0);
        if (pdws.empty()) {
            std::printf("no pulses in %s\n", pdw_path.c_str());
            return 0;
        }
        const double span_ms = (pdws.back().toa_us - pdws.front().toa_us) / 1000.0;
        std::printf("loaded   %zu PDWs from %s in %.1f ms (window %.1f ms)\n", pdws.size(),
                    pdw_path.c_str(), load_ms, span_ms);

        t0 = std::chrono::steady_clock::now();
        const auto labels = esm::cluster_pdws(pdws, cparams);
        const double cluster_ms = ms_since(t0);
        const int k = esm::cluster_count(labels);
        const auto n_noise = std::count(labels.begin(), labels.end(), -1);
        std::printf("cluster  %d clusters, %ld noise pulses in %.2f ms\n\n", k,
                    static_cast<long>(n_noise), cluster_ms);

        std::printf("%-4s %7s %10s %8s %7s%s\n", "id", "pulses", "RF MHz", "PW us", "AOA",
                    show_truth ? "   truth ids (count)" : "");
        for (int c = 0; c < k; ++c) {
            std::vector<double> aoas;
            double rf = 0.0;
            double pw = 0.0;
            std::map<int, int> truth;
            for (std::size_t i = 0; i < pdws.size(); ++i) {
                if (labels[i] != c) continue;
                rf += pdws[i].freq_mhz;
                pw += pdws[i].pw_us;
                aoas.push_back(pdws[i].aoa_deg);
                ++truth[pdws[i].truth_id];
            }
            const auto n = static_cast<double>(aoas.size());
            std::printf("C%-3d %7zu %10.1f %8.2f %7.1f", c, aoas.size(), rf / n, pw / n,
                        esm::circular_mean_deg(aoas));
            if (show_truth) {
                std::printf("  ");
                for (const auto& [id, cnt] : truth) std::printf(" E%d(%d)", id, cnt);
            }
            std::printf("\n");
        }
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
    return 0;
}
