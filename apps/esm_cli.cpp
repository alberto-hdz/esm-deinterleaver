// esm_cli: runs the whole chain on a PDW file and prints the emitters it found.
// Right now that's clustering + PRI deinterleaving; emitter ID comes in later.
#include <chrono>
#include <cstdio>
#include <exception>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>

#include "esm/pdw.hpp"
#include "esm/pipeline.hpp"

namespace {

void usage() {
    std::fprintf(stderr,
                 "usage: esm_cli --pdw <pdw.csv> [--no-aoa] [--show-truth] [--labels <out.csv>]\n"
                 "  --no-aoa      cluster on RF and PW only\n"
                 "  --show-truth  debug: show which simulator emitters ended up in each result\n"
                 "  --labels      write the emitter id assigned to every pulse (for scoring)\n");
}

double ms_since(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

}  // namespace

int main(int argc, char** argv) {
    std::string pdw_path;
    std::string labels_path;
    bool show_truth = false;
    esm::ClusterParams cparams;
    esm::PriParams pparams;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--pdw" && i + 1 < argc) {
            pdw_path = argv[++i];
        } else if (arg == "--labels" && i + 1 < argc) {
            labels_path = argv[++i];
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
        std::printf("loaded   %zu PDWs from %s in %.1f ms\n", pdws.size(), pdw_path.c_str(),
                    ms_since(t0));

        t0 = std::chrono::steady_clock::now();
        const auto res = esm::run_pipeline(pdws, cparams, pparams);
        long unassigned = 0;
        for (int l : res.pulse_label) unassigned += (l < 0);
        std::printf("pipeline %d clusters -> %zu emitters, %ld pulses unassigned, %.1f ms\n\n",
                    res.clusters, res.emitters.size(), unassigned, ms_since(t0));

        std::printf("%-4s %4s %7s %10s %8s %7s  %-30s%s\n", "id", "clus", "pulses", "RF MHz",
                    "PW us", "AOA", "PRI", show_truth ? " truth ids (count)" : "");
        for (const auto& e : res.emitters) {
            std::printf("T%-3d C%-3d %7zu %10.1f %8.2f %7.1f  ", e.id, e.cluster, e.pulses.size(),
                        e.freq_mhz, e.pw_us, e.aoa_deg);
            // PRI column: what kind of PRI it is and its numbers
            char pri[64];
            if (e.pri_us <= 0.0) {
                std::snprintf(pri, sizeof(pri), "none");
            } else if (e.pri_type == esm::PriType::Jitter) {
                std::snprintf(pri, sizeof(pri), "jitter %.1f +/-%.1f%%", e.pri_us, 100.0 * e.jitter);
            } else if (e.pri_type == esm::PriType::Stagger) {
                std::string s = "stagger";
                for (double v : e.stagger_us) s += " " + std::to_string(static_cast<int>(v + 0.5));
                std::snprintf(pri, sizeof(pri), "%s", s.c_str());
            } else {
                std::snprintf(pri, sizeof(pri), "constant %.2f", e.pri_us);
            }
            std::printf("%-30s", pri);
            if (show_truth) {
                std::map<int, int> truth;
                for (std::size_t i : e.pulses) ++truth[pdws[i].truth_id];
                std::printf(" ");
                for (const auto& [id, cnt] : truth) std::printf(" E%d(%d)", id, cnt);
            }
            std::printf("\n");
        }

        if (!labels_path.empty()) {
            std::ofstream out(labels_path);
            if (!out) throw std::runtime_error("cannot write " + labels_path);
            out << "pulse,emitter\n";
            for (std::size_t i = 0; i < res.pulse_label.size(); ++i) {
                out << i << ',' << res.pulse_label[i] << '\n';
            }
        }
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
    return 0;
}
