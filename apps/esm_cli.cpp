// esm_cli: runs the whole chain and prints the emitters it found:
// (I/Q pulse detection) -> clustering -> PRI deinterleaving -> (library identification).
// The input is either a PDW file or a raw I/Q capture.
#include <chrono>
#include <cmath>
#include <cstdio>
#include <exception>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "esm/emitter_library.hpp"
#include "esm/pdw.hpp"
#include "esm/pipeline.hpp"
#include "esm/pulse_detector.hpp"

namespace {

void usage() {
    std::fprintf(stderr,
                 "usage: esm_cli --pdw <pdw.csv> [options]\n"
                 "       esm_cli --iq <iq.c64> --fs-mhz <rate> --fc-mhz <center> [options]\n"
                 "options: [--library <lib.csv>] [--no-aoa] [--show-truth] [--labels <out.csv>]\n"
                 "         [--emitters <out.csv>] [--pdw-out <out.csv>]\n"
                 "  --iq          raw complex64 capture; pulses are detected first (implies --no-aoa)\n"
                 "  --pdw-out     write the PDWs the detector measured from the I/Q capture\n"
                 "  --library     identify each emitter against this library\n"
                 "  --no-aoa      cluster on RF and PW only\n"
                 "  --show-truth  debug: show which simulator emitters ended up in each result\n"
                 "  --labels      write the emitter id assigned to every pulse (for scoring)\n"
                 "  --emitters    write one row per emitter: parameters and library match\n");
}

double ms_since(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

}  // namespace

int main(int argc, char** argv) {
    std::string pdw_path;
    std::string iq_path;
    std::string pdw_out_path;
    esm::DetectorParams dparams;
    std::string labels_path;
    std::string library_path;
    std::string emitters_path;
    bool show_truth = false;
    esm::ClusterParams cparams;
    esm::PriParams pparams;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--pdw" && i + 1 < argc) {
            pdw_path = argv[++i];
        } else if (arg == "--iq" && i + 1 < argc) {
            iq_path = argv[++i];
        } else if (arg == "--fs-mhz" && i + 1 < argc) {
            dparams.fs_hz = std::stod(argv[++i]) * 1e6;
        } else if (arg == "--fc-mhz" && i + 1 < argc) {
            dparams.fc_mhz = std::stod(argv[++i]);
        } else if (arg == "--pdw-out" && i + 1 < argc) {
            pdw_out_path = argv[++i];
        } else if (arg == "--labels" && i + 1 < argc) {
            labels_path = argv[++i];
        } else if (arg == "--library" && i + 1 < argc) {
            library_path = argv[++i];
        } else if (arg == "--emitters" && i + 1 < argc) {
            emitters_path = argv[++i];
        } else if (arg == "--no-aoa") {
            cparams.use_aoa = false;
        } else if (arg == "--show-truth") {
            show_truth = true;
        } else {
            usage();
            return 2;
        }
    }
    if (pdw_path.empty() == iq_path.empty()) {  // need exactly one input
        usage();
        return 2;
    }

    try {
        auto t0 = std::chrono::steady_clock::now();
        std::vector<esm::Pdw> pdws;
        if (!iq_path.empty()) {
            const auto iq = esm::read_iq_file(iq_path);
            esm::DetectorStats st;
            pdws = esm::detect_pulses(iq, dparams, &st);
            cparams.use_aoa = false;  // one antenna, no AOA
            std::printf("detected %zu pulses in %zu samples (%.1f ms) from %s in %.1f ms, "
                        "noise floor %.1f dB\n", pdws.size(), iq.size(),
                        static_cast<double>(iq.size()) / dparams.fs_hz * 1e3, iq_path.c_str(),
                        ms_since(t0), 10.0 * std::log10(st.noise_power));
            if (!pdw_out_path.empty()) esm::write_pdw_csv(pdw_out_path, pdws);
        } else {
            pdws = esm::read_pdw_csv(pdw_path);
            std::printf("loaded   %zu PDWs from %s in %.1f ms\n", pdws.size(), pdw_path.c_str(),
                        ms_since(t0));
        }

        t0 = std::chrono::steady_clock::now();
        const auto res = esm::run_pipeline(pdws, cparams, pparams);
        long unassigned = 0;
        for (int l : res.pulse_label) unassigned += (l < 0);
        std::printf("pipeline %d clusters -> %zu emitters, %ld pulses unassigned, %.1f ms\n",
                    res.clusters, res.emitters.size(), unassigned, ms_since(t0));

        // Identification is optional: without a library every emitter is just UNKNOWN.
        std::vector<esm::IdResult> ids(res.emitters.size());
        if (!library_path.empty()) {
            const auto library = esm::read_library_csv(library_path);
            const esm::IdParams idp;
            for (std::size_t k = 0; k < res.emitters.size(); ++k) {
                ids[k] = esm::identify(res.emitters[k], library, idp);
            }
            std::printf("library  %zu entries from %s\n", library.size(), library_path.c_str());
        }
        std::printf("\n");

        const bool show_id = !library_path.empty();
        std::printf("%-4s %4s %7s %10s %8s %7s  %-30s%s%s\n", "id", "clus", "pulses", "RF MHz",
                    "PW us", "AOA", "PRI", show_id ? "ID (score)          " : "",
                    show_truth ? " truth ids (count)" : "");
        for (std::size_t k = 0; k < res.emitters.size(); ++k) {
            const auto& e = res.emitters[k];
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
            if (show_id) {
                // "?" marks an ambiguous match: the runner-up scored almost as high
                char id[64];
                std::snprintf(id, sizeof(id), "%s%s (%.2f)", ids[k].name.c_str(),
                              ids[k].ambiguous ? "?" : "", ids[k].score);
                std::printf("%-20s", id);
            }
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

        if (!emitters_path.empty()) {
            std::ofstream out(emitters_path);
            if (!out) throw std::runtime_error("cannot write " + emitters_path);
            out << "emitter,pulses,freq_mhz,pw_us,aoa_deg,pri_type,pri_us,jitter,name,score,ambiguous\n";
            for (std::size_t k = 0; k < res.emitters.size(); ++k) {
                const auto& e = res.emitters[k];
                out << e.id << ',' << e.pulses.size() << ',' << e.freq_mhz << ',' << e.pw_us << ','
                    << e.aoa_deg << ',' << (e.pri_us > 0.0 ? esm::pri_type_name(e.pri_type) : "none")
                    << ',' << e.pri_us << ',' << e.jitter << ',' << ids[k].name << ','
                    << ids[k].score << ',' << (ids[k].ambiguous ? 1 : 0) << '\n';
            }
        }
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
    return 0;
}
