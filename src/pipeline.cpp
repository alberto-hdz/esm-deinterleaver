#include "esm/pipeline.hpp"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <numeric>
#include <vector>

namespace esm {
namespace {

void fill_stats(Emitter& e, const std::vector<Pdw>& pdws) {
    double rf = 0.0;
    double pw = 0.0;
    std::vector<double> aoas;
    for (std::size_t i : e.pulses) {
        rf += pdws[i].freq_mhz;
        pw += pdws[i].pw_us;
        aoas.push_back(pdws[i].aoa_deg);
    }
    const auto n = static_cast<double>(e.pulses.size());
    e.freq_mhz = rf / n;
    e.pw_us = pw / n;
    e.aoa_deg = circular_mean_deg(aoas);
}

}  // namespace

PipelineResult run_pipeline(const std::vector<Pdw>& pdws, const ClusterParams& cp,
                            const PriParams& pp) {
    PipelineResult res;
    res.pulse_label.assign(pdws.size(), -1);
    using clock = std::chrono::steady_clock;
    auto t0 = clock::now();
    const auto labels = cluster_pdws(pdws, cp);
    res.clusters = cluster_count(labels);
    res.cluster_ms = std::chrono::duration<double, std::milli>(clock::now() - t0).count();
    t0 = clock::now();

    for (int c = 0; c < res.clusters; ++c) {
        // Pull this cluster's pulses out and sort them by TOA. I don't assume the
        // input file is sorted, even though my simulator writes it that way.
        std::vector<std::size_t> members;
        for (std::size_t i = 0; i < pdws.size(); ++i) {
            if (labels[i] == c) members.push_back(i);
        }
        std::sort(members.begin(), members.end(),
                  [&](std::size_t a, std::size_t b) { return pdws[a].toa_us < pdws[b].toa_us; });
        std::vector<double> toa(members.size());
        for (std::size_t k = 0; k < members.size(); ++k) toa[k] = pdws[members[k]].toa_us;

        std::vector<bool> assigned(members.size(), false);
        for (const auto& tr : deinterleave(toa, pp)) {
            Emitter e;
            e.id = static_cast<int>(res.emitters.size());
            e.cluster = c;
            e.pri_us = tr.pri_us;
            e.pri_type = tr.type;
            e.jitter = tr.jitter;
            e.stagger_us = tr.stagger_us;
            for (std::size_t k : tr.pulses) {
                e.pulses.push_back(members[k]);  // map back to the original PDW index
                assigned[k] = true;
            }
            std::sort(e.pulses.begin(), e.pulses.end());
            res.emitters.push_back(std::move(e));
        }

        // Whatever is left still looked like one emitter by its parameters,
        // so I report it rather than throwing it away.
        Emitter rest;
        rest.cluster = c;
        for (std::size_t k = 0; k < members.size(); ++k) {
            if (!assigned[k]) rest.pulses.push_back(members[k]);
        }
        if (static_cast<int>(rest.pulses.size()) >= pp.min_pulses) {
            rest.id = static_cast<int>(res.emitters.size());
            res.emitters.push_back(std::move(rest));
        }
    }

    res.pri_ms = std::chrono::duration<double, std::milli>(clock::now() - t0).count();

    for (auto& e : res.emitters) {
        fill_stats(e, pdws);
        for (std::size_t i : e.pulses) res.pulse_label[i] = e.id;
    }
    return res;
}

}  // namespace esm
