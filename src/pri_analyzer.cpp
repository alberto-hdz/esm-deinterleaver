#include "esm/pri_analyzer.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace esm {
namespace {

std::size_t bin_count(const PriParams& p) {
    return static_cast<std::size_t>(std::ceil(p.max_pri_us / p.bin_us));
}

double bin_center(std::size_t k, const PriParams& p) {
    return (static_cast<double>(k) + 0.5) * p.bin_us;
}

// Refine a bin center into an actual PRI by averaging the raw differences near it.
// Bins are 1 us wide, so the center alone could be off by half a microsecond, and
// that error adds up fast over a 100-pulse train.
double refine_pri(const std::vector<double>& toa, const std::vector<bool>& available, int level,
                  double tau, const PriParams& p) {
    std::vector<std::size_t> idx;
    for (std::size_t i = 0; i < toa.size(); ++i) {
        if (available[i]) idx.push_back(i);
    }
    const auto c = static_cast<std::size_t>(level);
    double sum = 0.0;
    int n = 0;
    for (std::size_t i = 0; i + c < idx.size(); ++i) {
        const double d = toa[idx[i + c]] - toa[idx[i]];
        if (std::abs(d - tau) <= 1.5 * p.bin_us) {
            sum += d;
            ++n;
        }
    }
    return n > 0 ? sum / n : tau;
}

}  // namespace

std::vector<double> difference_histogram(const std::vector<double>& toa,
                                         const std::vector<bool>& available, int level,
                                         const PriParams& p) {
    const std::size_t nb = bin_count(p);
    std::vector<double> raw(nb, 0.0);

    // Differences are taken between the remaining pulses only, so I collect them first.
    // Pulses already assigned to a train shouldn't create differences anymore.
    std::vector<double> t;
    for (std::size_t i = 0; i < toa.size(); ++i) {
        if (available[i]) t.push_back(toa[i]);
    }
    const auto c = static_cast<std::size_t>(level);
    for (std::size_t i = 0; i + c < t.size(); ++i) {
        const double d = t[i + c] - t[i];
        if (d < 0.0 || d >= p.max_pri_us) continue;
        ++raw[static_cast<std::size_t>(d / p.bin_us)];
    }

    // My Python prototype found a 310 us PRI split 70/44 across two bins,
    // which put both halves under the threshold. Summing each bin with its
    // neighbors fixes that.
    std::vector<double> h(nb, 0.0);
    for (std::size_t k = 0; k < nb; ++k) {
        h[k] = raw[k];
        if (k > 0) h[k] += raw[k - 1];
        if (k + 1 < nb) h[k] += raw[k + 1];
    }
    return h;
}

double sdif_threshold(std::size_t k, int level, std::size_t n, const PriParams& p) {
    // T(tau) = gamma * (N - c) * exp(-tau / (theta * B))
    // N - c is how many differences exist at this level. The exponential lowers
    // the bar for long PRIs, since those emitters give fewer pulses per window.
    const double nb = static_cast<double>(bin_count(p));
    const double count = static_cast<double>(n) - level;
    return p.gamma * std::max(count, 0.0) * std::exp(-bin_center(k, p) / p.bin_us / (p.theta * nb));
}

std::vector<PulseTrain> sequence_search(const std::vector<double>& toa,
                                        const std::vector<bool>& available, double pri_us,
                                        const PriParams& p) {
    std::vector<PulseTrain> trains;
    std::vector<bool> taken(toa.size(), false);  // so two trains in this call can't share a pulse
    if (toa.empty()) return trains;

    for (std::size_t start = 0; start < toa.size(); ++start) {
        if (!available[start] || taken[start]) continue;

        PulseTrain tr;
        tr.pri_us = pri_us;
        tr.pulses.push_back(start);
        double expected = toa[start] + pri_us;
        int misses = 0;
        int slots = 1;  // how many pulses this train should have had so far

        while (expected <= toa.back() + p.tol_us && misses <= p.max_misses) {
            // Search only the pulses inside the tolerance window, using binary search
            // since toa is sorted. I take the closest one if there's more than one.
            auto it = std::lower_bound(toa.begin(), toa.end(), expected - p.tol_us);
            std::size_t best = toa.size();
            double best_err = p.tol_us;
            for (; it != toa.end() && *it <= expected + p.tol_us; ++it) {
                const auto j = static_cast<std::size_t>(it - toa.begin());
                const double err = std::abs(*it - expected);
                if (available[j] && !taken[j] && err <= best_err) {
                    best = j;
                    best_err = err;
                }
            }
            ++slots;
            if (best < toa.size()) {
                tr.pulses.push_back(best);
                // Re-anchor on the pulse I actually found instead of adding PRIs to the
                // start time. Otherwise a tiny PRI error grows into a big one by pulse 100.
                expected = toa[best] + pri_us;
                misses = 0;
            } else {
                expected += pri_us;  // pulse was dropped, skip ahead and keep looking
                ++misses;
            }
        }

        // The misses at the very end are just where the train stopped, not dropped pulses.
        slots -= misses;
        const double fill = static_cast<double>(tr.pulses.size()) / slots;

        // Without the fill check (and with min_hits at 8), a jittered emitter gave me short
        // fake trains: allowing misses lets the search hop over off-time pulses until
        // enough happen to line up. A real train with 30% pulse loss still fills ~70%.
        if (static_cast<int>(tr.pulses.size()) >= p.min_hits && fill >= p.min_fill) {
            for (std::size_t j : tr.pulses) taken[j] = true;
            trains.push_back(std::move(tr));
        }
    }

    // If more than max_misses pulses in a row get dropped, one emitter's train breaks
    // into two pieces. I glue a piece onto an earlier one when the gap between them is
    // a whole number of PRIs, meaning it's the same emitter and same timing.
    // (A unit test with 10% drop caught this: it found 3 trains for 2 emitters.)
    std::vector<PulseTrain> merged;
    for (auto& tr : trains) {  // trains are already in order of their first pulse
        bool joined = false;
        for (auto& m : merged) {
            const double gap = toa[tr.pulses.front()] - toa[m.pulses.back()];
            if (gap <= 0.0) continue;
            const double k = std::round(gap / pri_us);
            if (std::abs(gap - k * pri_us) <= p.tol_us) {
                m.pulses.insert(m.pulses.end(), tr.pulses.begin(), tr.pulses.end());
                joined = true;
                break;
            }
        }
        if (!joined) merged.push_back(std::move(tr));
    }
    return merged;
}

std::vector<PulseTrain> deinterleave(const std::vector<double>& toa, const PriParams& p) {
    std::vector<PulseTrain> result;
    std::vector<bool> available(toa.size(), true);
    std::size_t remaining = toa.size();

    // The SDIF loop: look at gaps between neighbors (c = 1), then every second pulse
    // (c = 2), and so on. As soon as a PRI pulls out a train, I remove those pulses and
    // start over at c = 1, because the histogram is much cleaner without them.
    bool extracted = true;
    while (extracted && remaining >= static_cast<std::size_t>(p.min_pulses)) {
        extracted = false;
        for (int c = 1; c <= p.max_level && !extracted; ++c) {
            const auto h = difference_histogram(toa, available, c, p);

            // I go from the smallest tau up and try each bin over the threshold.
            // Taking the tallest peak instead is a trap: in the prototype it picked
            // 620 us, which is 2x the real 310 us PRI.
            for (std::size_t k = 0; k < h.size() && !extracted; ++k) {
                if (h[k] <= sdif_threshold(k, c, remaining, p)) continue;
                // only use local peaks, otherwise the neighbor of a peak gets tried too
                if (k > 0 && h[k - 1] > h[k]) continue;
                if (k + 1 < h.size() && h[k + 1] > h[k]) continue;

                const double pri = refine_pri(toa, available, c, bin_center(k, p), p);
                auto trains = sequence_search(toa, available, pri, p);
                for (auto& tr : trains) {
                    for (std::size_t j : tr.pulses) available[j] = false;
                    remaining -= tr.pulses.size();
                    result.push_back(std::move(tr));
                    extracted = true;
                }
            }
        }
    }
    return result;
}

}  // namespace esm
