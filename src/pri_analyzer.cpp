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

// Dropped pulses are random, so the holes in a real train shouldn't repeat on a
// schedule. A 300/400/500 us stagger searched at 800 us hits 2 slots out of every
// 3, so the fill check alone can't reject it, but its holes are always the same
// slot mod 3. I put every pulse on the PRI grid, list the empty slots, and for
// m = 2..6 count which slot-mod-m they fall in. If one bin is way over the 1/m a
// random process would give (z-score above 4), the train isn't real.
bool has_periodic_holes(const std::vector<double>& toa, const std::vector<std::size_t>& pulses,
                        double pri_us) {
    if (pulses.size() < 2) return false;
    std::vector<bool> present;
    for (std::size_t j : pulses) {
        const auto slot = static_cast<std::size_t>(std::llround((toa[j] - toa[pulses.front()]) / pri_us));
        if (slot >= present.size()) present.resize(slot + 1, false);
        present[slot] = true;
    }
    std::vector<int> holes;
    for (std::size_t k = 0; k < present.size(); ++k) {
        if (!present[k]) holes.push_back(static_cast<int>(k));
    }
    if (holes.size() < 6) return false;

    const auto n = static_cast<double>(holes.size());
    for (int m = 2; m <= 6; ++m) {
        std::vector<int> count(static_cast<std::size_t>(m), 0);
        for (int h : holes) ++count[static_cast<std::size_t>(h % m)];
        const double top = *std::max_element(count.begin(), count.end());
        const double p = 1.0 / m;
        const double z = (top - n * p) / std::sqrt(n * p * (1.0 - p));
        if (z > 4.0) return true;
    }
    return false;
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
    // T(tau) = gamma * (N - c) * exp(-tau / decay)
    // N - c is how many differences exist at this level. The exponential lowers
    // the bar for long PRIs, since those emitters give fewer pulses per window.
    // The paper writes the decay as theta * (number of bins); I keep it in
    // microseconds so raising max_pri_us doesn't change the threshold shape.
    const double count = static_cast<double>(n) - level;
    return p.gamma * std::max(count, 0.0) * std::exp(-bin_center(k, p) / p.decay_us);
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
            // For a jittered train each missed pulse adds another interval of
            // uncertainty, so the jitter pass lets the window grow after a miss.
            const double tol = p.tol_us * (1.0 + p.tol_growth * misses);
            auto it = std::lower_bound(toa.begin(), toa.end(), expected - tol);
            std::size_t best = toa.size();
            double best_err = tol;
            for (; it != toa.end() && *it <= expected + tol; ++it) {
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
    // into pieces. (A unit test with 10% drop caught this: it found 3 trains for 2
    // emitters.) I glue a piece onto another train when it sits a whole number of PRIs
    // from that train. Two things I learned the hard way: compare against the other
    // train's pulses nearest in time, since a 0.01 us PRI error over 140 PRIs is already
    // 1.4 us; and use medians, since either train can end or start with a stray pulse
    // from another emitter that the tolerance let in.
    auto wrap = [&](double gap) { return std::abs(gap - std::round(gap / pri_us) * pri_us); };
    auto offset_from = [&](const PulseTrain& m, std::size_t j) {
        auto it = std::lower_bound(m.pulses.begin(), m.pulses.end(), toa[j],
                                   [&](std::size_t a, double t) { return toa[a] < t; });
        const auto pos = static_cast<long>(it - m.pulses.begin());
        const long n = static_cast<long>(m.pulses.size());
        std::vector<double> off;
        for (long q = std::max(0L, pos - 3); q < std::min(n, pos + 3); ++q) {
            off.push_back(wrap(toa[j] - toa[m.pulses[static_cast<std::size_t>(q)]]));
        }
        std::nth_element(off.begin(), off.begin() + static_cast<long>(off.size() / 2), off.end());
        return off[off.size() / 2];
    };
    std::vector<PulseTrain> merged;
    for (auto& tr : trains) {
        bool joined = false;
        for (auto& m : merged) {
            std::vector<double> off;
            for (std::size_t j : tr.pulses) off.push_back(offset_from(m, j));
            std::nth_element(off.begin(), off.begin() + static_cast<long>(off.size() / 2), off.end());
            if (off[off.size() / 2] <= p.tol_us) {
                m.pulses.insert(m.pulses.end(), tr.pulses.begin(), tr.pulses.end());
                std::sort(m.pulses.begin(), m.pulses.end());  // keep time order for offset_from
                joined = true;
                break;
            }
        }
        if (!joined) merged.push_back(std::move(tr));
    }

    // If any train at this PRI has holes on a schedule, the PRI itself is a stagger
    // artifact. Short pieces at the same PRI can slip past the hole test (too few holes
    // to be sure), so I throw out everything found at this PRI.
    for (auto& m : merged) {
        std::sort(m.pulses.begin(), m.pulses.end());
        if (has_periodic_holes(toa, m.pulses, pri_us)) return {};
    }
    return merged;
}

PulseTrain merge_stagger(const std::vector<double>& toa, std::vector<PulseTrain> trains,
                         double frame_us) {
    // Each train is one "slot" of the stagger pattern, so all of them repeat at the
    // frame period and only their phase is different. I find each slot's phase
    // inside the frame with a circular mean (a phase near 0 and one near the frame
    // length are neighbors), sort them, and the gaps between phases are the intervals.
    constexpr double kTwoPi = 6.28318530717958647692;
    std::vector<double> phases;
    PulseTrain out;
    out.type = PriType::Stagger;
    out.pri_us = frame_us;
    for (const auto& tr : trains) {
        double s = 0.0;
        double c = 0.0;
        for (std::size_t j : tr.pulses) {
            const double ang = kTwoPi * std::fmod(toa[j], frame_us) / frame_us;
            s += std::sin(ang);
            c += std::cos(ang);
        }
        double ph = std::atan2(s, c) / kTwoPi * frame_us;
        if (ph < 0.0) ph += frame_us;
        phases.push_back(ph);
        out.pulses.insert(out.pulses.end(), tr.pulses.begin(), tr.pulses.end());
    }
    std::sort(phases.begin(), phases.end());
    for (std::size_t i = 0; i < phases.size(); ++i) {
        const double next = (i + 1 < phases.size()) ? phases[i + 1] : phases[0] + frame_us;
        out.stagger_us.push_back(next - phases[i]);
    }
    std::sort(out.pulses.begin(), out.pulses.end());
    return out;
}

JitterEstimate estimate_jittered_pri(const std::vector<double>& toa,
                                     const std::vector<bool>& available, const PriParams& p) {
    std::vector<double> t;
    for (std::size_t i = 0; i < toa.size(); ++i) {
        if (available[i]) t.push_back(toa[i]);
    }
    if (t.size() < static_cast<std::size_t>(p.min_pulses)) return {};

    // Gaps between neighboring pulses. For one jittered emitter most of these land
    // in [PRI(1 - j), PRI(1 + j)]; dropped pulses add a few at 2x and 3x.
    std::vector<double> d;
    for (std::size_t i = 1; i < t.size(); ++i) d.push_back(t[i] - t[i - 1]);
    std::sort(d.begin(), d.end());

    // Count how many gaps fall within +/- max_jitter of tau. The window scales
    // with tau, which a fixed-width histogram bin can't do. I step tau up by 0.5%
    // at a time and keep the tau that catches the most gaps.
    const double w = p.max_jitter;
    auto in_window = [&](double tau) {
        auto lo = std::lower_bound(d.begin(), d.end(), tau * (1.0 - w));
        auto hi = std::upper_bound(d.begin(), d.end(), tau * (1.0 + w));
        return std::make_pair(lo, hi);
    };
    double best_tau = 0.0;
    long best_count = 0;
    for (double tau = std::max(d.front(), p.bin_us); tau <= std::min(d.back(), p.max_pri_us);
         tau *= 1.005) {
        auto [lo, hi] = in_window(tau);
        const long count = hi - lo;
        if (count > best_count) {
            best_count = count;
            best_tau = tau;
        }
    }
    if (best_count < static_cast<long>(p.jitter_share * static_cast<double>(d.size()))) return {};

    // Re-center on the mean of the gaps in the window (twice, since the first
    // window can be off-center), then measure the spread with the 2nd and 98th
    // percentiles so one stray gap doesn't blow up the jitter estimate.
    double pri = best_tau;
    for (int pass = 0; pass < 2; ++pass) {
        auto [lo, hi] = in_window(pri);
        double sum = 0.0;
        for (auto it = lo; it != hi; ++it) sum += *it;
        pri = sum / static_cast<double>(hi - lo);
    }
    auto [lo, hi] = in_window(pri);
    const auto n = static_cast<std::size_t>(hi - lo);
    const double q02 = *(lo + static_cast<long>(0.02 * static_cast<double>(n - 1)));
    const double q98 = *(lo + static_cast<long>(0.98 * static_cast<double>(n - 1)));
    return {pri, (q98 - q02) / (2.0 * pri)};
}

std::vector<PulseTrain> deinterleave(const std::vector<double>& toa, const PriParams& p) {
    std::vector<PulseTrain> result;
    std::vector<bool> available(toa.size(), true);
    std::size_t remaining = toa.size();

    auto take = [&](PulseTrain tr) {
        for (std::size_t j : tr.pulses) available[j] = false;
        remaining -= tr.pulses.size();
        result.push_back(std::move(tr));
    };

    // Pass 1, the SDIF loop: look at gaps between neighbors (c = 1), then every
    // second pulse (c = 2), and so on. As soon as a PRI pulls out a train, I remove
    // those pulses and start over at c = 1, because the histogram is much cleaner
    // without them.
    bool extracted = true;
    while (extracted && remaining >= static_cast<std::size_t>(p.min_pulses)) {
        extracted = false;
        for (int c = 1; c <= p.max_level && !extracted; ++c) {
            const auto h = difference_histogram(toa, available, c, p);

            // Local peaks over the threshold, smallest tau first. Taking the tallest
            // peak instead is a trap: in the prototype it picked 620 us, which is 2x
            // the real 310 us PRI.
            std::vector<std::size_t> peaks;
            for (std::size_t k = 0; k < h.size(); ++k) {
                if (h[k] <= sdif_threshold(k, c, remaining, p)) continue;
                // ">=" on the left so a flat-topped peak only counts once
                if (k > 0 && h[k - 1] >= h[k]) continue;
                if (k + 1 < h.size() && h[k + 1] > h[k]) continue;
                peaks.push_back(k);
            }

            // Rule from the SDIF paper: at c = 1 only trust a peak if it's the only one.
            // I skipped this at first and a 300/400/500 us stagger got pulled out as a
            // "constant 400 us" train: from a 400 gap, skipping one slot lands on
            // 500 + 300 = 800 = 2 x 400, so the search hit 2 of every 3 pulses.
            if (c == 1 && peaks.size() > 1) continue;

            for (std::size_t k : peaks) {
                if (extracted) break;
                const double pri = refine_pri(toa, available, c, bin_center(k, p), p);
                auto trains = sequence_search(toa, available, pri, p);
                if (trains.empty()) continue;
                extracted = true;

                // A staggered emitter shows up here as several trains at its frame
                // period, one per slot in the pattern. My first version reported each
                // slot as its own emitter. Same PRI + same cluster = one emitter.
                // (Two constant emitters with the exact same PRI would also get merged.
                // From timing alone those look identical to a 2-level stagger.)
                if (trains.size() >= 2) {
                    auto st = merge_stagger(toa, std::move(trains), pri);
                    // If the slots are evenly spaced it was never a stagger. It's one
                    // constant train found at a multiple of its PRI (e.g. 620 us for a
                    // 310 us emitter, which shows up as two trains half a frame apart).
                    const auto [lo, hi] = std::minmax_element(st.stagger_us.begin(), st.stagger_us.end());
                    if (*hi - *lo <= 2.0 * p.tol_us) {
                        st.type = PriType::Constant;
                        st.pri_us = pri / static_cast<double>(st.stagger_us.size());
                        st.stagger_us.clear();
                    }
                    take(std::move(st));
                } else {
                    take(std::move(trains[0]));
                }
            }
        }
    }

    // Pass 2, jitter: SDIF with a 1.5 us tolerance can't follow a PRI that moves
    // around by several percent, so those pulses are still here. If most of the
    // leftover gaps agree on one PRI, I run the sequence search again with a
    // tolerance sized to the measured jitter.
    while (remaining >= static_cast<std::size_t>(p.min_pulses)) {
        const auto est = estimate_jittered_pri(toa, available, p);
        if (est.pri_us <= 0.0) break;

        PriParams q = p;
        q.tol_us = p.tol_us + 1.2 * est.jitter * est.pri_us;
        q.tol_growth = 1.0;
        auto trains = sequence_search(toa, available, est.pri_us, q);
        if (trains.empty()) break;

        // Breaks in a jittered train don't line up on whole PRIs, so the gap merge
        // in sequence_search can't rejoin them. Everything found at this PRI goes
        // into one emitter.
        PulseTrain tr;
        tr.type = est.jitter >= 0.01 ? PriType::Jitter : PriType::Constant;
        tr.pri_us = est.pri_us;
        tr.jitter = tr.type == PriType::Jitter ? est.jitter : 0.0;
        for (const auto& t : trains) tr.pulses.insert(tr.pulses.end(), t.pulses.begin(), t.pulses.end());
        std::sort(tr.pulses.begin(), tr.pulses.end());
        take(std::move(tr));
    }
    return result;
}

const char* pri_type_name(PriType t) {
    switch (t) {
        case PriType::Constant: return "constant";
        case PriType::Stagger: return "stagger";
        case PriType::Jitter: return "jitter";
    }
    return "?";
}

}  // namespace esm
